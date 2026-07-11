/**
 * @file findInteractingBlocks.cu
 * @brief Block-based (Verlet) neighbor-list construction for OpenMM's CUDA nonbonded solver.
 *
 * ## What this file produces
 *
 * OpenMM does not build a per-atom neighbor list. Atoms are partitioned into fixed
 * groups of TILE_SIZE (=32) consecutive atoms called *atom blocks*. A pair of blocks
 * (X, Y) forms a *tile*; the nonbonded kernel evaluates a tile as a 32x32 warp-parallel
 * interaction matrix. The job of the four kernels here is to decide, cheaply and
 * conservatively, *which* tiles can contain any atom pair within PADDED_CUTOFF, so the
 * expensive nonbonded kernel only visits those. The output is consumed by
 * CudaNonbondedUtilities.cpp (which drives these kernels) and by the nonbonded force
 * kernel that reads interactingTiles / interactingAtoms / singlePairs.
 *
 * ## Pipeline (one host call = one neighbor-list rebuild; see prepareInteractions())
 *
 *   1. findBlockBounds          -> axis-aligned bounding box (center + half-widths) per
 *                                  atom block, plus per-thread-block min/max box size.
 *   2. computeSortKeys          -> a 32-bit sort key per block = (size bin << BIN_SHIFT)
 *                                  | blockIndex. Bins are logarithmic in box size.
 *   3. blockSorter->sort()      -> HOST-side radix sort of the keys (ComputeSort). Groups
 *                                  blocks of similar size/extent adjacently so the pair
 *                                  search below has coherent, well-culled comparisons.
 *   4. sortBoxData              -> gather block center/box into sorted order (and, for big
 *                                  systems, build 32-block "large block" super-boxes); also
 *                                  the displacement test that decides whether a rebuild is
 *                                  actually needed this step.
 *   5. findBlocksWithInteractions -> the actual O(N_blocks^2 / warps) tile search. Warp-
 *                                  cooperative box-box culling (stage 1) then atom-box
 *                                  refinement (stage 2); appends surviving tiles to the
 *                                  neighbor list and, optionally, sparse survivors to a
 *                                  single-pair list.
 *
 * This is a *conservative* neighbor list: it may keep tiles that turn out empty, but it
 * never drops a tile that contains an in-cutoff pair. Bounding boxes are rounded OUTWARD
 * (see half3) and PADDING (the Verlet skin) lets the list survive several MD steps before
 * atoms drift far enough to force a rebuild.
 *
 * ## Precision types
 * `real`, `real2`, `real3`, `real4` are float or double depending on the build's precision
 * (from CudaContext compilation defines). `trimTo3` (vectorOps.cu) drops the .w lane.
 *
 * ============================================================================
 * COMPILE-TIME PARAMETERS (injected as -D defines by NVRTC; comments are stripped
 * before JIT by openmm/cmake_modules/strip_comments.py, so this documentation is free).
 * Unless noted, these come from CudaNonbondedUtilities::createKernelsForGroups()
 * (openmm/platforms/cuda/src/CudaNonbondedUtilities.cpp, the `defines` map ~line 507-528).
 * ============================================================================
 *
 * @par TILE_SIZE
 *   int, value 32 (= CudaContext::TileSize, one warp). Atoms per block; also the tile
 *   edge length. A block's atoms are posq[x*TILE_SIZE .. x*TILE_SIZE+31].
 *
 * @par NUM_BLOCKS
 *   int = context.getNumAtomBlocks() = ceil(NUM_ATOMS / TILE_SIZE). Number of atom blocks.
 *
 * @par NUM_ATOMS
 *   int = context.getNumAtoms(). Real atom count. The last block may be padded; atom index
 *   NUM_ATOMS is used as the "no atom / padding" sentinel written into interactingAtoms.
 *
 * @par PADDING
 *   real = paddedCutoff - maxCutoff (nm). The Verlet skin: extra distance added to the
 *   physical cutoff so the list stays valid while atoms move. Used ONLY in sortBoxData's
 *   rebuild test, where an atom that has moved more than PADDING/2 (compared via the
 *   0.25*PADDING*PADDING squared threshold) triggers a rebuild. paddedCutoff = padCutoff(maxCutoff).
 *
 * @par PADDED_CUTOFF / PADDED_CUTOFF_SQUARED
 *   real = maxCutoff + skin, and its square (nm, nm^2). The distance the searches actually
 *   test against (physical cutoff + skin). All box-box and atom-box culls use this.
 *
 * @par BIN_SHIFT
 *   int = smallest b such that (1<<b) > NUM_BLOCKS. Number of LOW bits in a sort key
 *   reserved to hold a block index; computed in C++ by
 *     `int binShift=1; while (1<<binShift <= getNumAtomBlocks()) binShift++;`
 *   The size bin is stored ABOVE these bits (key = bin<<BIN_SHIFT | blockIndex), so a plain
 *   ascending sort orders primarily by size bin, then by block index.
 *
 * @par BLOCK_INDEX_MASK
 *   int = (1<<BIN_SHIFT) - 1. AND a sort key with this to recover the original block index
 *   (`sortedBlocks[i] & BLOCK_INDEX_MASK`). The complementary high bits are the size bin.
 *
 * @par MAX_EXCLUSIONS
 *   int = max over all blocks of the number of blocks it is fully excluded against. Sizes
 *   the per-warp shared `warpExclusions` scratch. If >32, an extra __syncthreads() guards
 *   the shared load. Block-pair (tile) level exclusions only; per-atom exclusion masks are
 *   applied later, inside the nonbonded kernel.
 *
 * @par MAX_BITS_FOR_PAIRS
 *   int, one of {0, 2, 3}. 0 when the single-pair optimization is disabled (!canUsePairList);
 *   otherwise 2 for compute capability < 8.0, else 3. A Y-atom whose interaction bitmask
 *   against block X has popcount <= MAX_BITS_FOR_PAIRS is emitted as individual (atomX,atomY)
 *   entries in singlePairs instead of consuming a whole tile (a "mostly empty tile" saver).
 *   All `#if MAX_BITS_FOR_PAIRS > 0` / saveSinglePairs code compiles out when 0.
 *
 * @par LOG
 *   Global CudaContext define (CudaContext.cpp): "log" (double precision) or "logf"
 *   (single). Natural log applied to a block's total half-extent (box.x+box.y+box.z) so
 *   computeSortKeys spreads sizes over 20 LOGARITHMIC bins (few huge boxes, many tiny ones).
 *
 * @par USE_PERIODIC
 *   Defined ("1") iff usePeriodic. Enables the APPLY_PERIODIC_* minimum-image wrapping below.
 *
 * @par TRICLINIC
 *   Defined iff the box is triclinic. Guards extra "detailed check" fallbacks where the
 *   nearest-image shortcut can miss a periodic copy more than half a box width away.
 *
 * @par USE_LARGE_BLOCKS
 *   Defined iff useLargeBlocks (context.getNumAtoms() > 90000). Enables a coarse pre-cull:
 *   sortBoxData builds a super-box over each run of 32 sorted blocks, and
 *   findBlocksWithInteractions skips 32 candidate blocks at once when their super-box is
 *   outside cutoff. A memory/latency win only for very large systems.
 *
 * @par GROUP_SIZE / BUFFER_SIZE
 *   Both 256, #defined literally at the top of THIS file (not injected). GROUP_SIZE is the
 *   thread-block size findBlocksWithInteractions is launched with (executeKernel(..., 256))
 *   and its __launch_bounds__; there are GROUP_SIZE/32 = 8 warps per block. BUFFER_SIZE is
 *   the per-warp capacity (in atoms) of the shared neighbor staging buffer.
 *
 * ## Macros supplied globally by CudaContext.cpp (not by this file's C++ driver)
 *   - BALLOT(v)                        -> __ballot_sync(0xffffffff, v): 32-bit warp vote.
 *   - APPLY_PERIODIC_TO_POS(pos)       -> wrap a position into the primary box.
 *   - APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center) -> wrap pos to the image nearest center.
 *   - APPLY_PERIODIC_TO_DELTA(delta)   -> minimum-image a displacement vector.
 *   (Rectangular vs. triclinic forms are chosen there by box type.)
 *
 * ## Runtime kernel arguments worth flagging (NOT compile-time)
 *   - startBlockIndex / numBlocks: sub-range of blocks this launch owns (multi-GPU / domain
 *     split via setAtomBlockRange). Normally 0 / NUM_BLOCKS.
 *   - maxTiles      = min(20*NUM_BLOCKS, totalTiles): capacity of interactingTiles/Atoms;
 *     overflow is detected on the host (pinnedCountBuffer) and the buffers are grown + rebuilt.
 *   - maxSinglePairs = 5*NUM_ATOMS: capacity of the singlePairs buffer.
 */
#define GROUP_SIZE 256
#define BUFFER_SIZE 256

/**
 * To use half precision, we're supposed to include cuda_fp16.h.  Unfortunately,
 * it isn't included in the search path automatically, and there's no reliable
 * way to find where it's located on disk.  Instead we provide our own definitions
 * for the few symbols we need.
 */
struct __align__(2) __half {
    unsigned short x;
};
__device__ __half __float2half_ru(const float f) {
    __half h;
    asm("{cvt.rp.f16.f32 %0, %1;}" : "=h"(*reinterpret_cast<unsigned short *>(&h)) : "f"(f));
    return h;
}
__device__ float __half2float(const __half h) {
    float f;
    asm("{cvt.f32.f16 %0, %1;}" : "=f"(f) : "h"(*reinterpret_cast<const unsigned short *>(&h)));
    return f;
}
/**
 * @brief Three IEEE-754 half-precision floats, used to store block bounding-box half-widths
 *        compactly (6 bytes vs 12/24) since bounding boxes are read very often in the search.
 *
 * The float->half conversion always rounds toward +inf (__float2half_ru). For box half-widths
 * that means every box is stored slightly LARGER than reality, so the conservative-culling
 * invariant holds: a rounded box never excludes a pair the exact box would have kept.
 */
struct half3 {
    __device__ half3(real3 f) {
        // Round up so we'll err on the side of making the box a little too large.
        // This ensures interactions will never be missed.
        v[0] = __float2half_ru((float) f.x);
        v[1] = __float2half_ru((float) f.y);
        v[2] = __float2half_ru((float) f.z);
    }
    __device__ real3 toReal3() const {
        return make_real3(__half2float(v[0]), __half2float(v[1]), __half2float(v[2]));
    }
private:
    __half v[3];
};

/**
 * @brief Pass 1 of the pipeline: compute an axis-aligned bounding box for the TILE_SIZE
 *        atoms of each atom block, and the range of box sizes seen.
 *
 * For each block, one thread scans its up-to-TILE_SIZE atoms, forming the min/max corner
 * (under minimum-image wrapping when periodic), then derives the box center and half-widths.
 * It also computes `center.w` = the max distance from any atom to the box center (a bounding
 * *sphere* radius), which later searches use to tighten the sphere-based pre-cull. Threads
 * then cooperatively reduce, within their thread block, the min and max of the boxes' total
 * half-extent (x+y+z); those extremes seed the logarithmic binning in computeSortKeys.
 *
 * @param numAtoms            Number of real atoms (= NUM_ATOMS). Loop bound; last block may be short.
 * @param periodicBoxSize     Rectangular box lengths (nm), .xyz. Used by APPLY_PERIODIC_* if periodic.
 * @param invPeriodicBoxSize  Reciprocals of the box lengths.
 * @param periodicBoxVecX/Y/Z Box edge vectors (triclinic); consumed by the periodic macros.
 * @param posq                [in] Atom positions+charge, real4 per atom, indexed by atom id.
 * @param blockCenter         [out] real4 per block: .xyz = box center (nm), .w = bounding-sphere
 *                            radius about that center. Indexed by ORIGINAL (unsorted) block id.
 * @param blockBoundingBox    [out] real4 per block: .xyz = box HALF-widths (nm), .w unused.
 *                            Original block order.
 * @param rebuildNeighborList [out] int[>=1]. Element 0 is cleared to 0 here (by block 0 / thread 0)
 *                            as the running "rebuild needed" flag that sortBoxData may later set to 1.
 * @param blockSizeRange      [out] real2 per THREAD BLOCK: (min, max) total half-extent (x+y+z)
 *                            over the blocks this thread block processed. Length = grid size.
 *
 * @return Writes blockCenter, blockBoundingBox (all blocks), blockSizeRange (one per thread block),
 *         and resets rebuildNeighborList[0].
 *
 * @pre Launched with ThreadBlockSize (=64) threads per block (executeKernel default): the shared
 *      minBuffer/maxBuffer are sized [64] and the reduction spans exactly 64 lanes, so the launch
 *      block size MUST be 64. Work size = getNumAtomBlocks(); grid is capped at numThreadBlocks and
 *      each thread strides over multiple blocks. blockSizeRange must have >= gridDim.x entries
 *      (numBlockSizes = min(ceil(NUM_BLOCKS/64), numThreadBlocks)).
 * @note center.w carries the bounding-sphere radius forward; note it is written into blockCenter,
 *       NOT blockBoundingBox. The min-size reduction ignores the 1e38 sentinel only implicitly
 *       (a thread that processed no block keeps minSize=1e38, which min() discards downstream).
 */
extern "C" __global__ void findBlockBounds(int numAtoms, real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        const real4* __restrict__ posq, real4* __restrict__ blockCenter, real4* __restrict__ blockBoundingBox, int* __restrict__ rebuildNeighborList,
        real2* __restrict__ blockSizeRange) {
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int base = index*TILE_SIZE;
    real minSize = 1e38, maxSize = 0;
    while (base < numAtoms) {
        real4 pos = posq[base];
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_POS(pos)
#endif
        real4 minPos = pos;
        real4 maxPos = pos;
        int last = min(base+TILE_SIZE, numAtoms);
        for (int i = base+1; i < last; i++) {
            pos = posq[i];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos, center)
#endif
            minPos = make_real4(min(minPos.x,pos.x), min(minPos.y,pos.y), min(minPos.z,pos.z), 0);
            maxPos = make_real4(max(maxPos.x,pos.x), max(maxPos.y,pos.y), max(maxPos.z,pos.z), 0);
        }
        real4 blockSize = 0.5f*(maxPos-minPos);
        real4 center = 0.5f*(maxPos+minPos);
        center.w = 0;
        for (int i = base; i < last; i++) {
            pos = posq[i];
            real4 delta = posq[i]-center;
#ifdef USE_PERIODIC
            APPLY_PERIODIC_TO_DELTA(delta)
#endif
            center.w = max(center.w, delta.x*delta.x+delta.y*delta.y+delta.z*delta.z);
        }
        center.w = sqrt(center.w);
        blockBoundingBox[index] = blockSize;
        blockCenter[index] = center;
        real totalSize = blockSize.x+blockSize.y+blockSize.z;
        minSize = min(minSize, totalSize);
        maxSize = max(maxSize, totalSize);
        index += blockDim.x*gridDim.x;
        base = index*TILE_SIZE;
    }
    
    // Record the range of sizes seen by threads in this block.

    __shared__ real minBuffer[64], maxBuffer[64];
    minBuffer[threadIdx.x] = minSize;
    maxBuffer[threadIdx.x] = maxSize;
    __syncthreads();
    for (int step = 1; step < 64; step *= 2) {
        if (threadIdx.x+step < 64 && threadIdx.x%(2*step) == 0) {
            minBuffer[threadIdx.x] = min(minBuffer[threadIdx.x], minBuffer[threadIdx.x+step]);
            maxBuffer[threadIdx.x] = max(maxBuffer[threadIdx.x], maxBuffer[threadIdx.x+step]);
        }
        __syncthreads();
    }
    if (threadIdx.x == 0)
        blockSizeRange[blockIdx.x] = make_real2(minBuffer[0], maxBuffer[0]);
    if (blockIdx.x == 0 && threadIdx.x == 0)
        rebuildNeighborList[0] = 0;
}

/**
 * @brief Pass 2: assign each block a 32-bit sort key encoding (size bin, block index).
 *
 * Thread 0 first reduces the per-thread-block extremes from blockSizeRange into a single
 * global (min, max) total half-extent, and takes their LOG (natural log). Every block's total
 * half-extent is then mapped through LOG into one of 20 equal-width bins across [LOGmin, LOGmax]
 * (logarithmic in linear size, so tiny and huge boxes both bin well). The key packs the bin in
 * the high bits and the original block index in the low BIN_SHIFT bits:
 *     sortedBlocks[i] = (bin << BIN_SHIFT) | i.
 * The subsequent host-side radix sort (blockSorter->sort, key = the whole word) therefore orders
 * blocks primarily by size, adjacently grouping similar-extent blocks for coherent culling.
 *
 * @param blockBoundingBox [in] real4 per block, .xyz = box half-widths (original order).
 * @param sortedBlocks     [out] unsigned int per block = (bin<<BIN_SHIFT)|blockIndex. Written in
 *                         original order; the host sort permutes it in place afterward.
 * @param blockSizeRange   [in] real2 per thread block, the (min,max) half-extents from findBlockBounds.
 * @param numSizes         Number of valid blockSizeRange entries (= numBlockSizes on the host).
 *
 * @return Writes sortedBlocks (unsorted keys) for all NUM_BLOCKS blocks.
 *
 * @pre A size-.x of exactly 0 marks an unused blockSizeRange slot and is skipped in the min. The
 *      grid strides over NUM_BLOCKS. sizeRange lives in shared memory; only thread 0 fills it, so
 *      the __syncthreads() before use is required.
 * @note numSizeBins is hard-coded to 20. `scale = 20/(LOGmax-LOGmin)`; the bin is clamped to [0,19].
 */
extern "C" __global__ void computeSortKeys(const real4* __restrict__ blockBoundingBox, unsigned int* __restrict__ sortedBlocks, real2* __restrict__ blockSizeRange, int numSizes) {
    // Find the total range of sizes recorded by all blocks.

    __shared__ real2 sizeRange;
    if (threadIdx.x == 0) {
        sizeRange = blockSizeRange[0];
        for (int i = 1; i < numSizes; i++) {
            real2 size = blockSizeRange[i];
            if (size.x > 0)
                sizeRange.x = min(sizeRange.x, size.x);
            sizeRange.y = max(sizeRange.y, size.y);
        }
        sizeRange.x = LOG(sizeRange.x);
        sizeRange.y = LOG(sizeRange.y);
    }
    __syncthreads();

    // Sort keys store the bin in the high order part and the block in the low
    // order part.

    int numSizeBins = 20;
    real scale = numSizeBins/(sizeRange.y-sizeRange.x);
    for (unsigned int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_BLOCKS; i += blockDim.x*gridDim.x) {
        real4 box = blockBoundingBox[i];
        // LOG is injected by CudaContext compilationDefines as log (double build) / logf (single);
        // binning on the log of extent keeps both a few huge boxes and many tiny ones well separated.
        real size = LOG(box.x+box.y+box.z);
        int bin = (size-sizeRange.x)*scale;
        bin = max(0, min(bin, numSizeBins-1));
        // BIN_SHIFT / BLOCK_INDEX_MASK are computed host-side in createKernelsForGroups: binShift is
        // the smallest b with (1<<b) > numAtomBlocks, so the low BIN_SHIFT bits hold the original block
        // index intact. The following host radix sort (blockSorter->sort) sorts the whole word ascending,
        // hence blocks come out ordered by size bin first, original index second.
        sortedBlocks[i] = (((unsigned int) bin)<<BIN_SHIFT) + i;
    }
}

/**
 * @brief Pass 4 (runs after the host radix sort): gather block geometry into sorted order and
 *        decide whether the neighbor list actually needs rebuilding this step.
 *
 * The keys in sortedBlocks are now sorted, but blockCenter/blockBoundingBox are still in original
 * order. For each sorted position i this scatters block (sortedBlocks[i] & BLOCK_INDEX_MASK)'s
 * center and (half-precision) box into sortedBlockCenter[i] / sortedBlockBoundingBox[i], so
 * findBlocksWithInteractions reads geometry sequentially without an indirection per access.
 *
 * When USE_LARGE_BLOCKS, it additionally builds, for each sorted position i, a "large block"
 * super-box enclosing the 32 sorted blocks starting at i (union of their AABBs, min-image aware).
 * findBlocksWithInteractions uses these to reject 32 candidate blocks with one test.
 *
 * Finally, every thread compares each atom's current position to oldPositions (the coordinates the
 * current list was built on). If any atom moved more than PADDING/2 (tested as squared displacement
 * > 0.25*PADDING*PADDING) — or forceRebuild is set — it raises rebuildNeighborList[0] and zeroes the
 * interaction counters, which is what actually authorizes findBlocksWithInteractions to run.
 *
 * @param sortedBlocks          [in] sorted 32-bit keys; low BIN_SHIFT bits = original block index.
 * @param blockCenter           [in] per-block center(.xyz)+sphere radius(.w), original order.
 * @param blockBoundingBox      [in] per-block half-widths, original order.
 * @param sortedBlockCenter     [out] blockCenter gathered into sorted order.
 * @param sortedBlockBoundingBox[out] half3 half-widths gathered into sorted order (rounded outward).
 * @param largeBlockCenter      [out, USE_LARGE_BLOCKS] super-box center per sorted position.
 * @param largeBlockBoundingBox [out, USE_LARGE_BLOCKS] super-box half-widths (half3) per sorted position.
 * @param periodicBoxSize..VecZ [in, USE_LARGE_BLOCKS] box params for min-image union of super-boxes.
 * @param posq                  [in] current atom positions (for the displacement test).
 * @param oldPositions          [in] atom positions the current neighbor list was built on.
 * @param interactionCount      [out] unsigned int[2]; both zeroed when a rebuild is triggered
 *                              ([0]=tile count, [1]=single-pair count).
 * @param rebuildNeighborList   [in/out] int[1]; set to 1 if a rebuild is needed.
 * @param forceRebuild          [in] bool; unconditionally forces a rebuild (e.g. box/param change).
 *
 * @return Writes sortedBlockCenter/BoundingBox (+ large-block arrays), and possibly
 *         rebuildNeighborList[0]=1 and interactionCount[0]=interactionCount[1]=0.
 *
 * @pre Launched over getNumAtoms() work items; the block-gather loop strides NUM_BLOCKS and the
 *      displacement loop strides NUM_ATOMS, so a single launch covers both (they are independent
 *      strided loops, not fused). No __syncthreads(): each thread's rebuild vote hits the same flag
 *      via a benign write race (any "true" wins).
 * @note The super-box run length 32 (blocks) matches the 32-lane BALLOT pre-cull in the next kernel.
 */
extern "C" __global__ void sortBoxData(const unsigned int* __restrict__ sortedBlocks, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockBoundingBox, real4* __restrict__ sortedBlockCenter, half3* __restrict__ sortedBlockBoundingBox,
#ifdef USE_LARGE_BLOCKS
        real4* __restrict__ largeBlockCenter, half3* __restrict__ largeBlockBoundingBox, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
#endif
        const real4* __restrict__ posq, const real4* __restrict__ oldPositions,
        unsigned int* __restrict__ interactionCount, int* __restrict__ rebuildNeighborList, bool forceRebuild) {
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_BLOCKS; i += blockDim.x*gridDim.x) {
        unsigned int index = sortedBlocks[i] & BLOCK_INDEX_MASK;
        sortedBlockCenter[i] = blockCenter[index];
        sortedBlockBoundingBox[i] = half3(trimTo3(blockBoundingBox[index]));

#ifdef USE_LARGE_BLOCKS
        // Compute the sizes of large blocks (composed of 32 regular blocks) starting from each block.
    
        real4 minPos = blockCenter[index]-blockBoundingBox[index];
        real4 maxPos = blockCenter[index]+blockBoundingBox[index];
        int last = min(i+32, NUM_BLOCKS);
        for (int j = i+1; j < last; j++) {
            unsigned int index2 = sortedBlocks[j] & BLOCK_INDEX_MASK;
            real4 blockPos = blockCenter[index2];
            real4 width = blockBoundingBox[index2];
#ifdef USE_PERIODIC
            real4 center = 0.5f*(maxPos+minPos);
            APPLY_PERIODIC_TO_POS_WITH_CENTER(blockPos, center)
#endif
            minPos = make_real4(min(minPos.x, blockPos.x-width.x), min(minPos.y, blockPos.y-width.y), min(minPos.z, blockPos.z-width.z), 0);
            maxPos = make_real4(max(maxPos.x, blockPos.x+width.x), max(maxPos.y, blockPos.y+width.y), max(maxPos.z, blockPos.z+width.z), 0);
        }
        largeBlockCenter[i] = 0.5f*(maxPos+minPos);
        largeBlockBoundingBox[i] = half3(trimTo3(0.5f*(maxPos-minPos)));
#endif
    }

    // Also check whether any atom has moved enough so that we really need to rebuild the neighbor list.

    bool rebuild = forceRebuild;
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x) {
        real4 delta = oldPositions[i]-posq[i];
        
        // PADDING (host: paddedCutoff-maxCutoff) is the Verlet skin. Threshold is (PADDING/2)^2:
        // once two atoms have each drifted skin/2 they could have closed by a full skin, enough to
        // bring a pair that was outside PADDED_CUTOFF at list-build time inside it now, so the list
        // may be stale. forceRebuild (host forceRebuildNeighborList) covers box/parameter changes.
        if (delta.x*delta.x + delta.y*delta.y + delta.z*delta.z > 0.25f*PADDING*PADDING)
            rebuild = true;
    }
    if (rebuild) {
        rebuildNeighborList[0] = 1;
        interactionCount[0] = 0;
        interactionCount[1] = 0;
    }
}

/**
 * @brief Warp helper (only compiled when MAX_BITS_FOR_PAIRS>0): drain "sparse" survivors from a
 *        warp's staging buffer into the single-pair list, and compact what remains.
 *
 * The staging buffer holds, for block X, a list of candidate Y-atoms (`atoms`) each with a 32-bit
 * `flags` mask of which of X's 32 atoms it interacts with. A Y-atom whose mask has popcount
 * <= MAX_BITS_FOR_PAIRS interacts with very few X atoms; materializing a whole 32x32 tile for it
 * would be mostly wasted work, so each set bit is emitted as an explicit (atomX, atomY) pair into
 * `singlePairs`. Entries with more bits are left in the buffer, compacted to the front, and the new
 * length returned so the caller can flush them as full tiles.
 *
 * Reservation into the global singlePairs array uses a warp-scan (shuffle prefix sum) of each lane's
 * pair count, one atomicAdd on singlePairCount by lane 31, then per-lane offsets — one atomic per
 * warp rather than per pair. Writes are dropped (not the count) if they would exceed maxSinglePairs;
 * the host detects the overflow via the returned count and grows the buffer.
 *
 * @param x               Original index of block X (its atoms are x*TILE_SIZE + bit).
 * @param atoms           [in/out] shared per-warp buffer of candidate Y atom indices; compacted in place.
 * @param flags           [in/out] shared per-warp buffer of 32-bit interaction masks, parallel to atoms.
 * @param length          Number of valid entries in atoms/flags.
 * @param maxSinglePairs  Capacity of singlePairs (= host maxSinglePairs); writes past it are skipped.
 * @param singlePairCount [in/out] global counter (interactionCount[1]) atomically advanced.
 * @param singlePairs     [out] global int2 list; .x = Y atom, .y = X atom.
 * @param sumBuffer       Per-warp scratch (aliases posBuffer); unused by the current scan path.
 * @param pairStartIndex  [in/out] per-warp shared slot holding this warp's reserved base offset.
 * @return New compacted `length`: the count of dense entries left for tile emission.
 * @note `#pragma unroll 8` == GROUP_SIZE/TILE_SIZE. The compaction is a standard BALLOT + warp
 *       prefix-count stream compaction keeping only popcount>MAX_BITS_FOR_PAIRS entries.
 */
__device__ int saveSinglePairs(int x, int* atoms, int* flags, int length, unsigned int maxSinglePairs, unsigned int* singlePairCount, int2* singlePairs, int* sumBuffer, volatile unsigned int& pairStartIndex) {
    // Record interactions that should be computed as single pairs rather than in blocks.
    // NOTE: sumBuffer (the caller passes sumBuffer+warpStart, which aliases posBuffer) is dead here —
    // the prefix sum below runs entirely in registers via __shfl_up_sync. Safe to drop from the
    // signature if the aliasing with posBuffer in the caller is dropped too.

    const int indexInWarp = threadIdx.x%32;
    int sum = 0;
    #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
    for (int i = indexInWarp; i < length; i += 32) {
        int count = __popc(flags[i]);
        sum += (count <= MAX_BITS_FOR_PAIRS ? count : 0);
    }
    for (int i = 1; i < 32; i *= 2) {
        int n = __shfl_up_sync(0xffffffff, sum, i);
        if (indexInWarp >= i)
            sum += n;
    }
    if (indexInWarp == 31)
        pairStartIndex = atomicAdd(singlePairCount,(unsigned int) sum);
    __syncwarp();
    int prevSum = __shfl_up_sync(0xffffffff, sum, 1);
    unsigned int pairIndex = pairStartIndex + (indexInWarp > 0 ? prevSum : 0);
    for (int i = indexInWarp; i < length; i += 32) {
        int count = __popc(flags[i]);
        if (count <= MAX_BITS_FOR_PAIRS && pairIndex+count <= maxSinglePairs) {
            int f = flags[i];
            while (f != 0) {
                singlePairs[pairIndex] = make_int2(atoms[i], x*TILE_SIZE+__ffs(f)-1);
                f &= f-1;
                pairIndex++;
            }
        }
    }
    
    // Compact the remaining interactions.
    
    const int warpMask = (1<<indexInWarp)-1;
    int numCompacted = 0;
    for (int start = 0; start < length; start += 32) {
        int i = start+indexInWarp;
        int atom = atoms[i];
        int flag = flags[i];
        bool include = (i < length && __popc(flags[i]) > MAX_BITS_FOR_PAIRS);
        int includeFlags = BALLOT(include);
        if (include) {
            int index = numCompacted+__popc(includeFlags&warpMask);
            atoms[index] = atom;
            flags[index] = flag;
        }
        numCompacted += __popc(includeFlags);
    }
    return numCompacted;
}

/**
 * @brief Pass 5 (the core search): find, for each atom block X, the tiles/atoms it interacts with,
 *        and append them to the neighbor list. One warp owns one block X; the 32 lanes probe 32
 *        candidate blocks Y in parallel. Two nested stages, coarse then fine.
 *
 * ### STAGE 1 — coarse block-vs-block cull
 * The warp fixes block X (all lanes same X) and streams candidate blocks Y = block2Base + lane over
 * the sorted blocks *after* X (each pair is thus visited once; the tile list is upper-triangular in
 * sorted order). A lane keeps its Y iff the X/Y bounding boxes are within PADDED_CUTOFF, using two
 * escalating tests: a cheap sphere test (center distance vs PADDED_CUTOFF + the two sphere radii
 * carried in center.w) then the exact AABB gap test. Y is dropped if X excludes Y (block-level
 * exclusion). Under USE_LARGE_BLOCKS a super-box test first skips 32 Ys at a stroke. BALLOT turns the
 * per-lane keep/drop into a mask the warp then walks bit by bit.
 *
 * ### STAGE 2 — fine atom-vs-block refinement
 * For each surviving Y, lane j takes Y's atom j and tests it against all 32 atoms of X (held in shared
 * posBuffer). The single-periodic-copy fast path uses the precomputed .w = 0.5*|r|^2 so an in-cutoff
 * test becomes a cheap dot product (halfDist2 = posj.w+pos2.w - dot < 0.5*PADDED_CUTOFF_SQUARED);
 * otherwise it min-images each delta. The 32-bit result `interacts` records which X atoms this Y atom
 * is near. Surviving (Y atom, mask) pairs accumulate in the warp's shared buffer; when it fills past
 * BUFFER_SIZE-TILE_SIZE the warp flushes: optionally siphoning sparse survivors to singlePairs
 * (saveSinglePairs, only if MAX_BITS_FOR_PAIRS>0), then emitting whole tiles via one atomicAdd on
 * interactionCount[0] and coalesced writes into interactingTiles/interactingAtoms. A partial buffer is
 * flushed after the Y loop, padding the tail atoms with the sentinel NUM_ATOMS.
 *
 * ### Neighbor-list encoding written here
 *   - interactionCount[0] = number of tiles; interactionCount[1] = number of single pairs.
 *   - interactingTiles[t]              = X (original block index) for tile t.
 *   - interactingAtoms[t*TILE_SIZE+k]  = the k-th neighbor Y atom of tile t (or NUM_ATOMS if padding).
 *     (The X atoms of a tile are implicit: X*TILE_SIZE .. +31.)
 *   - singlePairs[p] = int2(Y atom, X atom) for interactions too sparse to justify a tile.
 * A tile is thus X paired with an arbitrary set of up-to-32 neighbor atoms — NOT necessarily a single
 * block Y — since stage 2 packs whichever Y atoms passed, across multiple Y blocks, into 32-atom rows.
 *
 * @param periodicBoxSize        [in] rectangular box lengths (nm).
 * @param invPeriodicBoxSize     [in] reciprocal box lengths.
 * @param periodicBoxVecX/Y/Z    [in] box edge vectors (triclinic); used by APPLY_PERIODIC_*.
 * @param interactionCount       [out] unsigned int[2]: [0]=tile count, [1]=single-pair count (both
 *                               atomically advanced; zeroed earlier by sortBoxData on rebuild).
 * @param interactingTiles       [out] per tile, the X block index. Capacity maxTiles.
 * @param interactingAtoms       [out] per tile, TILE_SIZE neighbor atom ids (NUM_ATOMS = padding).
 * @param singlePairs            [out] int2 sparse-pair list (.x=Y atom, .y=X atom). Capacity maxSinglePairs.
 * @param posq                   [in] atom positions+charge.
 * @param maxTiles               [in] capacity of interactingTiles/Atoms; over-capacity tiles are counted
 *                               but not written (host detects overflow and grows+rebuilds).
 * @param maxSinglePairs         [in] capacity of singlePairs.
 * @param startBlockIndex        [in] first sorted block this launch owns (domain/multi-GPU split).
 * @param numBlocks              [in] number of sorted blocks this launch owns.
 * @param sortedBlocks           [in] sorted keys; low BIN_SHIFT bits = original block index.
 * @param sortedBlockCenter      [in] centers gathered into sorted order (.w = sphere radius).
 * @param sortedBlockBoundingBox [in] half3 half-widths in sorted order.
 * @param largeBlockCenter       [in, USE_LARGE_BLOCKS] super-box centers (per sorted position).
 * @param largeBlockBoundingBox  [in, USE_LARGE_BLOCKS] super-box half-widths.
 * @param exclusionIndices       [in] CSR *data*: flat, concatenated list of block indices excluded
 *                               against each block. (Despite the name, this is the value array.)
 * @param exclusionRowIndices    [in] CSR *row pointers*, length NUM_BLOCKS+1: block x's excluded
 *                               blocks are exclusionIndices[exclusionRowIndices[x] .. exclusionRowIndices[x+1]-1].
 *     Worked example (three blocks; block 0 excludes blocks {3,5,6}, block 1 excludes {3,4},
 *     block 2 excludes {1,3,5,6}):
 *         exclusionRowIndices = [0, 3, 5, 9]                 (row pointers, index by block)
 *         exclusionIndices    = [3,5,6, 3,4, 1,3,5,6]        (concatenated excluded-block ids)
 *     NOTE: the names are counter-intuitive — exclusionRowIndices is the offset/pointer array and
 *     exclusionIndices is the packed value array. (A prior version of this comment had them swapped.)
 * @param oldPositions           [out] the atom positions this neighbor list was built on; written at the
 *                               end so sortBoxData's next-step displacement test has a reference.
 * @param rebuildNeighborList    [in] int[1]; if 0 the kernel returns immediately (list still valid).
 *
 * @pre __launch_bounds__(GROUP_SIZE=256, 3): 256 threads/block, min 3 resident blocks/SM. This is
 *      the occupancy floor the register allocator must honor; combined with the STATIC shared-memory
 *      budget below it caps SM residency, so both are hard occupancy inputs for any retune. Launched
 *      with block size 256 (executeKernel(..., getNumAtoms(), 256)) => 8 warps/block; grid =
 *      min(ceil(NUM_ATOMS/256), numThreadBlocks), so warps stride over blocks (block1 loop). No
 *      dynamic shared memory is passed (third cuLaunchKernel arg = 0); all shared is static:
 *        - workgroupBuffer       BUFFER_SIZE*(GROUP_SIZE/32) ints        = 256*8*4  = 8192 B
 *        - workgroupFlagsBuffer  BUFFER_SIZE*(GROUP_SIZE/32) ints        = 256*8*4  = 8192 B
 *        - warpExclusions        MAX_EXCLUSIONS*(GROUP_SIZE/32) ints     = 32*MAX_EXCLUSIONS B
 *        - posBuffer             GROUP_SIZE real4 (sumBuffer aliases it) = 256*sizeof(real4)
 *                                                                        = 4096 B (single) / 8192 B (double)
 *        - workgroupTileIndex + workgroupPairStartIndex  (GROUP_SIZE/32) uints each = 64 B
 *      Total ~= 20544 + 32*MAX_EXCLUSIONS bytes (single precision); double precision adds 4096 B via
 *      posBuffer. Each warp owns the slice indexed by warpStart/32 of every per-warp array.
 * @note Correctness rests on outward-rounded boxes (half3) + PADDED_CUTOFF skin => conservative: a
 *       kept-but-empty tile only wastes work; an in-cutoff pair is never dropped. The triclinic
 *       branches force-include when the nearest image could be more than half a box width away.
 */
extern "C" __global__ __launch_bounds__(GROUP_SIZE,3) void findBlocksWithInteractions(real4 periodicBoxSize, real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ,
        unsigned int* __restrict__ interactionCount, int* __restrict__ interactingTiles, unsigned int* __restrict__ interactingAtoms,
        int2* __restrict__ singlePairs, const real4* __restrict__ posq, unsigned int maxTiles, unsigned int maxSinglePairs, unsigned int startBlockIndex,
        unsigned int numBlocks, unsigned int* __restrict__ sortedBlocks, const real4* __restrict__ sortedBlockCenter, const half3* __restrict__ sortedBlockBoundingBox,
#ifdef USE_LARGE_BLOCKS
        real4* __restrict__ largeBlockCenter, half3* __restrict__ largeBlockBoundingBox,
#endif
        const unsigned int* __restrict__ exclusionIndices, const unsigned int* __restrict__ exclusionRowIndices,
        real4* __restrict__ oldPositions, const int* __restrict__ rebuildNeighborList) {

    // rebuildNeighborList[0] is the cross-kernel handshake: cleared to 0 by findBlockBounds every pass
    // and raised to 1 only by sortBoxData (or host forceRebuild). All four kernels launch unconditionally
    // each step (see prepareInteractions); this guard is what makes the search a no-op when the skin still
    // holds, so the frequent case here is an early return, not a full O(N_blocks^2) search.
    if (rebuildNeighborList[0] == 0)
        return; // The neighbor list doesn't need to be rebuilt.

    const int indexInWarp = threadIdx.x%32;
    const int warpStart = threadIdx.x-indexInWarp;
    const int totalWarps = blockDim.x*gridDim.x/32;
    const int warpIndex = (blockIdx.x*blockDim.x+threadIdx.x)/32;
    const int warpMask = (1<<indexInWarp)-1;
    // These static allocations (~20 KB single / ~24 KB double per block, plus 32*MAX_EXCLUSIONS) are
    // the binding constraint behind __launch_bounds__(256,3): the ",3" promises >=3 resident blocks/SM,
    // so on a 48 KB-shared SM this budget is what actually gates occupancy — retuning either number
    // must be checked against the other. No dynamic shared memory is requested at the launch site.
    __shared__ int workgroupBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int workgroupFlagsBuffer[BUFFER_SIZE*(GROUP_SIZE/32)];
    __shared__ int warpExclusions[MAX_EXCLUSIONS*(GROUP_SIZE/32)];
    __shared__ real4 posBuffer[GROUP_SIZE];
    __shared__ volatile unsigned int workgroupTileIndex[GROUP_SIZE/32];
    __shared__ unsigned int workgroupPairStartIndex[GROUP_SIZE/32];
    // Per-warp views into the shared scratch (each warp owns slice warpStart/32).
    int* sumBuffer = (int*) posBuffer; // Reuse the same buffer to save memory
    int* buffer = workgroupBuffer+BUFFER_SIZE*(warpStart/32);
    int* flagsBuffer = workgroupFlagsBuffer+BUFFER_SIZE*(warpStart/32);
    int* exclusionsForX = warpExclusions+MAX_EXCLUSIONS*(warpStart/32);
    volatile unsigned int& tileStartIndex = workgroupTileIndex[warpStart/32];
    volatile unsigned int& pairStartIndex = workgroupPairStartIndex[warpStart/32];

    // Loop over blocks.
    
    for (int block1 = startBlockIndex+warpIndex; block1 < startBlockIndex+numBlocks; block1 += totalWarps) {
        // Load data for this block.  Note that all threads in a warp are processing the same block.
        
        int x = sortedBlocks[block1] & BLOCK_INDEX_MASK;
        real4 blockCenterX = sortedBlockCenter[block1];
        real3 blockSizeX = sortedBlockBoundingBox[block1].toReal3();
        int neighborsInBuffer = 0;
        real4 pos1 = posq[x*TILE_SIZE+indexInWarp];
#ifdef USE_PERIODIC
        const bool singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.y-blockSizeX.y >= PADDED_CUTOFF &&
                                         0.5f*periodicBoxSize.z-blockSizeX.z >= PADDED_CUTOFF);
        if (singlePeriodicCopy) {
            // The box is small enough that we can just translate all the atoms into a single periodic
            // box, then skip having to apply periodic boundary conditions later.
            
            APPLY_PERIODIC_TO_POS_WITH_CENTER(pos1, blockCenterX)
        }
#endif
        // Stash 0.5*|pos1|^2 in .w so the single-copy fast path can test distance with a dot
        // product only: 0.5*|a-b|^2 = a.w + b.w - dot(a,b) < 0.5*PADDED_CUTOFF_SQUARED.
        pos1.w = 0.5f * (pos1.x * pos1.x + pos1.y * pos1.y + pos1.z * pos1.z);
        posBuffer[threadIdx.x] = pos1;

        // Load exclusion data for block x.
        // CSR, and the parameter names are the reverse of what they suggest: exclusionRowIndices is the
        // row-pointer array (host exclusionRowIndicesVec, length numAtomBlocks+1) and exclusionIndices is
        // the concatenated value list (host exclusionIndicesVec). So block x's fully-excluded partner
        // blocks are exclusionIndices[exclusionRowIndices[x] .. exclusionRowIndices[x+1]-1].
        const int exclusionStart = exclusionRowIndices[x];
        const int exclusionEnd = exclusionRowIndices[x+1];
        const int numExclusions = exclusionEnd-exclusionStart;
        #pragma unroll 4 // (MAX_EXCLUSIONS)
        for (int j = indexInWarp; j < numExclusions; j += 32)
            exclusionsForX[j] = exclusionIndices[exclusionStart+j];
        if (MAX_EXCLUSIONS > 32)
            __syncthreads();
        
        // Loop over atom blocks to search for neighbors.  The threads in a warp compare block1 against 32
        // other blocks in parallel.

#ifdef USE_LARGE_BLOCKS
        int largeBlockFlags = 0;
        int loadedLargeBlocks = 0;
#endif
        for (int block2Base = block1+1; block2Base < NUM_BLOCKS; block2Base += 32) {
#ifdef USE_LARGE_BLOCKS
            if (loadedLargeBlocks == 0) {
                // Check the next set of large blocks.

                int largeBlockIndex = block2Base + 32*indexInWarp;
                bool includeLargeBlock = false;
                if (largeBlockIndex < NUM_BLOCKS) {
                    real4 largeCenter = largeBlockCenter[largeBlockIndex];
                    real3 largeSize = largeBlockBoundingBox[largeBlockIndex].toReal3();
                    real4 blockDelta = blockCenterX-largeCenter;
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                    blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-largeSize.x);
                    blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-largeSize.y);
                    blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-largeSize.z);
                    includeLargeBlock = (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                    // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                    // If there's any possibility we might have missed it, do a detailed check.

                    if (periodicBoxSize.z/2-blockSizeX.z-largeSize.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-largeSize.y < PADDED_CUTOFF)
                        includeLargeBlock = true;
#endif
                }
                largeBlockFlags = BALLOT(includeLargeBlock);
                loadedLargeBlocks = 32;
            }
            loadedLargeBlocks--;
            if ((largeBlockFlags&1) == 0) {
                // None of the next 32 blocks interact with block 1.

                largeBlockFlags >>= 1;
                continue;
            }
            largeBlockFlags >>= 1;
#endif
            int block2 = block2Base+indexInWarp;
            bool includeBlock2 = (block2 < NUM_BLOCKS);
            bool forceInclude = false;
            if (includeBlock2) {
                real4 blockCenterY = sortedBlockCenter[block2];
                real3 blockSizeY = sortedBlockBoundingBox[block2].toReal3();
                real4 blockDelta = blockCenterX-blockCenterY;
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(blockDelta)
#endif
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < (PADDED_CUTOFF+blockCenterX.w+blockCenterY.w)*(PADDED_CUTOFF+blockCenterX.w+blockCenterY.w));
                blockDelta.x = max(0.0f, fabs(blockDelta.x)-blockSizeX.x-blockSizeY.x);
                blockDelta.y = max(0.0f, fabs(blockDelta.y)-blockSizeX.y-blockSizeY.y);
                blockDelta.z = max(0.0f, fabs(blockDelta.z)-blockSizeX.z-blockSizeY.z);
                includeBlock2 &= (blockDelta.x*blockDelta.x+blockDelta.y*blockDelta.y+blockDelta.z*blockDelta.z < PADDED_CUTOFF_SQUARED);
#ifdef TRICLINIC
                // The calculation to find the nearest periodic copy is only guaranteed to work if the nearest copy is less than half a box width away.
                // If there's any possibility we might have missed it, do a detailed check.

                if (periodicBoxSize.z/2-blockSizeX.z-blockSizeY.z < PADDED_CUTOFF || periodicBoxSize.y/2-blockSizeX.y-blockSizeY.y < PADDED_CUTOFF)
                    includeBlock2 = forceInclude = true;
#endif
                if (includeBlock2) {
                    int y = sortedBlocks[block2] & BLOCK_INDEX_MASK;
                    #pragma unroll 4 // (MAX_EXCLUSIONS)
                    for (int k = 0; k < numExclusions; k++)
                        includeBlock2 &= (exclusionsForX[k] != y);
                }
            }
            
            // Loop over any blocks we identified as potentially containing neighbors.
            
            int includeBlockFlags = BALLOT(includeBlock2);
            int forceIncludeFlags = BALLOT(forceInclude);
            while (includeBlockFlags != 0) {
                int i = __ffs(includeBlockFlags)-1;
                includeBlockFlags &= includeBlockFlags-1;
                forceInclude = (forceIncludeFlags>>i) & 1;
                int y = sortedBlocks[block2Base+i] & BLOCK_INDEX_MASK;

                // Check each atom in block Y for interactions.

                int atom2 = y*TILE_SIZE+indexInWarp;
                real4 pos2 = posq[atom2];
#ifdef USE_PERIODIC
                if (singlePeriodicCopy) {
                    APPLY_PERIODIC_TO_POS_WITH_CENTER(pos2, blockCenterX)
                }
#endif
                pos2.w = 0.5f * (pos2.x * pos2.x + pos2.y * pos2.y + pos2.z * pos2.z);

                real4 blockCenterY = sortedBlockCenter[block2Base+i];
                real3 atomDelta = trimTo3(posBuffer[warpStart+indexInWarp])-trimTo3(blockCenterY);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(atomDelta)
#endif
                int atomFlags = BALLOT(forceInclude || atomDelta.x*atomDelta.x+atomDelta.y*atomDelta.y+atomDelta.z*atomDelta.z < (PADDED_CUTOFF+blockCenterY.w)*(PADDED_CUTOFF+blockCenterY.w));
                int interacts = 0;
                if (atom2 < NUM_ATOMS && atomFlags != 0) {
#ifdef USE_PERIODIC
                    if (!singlePeriodicCopy) {
                        int first = __ffs(atomFlags)-1;
                        int last = 32-__clz(atomFlags);
                        for (int j = first; j < last; j++) {
                            real3 delta = trimTo3(pos2)-trimTo3(posBuffer[warpStart+j]);
                            APPLY_PERIODIC_TO_DELTA(delta)
                            interacts |= (delta.x*delta.x+delta.y*delta.y+delta.z*delta.z < PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                        }
                    }
                    else {
#endif
                        #pragma unroll
                        for (int j = 0; j < 32; j++) {
                            real4 posj = posBuffer[warpStart+j];
                            real halfDist2 = posj.w + pos2.w - posj.x*pos2.x - posj.y*pos2.y - posj.z*pos2.z;
                            interacts |= (halfDist2 < 0.5f * PADDED_CUTOFF_SQUARED ? 1<<j : 0);
                        }
#ifdef USE_PERIODIC
                    }
#endif
                }
                
                // Add any interacting atoms to the buffer.
                
                int includeAtomFlags = BALLOT(interacts);
                if (interacts) {
                    int index = neighborsInBuffer+__popc(includeAtomFlags&warpMask);
                    buffer[index] = atom2;
                    flagsBuffer[index] = interacts;
                }
                neighborsInBuffer += __popc(includeAtomFlags);
                if (neighborsInBuffer > BUFFER_SIZE-TILE_SIZE) {
                    // Store the new tiles to memory.
                    
#if MAX_BITS_FOR_PAIRS > 0
                    neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs, sumBuffer+warpStart, pairStartIndex);
#endif
                    unsigned int tilesToStore = neighborsInBuffer/TILE_SIZE;
                    if (tilesToStore > 0) {
                        if (indexInWarp == 0)
                            tileStartIndex = atomicAdd(&interactionCount[0], tilesToStore);
                        unsigned int newTileStartIndex = tileStartIndex;
                        // interactionCount[0] keeps counting past maxTiles but writes are suppressed here;
                        // the host reads interactionCount via pinnedCountBuffer, and when it exceeds
                        // maxTiles (init 20*numAtomBlocks) it grows the arrays to 1.2x and forces a full
                        // rebuild, so a truncated list this pass is discarded rather than consumed.
                        if (newTileStartIndex+tilesToStore <= maxTiles) {
                            if (indexInWarp < tilesToStore)
                                interactingTiles[newTileStartIndex+indexInWarp] = x;
                            #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                            for (int j = 0; j < tilesToStore; j++)
                                interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = buffer[indexInWarp+j*TILE_SIZE];
                        }
                        if (indexInWarp+TILE_SIZE*tilesToStore < BUFFER_SIZE)
                            buffer[indexInWarp] = buffer[indexInWarp+TILE_SIZE*tilesToStore];
                        neighborsInBuffer -= TILE_SIZE*tilesToStore;
                    }
                }
            }
        }
        
        // If we have a partially filled buffer,  store it to memory.
        
#if MAX_BITS_FOR_PAIRS > 0
        if (neighborsInBuffer > 32)
            neighborsInBuffer = saveSinglePairs(x, buffer, flagsBuffer, neighborsInBuffer, maxSinglePairs, &interactionCount[1], singlePairs, sumBuffer+warpStart, pairStartIndex);
#endif
        if (neighborsInBuffer > 0) {
            unsigned int tilesToStore = (neighborsInBuffer+TILE_SIZE-1)/TILE_SIZE;
            if (indexInWarp == 0)
                tileStartIndex = atomicAdd(&interactionCount[0], tilesToStore);
            unsigned int newTileStartIndex = tileStartIndex;
            if (newTileStartIndex+tilesToStore <= maxTiles) {
                if (indexInWarp < tilesToStore)
                    interactingTiles[newTileStartIndex+indexInWarp] = x;
                #pragma unroll 8 // (GROUP_SIZE / TILE_SIZE)
                for (int j = 0; j < tilesToStore; j++)
                    interactingAtoms[(newTileStartIndex+j)*TILE_SIZE+indexInWarp] = (indexInWarp+j*TILE_SIZE < neighborsInBuffer ? buffer[indexInWarp+j*TILE_SIZE] : NUM_ATOMS);
            }
        }
    }
    
    // Record the positions the neighbor list is based on.
    
    for (int i = threadIdx.x+blockIdx.x*blockDim.x; i < NUM_ATOMS; i += blockDim.x*gridDim.x)
        oldPositions[i] = posq[i];
}
