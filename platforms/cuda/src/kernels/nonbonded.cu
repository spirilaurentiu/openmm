/**
 * @file nonbonded.cu
 * @brief Tile-based O(N^2) nonbonded interaction kernel (`computeNonbonded`) for the OpenMM CUDA
 *        platform, vendored into Robosample from OpenMM 8.5.
 *
 * This translation unit is a *template*, not a standalone kernel. It is JIT-compiled once per
 * force group by `CudaNonbondedUtilities::createInteractionKernel` (in
 * `../CudaNonbondedUtilities.cpp`). Two independent substitution mechanisms fill it in:
 *
 *   1. String replacement (`context.replaceStrings(kernelSource, replacements)`): the ALL-CAPS
 *      CODE-SNIPPET macros below are replaced textually with C++/CUDA source that
 *      `createInteractionKernel` assembles per force from that force's parameter list and its
 *      Lepton-derived energy/force expression, BEFORE the source is handed to NVRTC.
 *   2. Preprocessor `#define`s (the `defines` map): the SCALAR VALUES and `USE_*` feature flags
 *      below become ordinary `-D` macros seen by NVRTC.
 *
 * Additional macros (`real`, `real3`, `real4`, `mixed`, `tileflags`, `make_real3`, `RSQRT`,
 * `SQRT`, `SHFL`, `APPLY_PERIODIC_TO_DELTA`, `APPLY_PERIODIC_TO_POS_WITH_CENTER`) come from the
 * common preamble that `CudaContext` prepends (`CudaContext.cpp` compilationDefines) and from
 * `CudaKernelSources::vectorOps`; `realToFixedPoint` is defined in `common.cu`. See the NOTE at
 * the bottom of this block.
 *
 * IMPORTANT: comments here are stripped by `openmm/cmake_modules/strip_comments.py` before NVRTC
 * sees the source, so this documentation carries ZERO runtime cost.
 *
 * ============================================================================================
 * (1) COMPILE-TIME PARAMETERS -- SCALAR VALUES  (plain `#define`s injected as NVRTC macros)
 * ============================================================================================
 * All are set in `createInteractionKernel`'s `defines` map (or `createKernelsForGroups` /
 * `CudaContext`) unless noted. "int literal" means the value is a decimal integer literal.
 *
 *   PADDED_NUM_ATOMS   int. context.getPaddedNumAtoms(). The SoA row stride: force/pos buffers
 *                      have three planes at [atom], [atom+PADDED_NUM_ATOMS], [atom+2*PADDED_NUM_ATOMS].
 *                      Padded up to a multiple of TILE_SIZE so partial tiles are always full width.
 *   NUM_ATOMS          int. context.getNumAtoms(). The real atom count; indices in
 *                      [NUM_ATOMS, PADDED_NUM_ATOMS) are padding and are force-masked / excluded.
 *   NUM_BLOCKS         int. context.getNumAtomBlocks() = ceil(NUM_ATOMS/TILE_SIZE). Number of
 *                      atom blocks; used to enumerate the lower-triangular tile grid.
 *   TILE_SIZE          int. CudaContext::TileSize (== 32, the warp size). A tile is a
 *                      TILE_SIZE x TILE_SIZE block-vs-block atom pair set; one warp owns one tile.
 *   THREAD_BLOCK_SIZE  int. forceThreadBlockSize (128 on CC<2.0 else 256). CUDA block width;
 *                      WARPS_PER_GROUP = THREAD_BLOCK_SIZE/TILE_SIZE warps per block.
 *   NUM_TILES_WITH_EXCLUSIONS  int. exclusionTiles.getSize(); length of the exclusionTiles list.
 *   FIRST_EXCLUSION_TILE / LAST_EXCLUSION_TILE  int. This context's [start,end) slice of the
 *                      exclusion-tile list (multi-GPU partition by context index). The first loop
 *                      walks [FIRST_EXCLUSION_TILE, LAST_EXCLUSION_TILE).
 *   MAX_CUTOFF         real (double literal). Largest cutoff over the force groups in this kernel;
 *                      used for the singlePeriodicCopy fast-path test.
 *   CUTOFF_<g>, CUTOFF_<g>_SQUARED  real. Per-force-group cutoff and its square; the force's
 *                      COMPUTE_INTERACTION references CUTOFF / CUTOFF_SQUARED, which addInteraction
 *                      rewrote to the group-specific names.
 *   PARAMETER_SIZE_IS_EVEN  Defined (=1) iff (localDataSize/4) is even and single precision. Not
 *                      referenced by this file; consumed by other shared kernels.
 *
 *   Feature flags (defined => 1, else undefined; guard `#ifdef`/`#if`):
 *     USE_CUTOFF        useCutoff. Adds the neighbor-list argument block to the signature.
 *     USE_PERIODIC      usePeriodic. Enables minimum-image via APPLY_PERIODIC_TO_DELTA / center wrap.
 *     USE_EXCLUSIONS    useExclusions (true here). Enables the 1024-bit per-tile exclusion mask path.
 *     USE_SYMMETRIC     isSymmetric (true here). Force is central: dEdR scalar along `delta`, with
 *                       Newton's third law giving atom2 the negation. When undefined, the snippet
 *                       instead fills independent dEdR1/dEdR2 vectors (non-central forces).
 *     USE_NEIGHBOR_LIST useNeighborList. Second loop reads the interacting-tile list instead of
 *                       enumerating the full lower triangle; enables the third (single-pair) loop.
 *     INCLUDE_FORCES / INCLUDE_ENERGY  whether this JIT variant accumulates forces / energy
 *                       (createInteractionKernel is called separately for force / energy / both).
 *     ENABLE_SHUFFLE    always 1 (defines map in createInteractionKernel). The kernel has no
 *                       non-shuffle fallback.
 *
 *     NOT defined for THIS module: USE_LARGE_BLOCKS and TRICLINIC are added by
 *     createKernelsForGroups only to the SEPARATE findInteractingBlocks module (neighbor-list
 *     build), not to computeNonbonded. This kernel's rectangular-vs-triclinic behavior is instead
 *     baked into the APPLY_PERIODIC_* macros by CudaContext (compilationDefines, selected by
 *     boxIsTriclinic), so it never branches on TRICLINIC directly.
 *
 * ============================================================================================
 * (2) COMPILE-TIME PARAMETERS -- CODE-SNIPPET MACROS  (textual `replaceStrings` substitution)
 * ============================================================================================
 * Each expands to a BLOCK of generated CUDA source. Exact content is FORCE-DEPENDENT: it is built
 * by looping over this force's `params` (each a ComputeParameterInfo: name, type, #components) and,
 * for COMPUTE_INTERACTION, from the force's Lepton energy expression. INVARIANT (established at the
 * createInteractionKernel dispatch boundary, unverifiable by the compiler after substitution): all
 * of these snippets are generated from the SAME parameter set in the SAME order, so the register
 * names (P1, P2, shflP) and the global_<param> formals agree in name, type, and component count
 * across LOAD_ATOM1/2, DECLARE/CLEAR/LOAD_LOCAL, SHUFFLE/BROADCAST, and COMPUTE_INTERACTION. Mixing
 * snippets from different forces is undefined behavior. The templates (from createInteractionKernel)
 * are:
 *
 *   PARAMETER_ARGUMENTS  Trailing kernel formal parameters. For each per-atom param P:
 *                        ", const T* __restrict__ global_P" (SoA, one array per param, indexed by
 *                        atom). Then any extra `arguments`, then optionally
 *                        ", mixed* __restrict__ energyParamDerivs". Appended to the signature.
 *   INIT_DERIVATIVES     One "mixed energyParamDerivN = 0;" per energy-parameter derivative (empty
 *                        if none). Declares per-thread derivative accumulators.
 *   LOAD_ATOM1_PARAMETERS  For each param P: "T P1 = global_P[atom1];". Loads the row-atom's params
 *                        into registers named P1 (referenced by COMPUTE_INTERACTION).
 *   DECLARE_LOCAL_PARAMETERS  For each P: "T shflP;". Register copies of atom-j params that get
 *                        rotated around the warp by SHUFFLE_WARP_DATA.
 *   LOAD_LOCAL_PARAMETERS_FROM_GLOBAL  For each P: "shflP = global_P[j];". Fills the shfl* copies
 *                        from global memory for the off-diagonal / neighbor-list column atom.
 *   CLEAR_LOCAL_PARAMETERS  For each P: "shflP = 0;" (or make_T(0)). Zeroes shfl* for padding atoms
 *                        (j >= PADDED_NUM_ATOMS).
 *   LOAD_ATOM2_PARAMETERS  For each P: "T P2 = shflP;". Materializes atom-j params as P2 for the
 *                        current warp lane (referenced by COMPUTE_INTERACTION).
 *   LOAD_ATOM2_PARAMETERS_FROM_GLOBAL  For each P: "T P2 = global_P[atom2];". Direct load, used in
 *                        the single-pair loop where there is no warp rotation.
 *   BROADCAST_WARP_DATA  Diagonal-tile fast path. Emits "posq2.{x,y,z,w} = real_shfl(shflPosq.*, j);"
 *                        and, per P, "T shflP; shflP = real_shfl(P1, j);" -- i.e. broadcasts lane j's
 *                        atom1 data to the whole warp so no second global load is needed on the
 *                        symmetric diagonal tile.
 *   SHUFFLE_WARP_DATA    Off-diagonal / neighbor rotation. Emits real_shfl(..., tgx+1) for
 *                        shflPosq.{x,y,z,w}, shflForce.{x,y,z} and every shflP, advancing atom-j
 *                        data (and its accumulated reaction force) by one lane each inner iteration.
 *   COMPUTE_INTERACTION  THE FORCE ITSELF. This is `source` = the force's kernel snippet (the Lepton
 *                        energy expression compiled to CUDA by CudaExpressionUtilities, concatenated
 *                        across all forces in the group). It reads posq1/posq2, r, r2, invR, the P1/P2
 *                        params, isExcluded, hasExclusions, interactionScale, and writes tempEnergy
 *                        and either dEdR (USE_SYMMETRIC) or dEdR1/dEdR2. Content varies entirely by
 *                        force; described here as a contract, not a fixed instance.
 *   SAVE_DERIVATIVES     Per energy-parameter derivative: "energyParamDerivs[GLOBAL_ID*numDerivs+idx]
 *                        += energyParamDerivN;". Flushes the accumulators at kernel end.
 *
 * NOTE: the periodic macros expand differently for rectangular vs. triclinic boxes (see
 * CudaContext.cpp): rectangular does per-axis minimum image; triclinic subtracts box vectors in
 * z,y,x order. APPLY_PERIODIC_TO_DELTA rounds (…+0.5) to nearest image; the WITH_CENTER form wraps a
 * position into the image nearest a block center. SHFL expands to __shfl_sync(0xffffffff, …) on
 * modern archs. realToFixedPoint scales by 0x100000000 (2^32) -- see the force-accumulation NOTE on
 * computeNonbonded.
 */
#define WARPS_PER_GROUP (THREAD_BLOCK_SIZE/TILE_SIZE)

/**
 * @brief Warp shuffle helpers: return the value of `var` held by lane `srcLane` of the caller's warp.
 *
 * Four overloads select on the element width. The 32-bit float/int overloads forward one SHFL; the
 * 64-bit double/long long overloads split the value into two 32-bit halves (the hardware SHFL moves
 * 32 bits at a time), shuffle each half, and recombine -- hence "support for 64 bit shuffles". These
 * are the sole intra-warp transport in computeNonbonded (broadcast on diagonal tiles, one-lane
 * rotation on off-diagonal / neighbor tiles) and touch neither shared nor global memory.
 *
 * @par Participation (contract): SHFL expands to __shfl_sync with a full 0xffffffff mask, so EVERY
 *      lane of the warp must call this in lockstep; a partially-diverged warp yields undefined lane
 *      values. `srcLane` is reduced modulo warpSize by the intrinsic (the rotation callers pass
 *      tgx+1, which wraps at lane 31). For 64-bit values both halves must reach their SHFL, which
 *      the full-mask, non-divergent call site guarantees.
 *
 * @param var     [in] the per-lane value to be read from `srcLane`.
 * @param srcLane [in] source lane index within the warp (taken mod warpSize).
 * @return The value of `var` on lane `srcLane`.
 * @note The long long overload recombines via an int2 bit-cast because there is no
 *       __hiloint2longlong intrinsic; it is a bit reinterpretation, not an arithmetic conversion.
 */
//support for 64 bit shuffles
static __inline__ __device__ float real_shfl(float var, int srcLane) {
    return SHFL(var, srcLane);
}

static __inline__ __device__ float real_shfl(int var, int srcLane) {
    return SHFL(var, srcLane);
}

static __inline__ __device__ double real_shfl(double var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "d"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    return __hiloint2double( hi, lo );
}

static __inline__ __device__ long long real_shfl(long long var, int srcLane) {
    int hi, lo;
    asm volatile("mov.b64 { %0, %1 }, %2;" : "=r"(lo), "=r"(hi) : "l"(var));
    hi = SHFL(hi, srcLane);
    lo = SHFL(lo, srcLane);
    // unforunately there isn't an __nv_hiloint2long(hi,lo) intrinsic cast
    int2 fuse; fuse.x = lo; fuse.y = hi;
    return *reinterpret_cast<long long*>(&fuse);
}

/**
 * @brief Atomically accumulate a real-valued force vector onto one atom's fixed-point force buffer.
 *
 * Each nonzero Cartesian component is converted to 64-bit signed fixed-point (realToFixedPoint,
 * scale 2^32) and added with atomicAdd. Because integer atomicAdd is associative and commutative,
 * the resulting force sum is ORDER-INDEPENDENT and BIT-EXACT REPRODUCIBLE across runs irrespective
 * of how threads race -- this is the determinism guarantee of the whole force pipeline, not an
 * incidental detail. Exactly-zero components are skipped to elide the atomic; this changes nothing
 * numerically (adding 0 is identity). Called per thread by the phase-3 single-pair loop of
 * computeNonbonded; imposes no warp-participation or synchronization requirement of its own.
 *
 * @param atom          [in] global atom index; row into each SoA force plane. Must be < the buffer
 *                      capacity (callers pass real atom indices from `singlePairs`, always valid).
 * @param force         [in] force to add for `atom`, in the engine's real force units.
 * @param forceBuffers  [inout] device SoA fixed-point accumulator; the x/y/z planes live at
 *                      [atom], [atom+PADDED_NUM_ATOMS], [atom+2*PADDED_NUM_ATOMS]. Read-modify-write
 *                      via atomicAdd; borrowed (not owned here).
 */
__device__ void saveSingleForce(int atom, real3 force, unsigned long long* forceBuffers) {
    if (force.x != 0)
        atomicAdd(&forceBuffers[atom], static_cast<unsigned long long>(realToFixedPoint(force.x)));
    if (force.y != 0)
        atomicAdd(&forceBuffers[atom+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
    if (force.z != 0)
        atomicAdd(&forceBuffers[atom+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
}

/**
 * @brief Compute all pairwise nonbonded interactions (energy, forces, energy-parameter
 *        derivatives) for one JIT-specialized force group via a tile decomposition of the
 *        atom-vs-atom interaction matrix.
 *
 * This is the OpenMM CUDA platform's central force-evaluation hotspot. It is JIT-built per force
 * group by CudaNonbondedUtilities::createInteractionKernel and launched (not via `<<< >>>`) by
 * CudaNonbondedUtilities::computeInteractions through CudaContext::executeKernel -> cuLaunchKernel.
 * Everything below is the CONTRACT an optimizer or profiler may rely on; it survives any
 * behavior-preserving rewrite. The functional form of each pair interaction is not part of this
 * file -- it is the injected COMPUTE_INTERACTION snippet (see the file-header section 2).
 *
 * @par Decomposition and warp-per-tile mapping (contract):
 * Atoms are grouped into NUM_BLOCKS blocks of TILE_SIZE (== warp size, 32) atoms. A tile is the
 * interaction of block x against block y -- a TILE_SIZE x TILE_SIZE sub-matrix. ONE WARP OWNS ONE
 * TILE: lane `tgx` (0..31) holds the row atom atom1 = x*TILE_SIZE+tgx for the whole tile, and the
 * 32 column atoms of block y are streamed past the warp lane by lane. `warp` is the global warp
 * index (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; `totalWarps` is the launch-wide warp count.
 * Each warp is statically assigned a contiguous half-open slice of every phase's tile range by
 * proportional partition (warp*N/totalWarps .. (warp+1)*N/totalWarps). Only the lower triangle
 * x >= y is enumerated; Newton's third law supplies the upper triangle by also accumulating the
 * equal-and-opposite reaction onto the column atom. This is a load-balanced grid-stride-like
 * partition, not a 1:1 thread:atom map -- an optimizer must preserve the warp-per-tile invariant
 * and the 32-wide lane assignment, not any particular tile-to-warp assignment.
 *
 * @par Three-phase structure (contract):
 *   Phase 1 -- exclusion tiles: linear indices [FIRST_EXCLUSION_TILE, LAST_EXCLUSION_TILE) into
 *     `exclusionTiles`. Every such tile computes the FULL TILE_SIZE^2 pair set and applies the
 *     per-pair exclusion bitmask. This phase runs for every configuration (with or without cutoff).
 *   Phase 2 -- non-exclusion tiles: with USE_NEIGHBOR_LIST, the interacting tiles produced by
 *     findInteractingBlocks (`tiles` + `interactingAtoms`, count `interactionCount[0]`); without a
 *     neighbor list, an enumeration of the entire lower triangle owned by this context
 *     ([startTileIndex, startTileIndex+numTileIndices)), skipping tiles already handled in phase 1
 *     via the `skipTiles` shared window. Each tile computes TILE_SIZE^2 ordered pairs once.
 *   Phase 3 -- single pairs (USE_NEIGHBOR_LIST only): leftover close pairs that findInteractingBlocks
 *     emitted individually into `singlePairs` (count `interactionCount[1]`) because their tile had
 *     too few interacting atoms to justify a full warp-tile. Here the mapping is 1 THREAD : 1 PAIR
 *     in a grid-stride loop; no shuffle, no exclusion mask (these pairs are never excluded).
 *
 * @par Diagonal vs off-diagonal tiles and the double-counting convention (contract):
 * On a DIAGONAL tile (x == y) atom1 and atom2 belong to the same block, so column data is obtained
 * by broadcasting lane j's registers across the warp (BROADCAST_WARP_DATA) rather than a second
 * global load, and all 32x32 ordered pairs (including i==j self-pairs, masked by exclusions) are
 * visited. Each unordered pair is therefore seen twice, so interactionScale = 0.5 and the energy
 * contribution is explicitly halved (energy += 0.5*tempEnergy); only atom1's half of the force is
 * kept (the reaction is the same warp's other visit). On an OFF-DIAGONAL tile (x > y) each lane
 * loads its own column atom once (LOAD_LOCAL_PARAMETERS_FROM_GLOBAL) and the warp rotates that
 * column data plus its accumulating reaction force by one lane each inner step (SHUFFLE_WARP_DATA);
 * all 32x32 ordered pairs are covered exactly once with interactionScale = 1.0. The rotation walks
 * the matrix diagonally so that at every step the 32 lanes touch 32 distinct column atoms -- no two
 * lanes accumulate onto the same `shflForce` simultaneously:
 *
 *        threads (lane tgx, holds atom1)
 *     0 1 2 3 4 5 6 7
 *         atom1
 * L    a b c d e f g h
 * o  i 1 2 3 4 5 6 7 8
 * c  j 8 1 2 3 4 5 6 7
 * a  k 7 8 1 2 3 4 5 6
 * l  l 6 7 8 1 2 3 4 5
 * D  m 5 6 7 8 1 2 3 4
 * a  n 4 5 6 7 8 1 2 3
 * t  o 3 4 5 6 7 8 1 2
 * a  p 2 3 4 5 6 7 8 1
 *
 * ([a-h] = block x row atoms, one per lane; [i-p] = block y column atoms, rotated; the digit is the
 * inner-loop step at which the owning lane computes that pair; 8-wide shown, real width 32.) atom1's
 * force accumulates in lane-local `force`; the column atom's reaction rides in `shflForce`, shuffled
 * along with the column data so it stays attached to its atom, and is written out once at tile end.
 *
 * @par Launch configuration:
 * Launched by computeInteractions as gridDim = min(numForceThreadBlocks, context numThreadBlocks),
 * blockDim = forceThreadBlockSize, 1-D grid, on the context's current stream. numForceThreadBlocks
 * = 4 * multiprocessor count; forceThreadBlockSize = THREAD_BLOCK_SIZE = 128 on compute capability
 * < 2.0 else 256 (set in the CudaNonbondedUtilities constructor). The kernel is oversubscription-
 * tolerant: it does NOT assume one warp per tile globally -- warps loop over their tile slice -- so
 * the grid is sized to fill the device, not to the problem. There is NO __launch_bounds__ on this
 * kernel: it makes no occupancy promise to the compiler, and register pressure (which grows with
 * the number of injected per-atom parameters and the complexity of COMPUTE_INTERACTION) is the
 * dominant occupancy lever an optimizer controls. THREAD_BLOCK_SIZE must be a multiple of TILE_SIZE.
 *
 * @par Shared memory:
 * Static only; the launch passes dynamic shared size 0 (executeKernel default). `atomIndices`
 * [THREAD_BLOCK_SIZE ints] is always allocated (holds each lane's phase-2 column atom index for the
 * force write-back). In the no-neighbor-list path an additional `skipTiles`
 * [THREAD_BLOCK_SIZE volatile ints] caches the sorted exclusion-tile-index window used to skip
 * phase-1 tiles during phase 2. Total static shared = THREAD_BLOCK_SIZE*4 bytes (neighbor-list
 * path) or *8 bytes (no-neighbor-list path). No dynamic shared memory formula scales with the
 * injected parameter set -- per-atom params live in registers (the shfl* copies), not shared memory.
 *
 * @par Memory spaces and SoA layout:
 * All pointer parameters are device global memory, borrowed for the launch (owned by
 * CudaNonbondedUtilities / CudaContext; none allocated or freed here). `posq` is AoS real4
 * (xyz+charge), length PADDED_NUM_ATOMS. `forceBuffers` is a fixed-point SoA accumulator with three
 * Cartesian planes at strides 0, PADDED_NUM_ATOMS, 2*PADDED_NUM_ATOMS -- PADDED_NUM_ATOMS is the row
 * stride and is a multiple of TILE_SIZE so partial tiles are always full-width. Each injected
 * global_<param> is its own SoA array indexed by atom (see PARAMETER_ARGUMENTS). Atoms in
 * [NUM_ATOMS, PADDED_NUM_ATOMS) are padding: masked out via isExcluded / the j < PADDED_NUM_ATOMS
 * guard, never written with meaningful force. Alignment is that of the element types (real4 => 16 B
 * for posq); no wider reinterpretation occurs.
 *
 * @par Warp shuffle, participation, and synchronization scope:
 * All intra-tile data movement is via warp shuffle (real_shfl / SHFL => __shfl_sync with a full
 * 0xffffffff mask on modern archs). Correctness REQUIRES FULL-WARP PARTICIPATION: every lane of the
 * warp must execute the inner tile loop in lockstep -- the includeTile branch and the padding
 * guards are warp-uniform so the warp does not diverge across a shuffle. There is NO __syncthreads
 * inside a tile and NO cross-block or cross-warp synchronization anywhere in this kernel; each warp
 * is independent and each block's shared arrays are private. Correctness of the shuffle rotation
 * rests on TILE_SIZE == warpSize (32); this is not portable to other warp widths without change.
 *
 * @par Cutoff and periodic boundary handling:
 * The distance cutoff test lives entirely inside the injected COMPUTE_INTERACTION snippet (it
 * compares r2 against CUTOFF_<g>_SQUARED); this kernel always computes r2, invR = RSQRT(r2), and
 * r = r2*invR for every visited pair and lets the snippet decide. Under USE_PERIODIC the per-pair
 * minimum image is applied to `delta` via APPLY_PERIODIC_TO_DELTA (rectangular per-axis rounding or
 * triclinic box-vector subtraction, selected by CudaContext boxIsTriclinic). The singlePeriodicCopy
 * fast path (phase 2, when the box half-width minus the block extent is >= MAX_CUTOFF on all axes)
 * pre-wraps atom1 and the column atoms once into the image nearest the block center
 * (APPLY_PERIODIC_TO_POS_WITH_CENTER) and then skips per-pair PBC entirely.
 *
 * @par Fixed-point accumulation and determinism (guarantee):
 * Forces are accumulated into `forceBuffers` as 64-bit signed fixed-point (scale 2^32, via
 * realToFixedPoint) with atomicAdd -- directly for tile forces, and via saveSingleForce for phase 3.
 * Integer atomicAdd is associative and commutative, so the force reduction is ORDER-INDEPENDENT and
 * BIT-EXACT REPRODUCIBLE across runs regardless of warp scheduling or neighbor-list ordering. This
 * is the determinism guarantee an optimizer must not break: any reordering of atomic contributions
 * is safe, but switching to floating-point accumulation would forfeit it. Energy is accumulated per
 * thread in `mixed` (energyBuffer), reduced elsewhere; that reduction is not fixed-point.
 *
 * @par Performance flags (compile-time, from createInteractionKernel):
 *   USE_CUTOFF         adds the neighbor-list argument block; without it, phase 2 enumerates the
 *                      full lower triangle and phase 3 is compiled out.
 *   USE_PERIODIC       enables minimum-image / singlePeriodicCopy handling.
 *   USE_EXCLUSIONS     enables the per-tile TILE_SIZE-bit exclusion mask path (always on here).
 *   USE_NEIGHBOR_LIST  phase 2 reads the interacting-tile list and enables phase 3.
 *   USE_SYMMETRIC      central-force fast path: COMPUTE_INTERACTION yields a scalar dEdR along
 *                      `delta`; the reaction is -dEdR*delta. When undefined, the snippet fills
 *                      independent dEdR1/dEdR2 vectors (non-central forces).
 *   INCLUDE_FORCES / INCLUDE_ENERGY  compiled separately for the force / energy / force+energy
 *                      variants (computeInteractions selects the CUfunction); each gates its writes.
 *   ENABLE_SHUFFLE     always 1. TILE_SIZE == 32, THREAD_BLOCK_SIZE, PADDED_NUM_ATOMS, NUM_ATOMS,
 *                      NUM_BLOCKS, NUM_TILES_WITH_EXCLUSIONS, FIRST/LAST_EXCLUSION_TILE, MAX_CUTOFF,
 *                      CUTOFF_<g>[_SQUARED]  are scalar `#define`s (see file-header section 1).
 *
 * @param forceBuffers    [out] device SoA fixed-point force accumulator, 3*PADDED_NUM_ATOMS entries
 *                        (planes at offsets 0, PADDED_NUM_ATOMS, 2*PADDED_NUM_ATOMS). Read-modify-
 *                        write via atomicAdd; written only if INCLUDE_FORCES. Borrowed.
 * @param energyBuffer    [out] device per-thread `mixed` energy accumulator (one slot per global
 *                        thread; reduced elsewhere). Read-modify-write; written only if
 *                        INCLUDE_ENERGY. Borrowed.
 * @param posq            [in]  device AoS real4 per atom (xyz position + w charge), length
 *                        PADDED_NUM_ATOMS. Read-only, borrowed; must outlive the launch.
 * @param exclusions      [in]  device tileflags (32-bit) column masks, TILE_SIZE per exclusion tile:
 *                        exclusions[pos*TILE_SIZE + tgx] is row atom tgx's mask (bit k set => the
 *                        pair with column atom k is INCLUDED). Read-only, borrowed.
 * @param exclusionTiles  [in]  device int2 {x,y} block-index pairs (x >= y) each containing at least
 *                        one exclusion. Length NUM_TILES_WITH_EXCLUSIONS. Read-only, borrowed.
 * @param startTileIndex  [in]  first lower-triangle linear tile index this context owns (no-neighbor-
 *                        list path); de-facto constant per launch.
 * @param numTileIndices  [in]  count of tiles this context owns (no-neighbor-list path).
 * @param tiles           [in]  (USE_CUTOFF) device block index x per interacting tile; column atoms
 *                        in `interactingAtoms`. Read-only, borrowed.
 * @param interactionCount[in]  (USE_CUTOFF) device [numInteractingTiles, numSinglePairs] filled by
 *                        findInteractingBlocks. If numInteractingTiles > maxTiles or
 *                        numSinglePairs > maxSinglePairs the affected phase returns early (partial or
 *                        no work) and the caller (updateNeighborListSize) grows the arrays and marks
 *                        forces invalid to retry -- rare after the first step. Read-only, borrowed.
 * @param periodicBoxSize / invPeriodicBoxSize / periodicBoxVecX,Y,Z  [in] (USE_CUTOFF) box metrics
 *                        (real4; w unused for the size vectors). By-value; consumed by the periodic
 *                        macros.
 * @param maxTiles        [in]  (USE_CUTOFF) capacity of tiles/interactingAtoms; overflow guard.
 * @param blockCenter     [in]  (USE_CUTOFF) device real4 per block: geometric center (xyz). Read-only.
 * @param blockSize       [in]  (USE_CUTOFF) device real4 per block: half-extents (xyz; w unused),
 *                        used for the singlePeriodicCopy test. Read-only, borrowed.
 * @param interactingAtoms[in]  (USE_CUTOFF) device column atom indices, TILE_SIZE per tile:
 *                        interactingAtoms[pos*TILE_SIZE + tgx] is this lane's column atom for tile
 *                        `pos`. Read-only, borrowed.
 * @param maxSinglePairs  [in]  (USE_CUTOFF) capacity of singlePairs; overflow guard.
 * @param singlePairs     [in]  (USE_CUTOFF) device int2 {atom1,atom2} pairs handled in phase 3.
 *                        Read-only, borrowed.
 * @param PARAMETER_ARGUMENTS  [in] injected trailing formals (see file-header section 2): one
 *                        device SoA global_<param> array per per-atom force parameter, then extra
 *                        `arguments`, then optionally `mixed* energyParamDerivs` [out]. All borrowed.
 *
 * @pre The injected snippets (COMPUTE_INTERACTION and the LOAD_/DECLARE_/SHUFFLE_/BROADCAST_/
 *      CLEAR_/SAVE_ macros) MUST have been generated by createInteractionKernel from EXACTLY this
 *      force group's parameter set and Lepton expression: the register names P1/P2/shflP and the
 *      global_<param> formals produced by the ComputeParameterInfo loop must agree in name, type,
 *      and component count. This tag/type consistency is established at the createInteractionKernel
 *      dispatch boundary and is not verifiable by the compiler after substitution; a mismatch is
 *      undefined behavior.
 * @pre THREAD_BLOCK_SIZE is a multiple of TILE_SIZE and TILE_SIZE == warp size (32).
 * @pre With USE_NEIGHBOR_LIST, prepareInteractions has already run findInteractingBlocks for this
 *      step (interactionCount, tiles, interactingAtoms, singlePairs are current).
 * @pre `posq` holds current atom positions; forceBuffers/energyBuffer were zeroed before launch
 *      (the accumulation is additive).
 *
 * @post forceBuffers is incremented by this group's fixed-point forces (if INCLUDE_FORCES),
 *       energyBuffer by its energy (if INCLUDE_ENERGY), and energyParamDerivs by any declared
 *       energy-parameter derivatives (SAVE_DERIVATIVES). No global synchronization is performed; the
 *       results are complete only after the caller synchronizes the launch stream.
 *
 * @note Determinism: fixed-point 2^32 force atomics are order-independent and bitwise reproducible.
 * @note This kernel never allocates or frees device memory and performs no host<->device transfer.
 * @see saveSingleForce, real_shfl, CudaNonbondedUtilities::createInteractionKernel,
 *      CudaNonbondedUtilities::computeInteractions.
 */
extern "C" __global__ void computeNonbonded(
        unsigned long long* __restrict__ forceBuffers, mixed* __restrict__ energyBuffer, const real4* __restrict__ posq, const tileflags* __restrict__ exclusions,
        const int2* __restrict__ exclusionTiles, unsigned int startTileIndex, unsigned long long numTileIndices
#ifdef USE_CUTOFF
        , const int* __restrict__ tiles, const unsigned int* __restrict__ interactionCount, real4 periodicBoxSize, real4 invPeriodicBoxSize, 
        real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, unsigned int maxTiles, const real4* __restrict__ blockCenter,
        const real4* __restrict__ blockSize, const unsigned int* __restrict__ interactingAtoms, unsigned int maxSinglePairs,
        const int2* __restrict__ singlePairs
#endif
        PARAMETER_ARGUMENTS) {
    const unsigned int totalWarps = (blockDim.x*gridDim.x)/TILE_SIZE;
    const unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/TILE_SIZE; // global warpIndex
    const unsigned int tgx = threadIdx.x & (TILE_SIZE-1); // index within the warp
    const unsigned int tbx = threadIdx.x - tgx;           // block warpIndex
    mixed energy = 0;
    INIT_DERIVATIVES

    // First loop: process tiles that contain exclusions.
    // [FIRST_EXCLUSION_TILE, LAST_EXCLUSION_TILE) is THIS context's slice of the global
    // exclusionTiles list; createInteractionKernel derives it by contextIndex/numContexts so
    // multi-GPU runs partition the exclusion tiles disjointly (single-GPU => the whole list).

    const unsigned int firstExclusionTile = FIRST_EXCLUSION_TILE+warp*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    const unsigned int lastExclusionTile = FIRST_EXCLUSION_TILE+(warp+1)*(LAST_EXCLUSION_TILE-FIRST_EXCLUSION_TILE)/totalWarps;
    for (int pos = firstExclusionTile; pos < lastExclusionTile; pos++) {
        const int2 tileIndices = exclusionTiles[pos];
        const unsigned int x = tileIndices.x;
        const unsigned int y = tileIndices.y;
        real3 force = make_real3(0);
        unsigned int atom1 = x*TILE_SIZE + tgx;
        real4 posq1 = posq[atom1];
        // First injection point of the per-force snippet family. LOAD_ATOM1_PARAMETERS (and its
        // siblings COMPUTE_INTERACTION / BROADCAST_/SHUFFLE_/LOAD_/DECLARE_/CLEAR_LOCAL_PARAMETERS)
        // are force-dependent source built by createInteractionKernel from this group's parameter
        // list + Lepton expression; here it declares the P1 registers COMPUTE_INTERACTION reads.
        LOAD_ATOM1_PARAMETERS
#ifdef USE_EXCLUSIONS
        tileflags excl = exclusions[pos*TILE_SIZE+tgx];
#endif
        // Read by the injected COMPUTE_INTERACTION: signals it that this tile may carry excluded
        // pairs, so it must honor `isExcluded`. Phase 2/3 set it false (no exclusions possible).
        const bool hasExclusions = true;
        if (x == y) {
            // This tile is on the diagonal.
            real4 shflPosq = posq1;

            // we do not need to fetch parameters from global since this is a symmetric tile
            // instead we can broadcast the values using shuffle
            for (unsigned int j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+j;
                real4 posq2;
                BROADCAST_WARP_DATA
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = r2*invR;
                LOAD_ATOM2_PARAMETERS
                // atom2 was the warp-lane address (tbx+j) used to fetch posq2/params via shuffle;
                // now reassign it to the true global atom index for exclusion/padding tests.
                atom2 = y*TILE_SIZE+j;
#ifdef USE_SYMMETRIC
                real dEdR = 0.0f;
#else
                real3 dEdR1 = make_real3(0);
                real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
#endif
                real tempEnergy = 0.0f;
                // Diagonal tile visits every unordered pair twice (both i,j and j,i); interactionScale
                // = 0.5 tells the snippet to halve any per-pair scalar it forms, and energy is halved
                // here for the same reason. Off-diagonal tiles use 1.0 (each ordered pair seen once).
                const real interactionScale = 0.5f;
                COMPUTE_INTERACTION
                energy += 0.5f*tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                force.x -= delta.x*dEdR;
                force.y -= delta.y*dEdR;
                force.z -= delta.z*dEdR;
#else
                force.x -= dEdR1.x;
                force.y -= dEdR1.y;
                force.z -= dEdR1.z;
#endif
#endif
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
            }
        }
        else {
            // This is an off-diagonal tile.
            unsigned int j = y*TILE_SIZE + tgx;
            real4 shflPosq = posq[j];
            real3 shflForce;
            shflForce.x = 0.0f;
            shflForce.y = 0.0f;
            shflForce.z = 0.0f;
            DECLARE_LOCAL_PARAMETERS
            LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
#ifdef USE_EXCLUSIONS
            // Cyclically pre-rotate this lane's column mask by tgx so that bit 0 lines up with the
            // first atom2 this lane sees (tj starts at tgx). Then `excl >>= 1` each inner step keeps
            // the low bit aligned with the current atom2 as the warp rotates.
            excl = (excl >> tgx) | (excl << (TILE_SIZE - tgx));
#endif
            unsigned int tj = tgx;
            for (j = 0; j < TILE_SIZE; j++) {
                int atom2 = tbx+tj;
                real4 posq2 = shflPosq;
                real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                APPLY_PERIODIC_TO_DELTA(delta)
#endif
                real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                real invR = RSQRT(r2);
                real r = r2*invR;
                LOAD_ATOM2_PARAMETERS
                atom2 = y*TILE_SIZE+tj;
#ifdef USE_SYMMETRIC
                real dEdR = 0.0f;
#else
                real3 dEdR1 = make_real3(0);
                real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS || !(excl & 0x1));
#endif
                real tempEnergy = 0.0f;
                const real interactionScale = 1.0f;
                COMPUTE_INTERACTION
                energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                delta *= dEdR;
                force.x -= delta.x;
                force.y -= delta.y;
                force.z -= delta.z;
                shflForce.x += delta.x;
                shflForce.y += delta.y;
                shflForce.z += delta.z;
#else // !USE_SYMMETRIC
                force.x -= dEdR1.x;
                force.y -= dEdR1.y;
                force.z -= dEdR1.z;
                shflForce.x += dEdR2.x;
                shflForce.y += dEdR2.y;
                shflForce.z += dEdR2.z;
#endif // end USE_SYMMETRIC
#endif
                SHUFFLE_WARP_DATA
#ifdef USE_EXCLUSIONS
                excl >>= 1;
#endif
                // cycles the indices
                // 0 1 2 3 4 5 6 7 -> 1 2 3 4 5 6 7 0
                tj = (tj + 1) & (TILE_SIZE - 1);
            }
            const unsigned int offset = y*TILE_SIZE + tgx;
            // write results for off diagonal tiles
            // realToFixedPoint scales by 2^32 into signed 64-bit; forceBuffers is context.getForce(),
            // a fixed-point SoA buffer the host reduces/interprets as fixed-point. Keeping this an
            // INTEGER atomicAdd is what makes the force sum order-independent and bit-reproducible --
            // an optimizer must not swap it for float accumulation.
#ifdef INCLUDE_FORCES
            atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>(realToFixedPoint(shflForce.x)));
            atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.y)));
            atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.z)));
#endif
        }
        // Write results for on and off diagonal tiles
#ifdef INCLUDE_FORCES
        const unsigned int offset = x*TILE_SIZE + tgx;
        atomicAdd(&forceBuffers[offset], static_cast<unsigned long long>(realToFixedPoint(force.x)));
        atomicAdd(&forceBuffers[offset+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
        atomicAdd(&forceBuffers[offset+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
#endif
    }

    // Second loop: tiles without exclusions, either from the neighbor list (with cutoff) or just enumerating all
    // of them (no cutoff).

#ifdef USE_NEIGHBOR_LIST
    const unsigned int numTiles = interactionCount[0];
    // Overflow: findInteractingBlocks found more tiles than the arrays hold. Bail out (forces stay
    // partial); the host updateNeighborListSize sees interactionCount > maxTiles, grows the arrays
    // ~1.2x, calls setForcesValid(false), and reruns the step -- so this early return is recoverable.
    if (numTiles > maxTiles)
        return; // There wasn't enough memory for the neighbor list.
    int pos = (int) (warp*(long long)numTiles/totalWarps);
    int end = (int) ((warp+1)*(long long)numTiles/totalWarps);
#else
    int pos = (int) (startTileIndex+warp*numTileIndices/totalWarps);
    int end = (int) (startTileIndex+(warp+1)*numTileIndices/totalWarps);
    int skipBase = 0;
    int currentSkipIndex = tbx;
    __shared__ volatile int skipTiles[THREAD_BLOCK_SIZE];
    skipTiles[threadIdx.x] = -1;
#endif
    // atomIndices can probably be shuffled as well
    // but it probably wouldn't make things any faster
    __shared__ int atomIndices[THREAD_BLOCK_SIZE];
    
    while (pos < end) {
        const bool hasExclusions = false;
        real3 force = make_real3(0);
        bool includeTile = true;

        // Extract the coordinates of this tile.
        int x, y;
        bool singlePeriodicCopy = false;
#ifdef USE_NEIGHBOR_LIST
        x = tiles[pos];
        // The `blockSize` formal is bound (in initialize's forceArgs) to the blockBoundingBox array,
        // so blockSizeX holds block x's half-extents (xyz), not a scalar size.
        real4 blockSizeX = blockSize[x];
        singlePeriodicCopy = (0.5f*periodicBoxSize.x-blockSizeX.x >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.y-blockSizeX.y >= MAX_CUTOFF &&
                              0.5f*periodicBoxSize.z-blockSizeX.z >= MAX_CUTOFF);
#else
        // No neighbor list: invert the linear tile index `pos` into (x,y) lower-triangle block
        // coordinates. The float sqrt can be off by one at boundaries, hence the correction below.
        y = (int) floor(NUM_BLOCKS+0.5f-SQRT((NUM_BLOCKS+0.5f)*(NUM_BLOCKS+0.5f)-2*pos));
        x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        if (x < y || x >= NUM_BLOCKS) { // Occasionally happens due to roundoff error.
            y += (x < y ? -1 : 1);
            x = (pos-y*NUM_BLOCKS+y*(y+1)/2);
        }

        // Skip over tiles that have exclusions, since they were already processed in phase 1.
        // skipTiles caches, per warp, a sorted window of linear indices of exclusion tiles; the
        // warp advances `currentSkipIndex` until it reaches/passes `pos` and drops the tile if it
        // matches. Each lane loads one exclusion tile's linear index (x + y*NUM_BLOCKS - y(y+1)/2).

        while (skipTiles[tbx+TILE_SIZE-1] < pos) {
            if (skipBase+tgx < NUM_TILES_WITH_EXCLUSIONS) {
                int2 tile = exclusionTiles[skipBase+tgx];
                skipTiles[threadIdx.x] = tile.x + tile.y*NUM_BLOCKS - tile.y*(tile.y+1)/2;
            }
            else
                skipTiles[threadIdx.x] = end;
            skipBase += TILE_SIZE;            
            currentSkipIndex = tbx;
        }
        while (skipTiles[currentSkipIndex] < pos)
            currentSkipIndex++;
        includeTile = (skipTiles[currentSkipIndex] != pos);
#endif
        if (includeTile) {
            unsigned int atom1 = x*TILE_SIZE + tgx;
            // Load atom data for this tile.
            real4 posq1 = posq[atom1];
            LOAD_ATOM1_PARAMETERS
#ifdef USE_NEIGHBOR_LIST
            unsigned int j = interactingAtoms[pos*TILE_SIZE+tgx];
#else
            unsigned int j = y*TILE_SIZE + tgx;
#endif
            atomIndices[threadIdx.x] = j;
            DECLARE_LOCAL_PARAMETERS
            real4 shflPosq;
            real3 shflForce;
            shflForce.x = 0.0f;
            shflForce.y = 0.0f;
            shflForce.z = 0.0f;
            if (j < PADDED_NUM_ATOMS) {
                // Load position of atom j from from global memory
                shflPosq = posq[j];
                LOAD_LOCAL_PARAMETERS_FROM_GLOBAL
            }
            else {
                shflPosq = make_real4(0, 0, 0, 0);
                CLEAR_LOCAL_PARAMETERS
            }
#ifdef USE_PERIODIC
            if (singlePeriodicCopy) {
                // The box is small enough that we can just translate all the atoms into a single periodic
                // box, then skip having to apply periodic boundary conditions later.
                real4 blockCenterX = blockCenter[x];
                APPLY_PERIODIC_TO_POS_WITH_CENTER(posq1, blockCenterX)
                APPLY_PERIODIC_TO_POS_WITH_CENTER(shflPosq, blockCenterX)
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real4 posq2 = shflPosq; 
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
#ifdef USE_SYMMETRIC
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
#endif
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                    delta *= dEdR;
                    force.x -= delta.x;
                    force.y -= delta.y;
                    force.z -= delta.z;
                    shflForce.x += delta.x;
                    shflForce.y += delta.y;
                    shflForce.z += delta.z;
#else // !USE_SYMMETRIC
                    force.x -= dEdR1.x;
                    force.y -= dEdR1.y;
                    force.z -= dEdR1.z;
                    shflForce.x += dEdR2.x;
                    shflForce.y += dEdR2.y;
                    shflForce.z += dEdR2.z;
#endif // end USE_SYMMETRIC
#endif
                    SHUFFLE_WARP_DATA
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }
            else
#endif
            {
                // We need to apply periodic boundary conditions separately for each interaction.
                unsigned int tj = tgx;
                for (j = 0; j < TILE_SIZE; j++) {
                    int atom2 = tbx+tj;
                    real4 posq2 = shflPosq;
                    real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
                    APPLY_PERIODIC_TO_DELTA(delta)
#endif
                    real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
                    real invR = RSQRT(r2);
                    real r = r2*invR;
                    LOAD_ATOM2_PARAMETERS
                    atom2 = atomIndices[tbx+tj];
#ifdef USE_SYMMETRIC
                    real dEdR = 0.0f;
#else
                    real3 dEdR1 = make_real3(0);
                    real3 dEdR2 = make_real3(0);
#endif
#ifdef USE_EXCLUSIONS
                    bool isExcluded = (atom1 >= NUM_ATOMS || atom2 >= NUM_ATOMS);
#endif
                    real tempEnergy = 0.0f;
                    const real interactionScale = 1.0f;
                    COMPUTE_INTERACTION
                    energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
                    delta *= dEdR;
                    force.x -= delta.x;
                    force.y -= delta.y;
                    force.z -= delta.z;
                    shflForce.x += delta.x;
                    shflForce.y += delta.y;
                    shflForce.z += delta.z;
#else // !USE_SYMMETRIC
                    force.x -= dEdR1.x;
                    force.y -= dEdR1.y;
                    force.z -= dEdR1.z;
                    shflForce.x += dEdR2.x;
                    shflForce.y += dEdR2.y;
                    shflForce.z += dEdR2.z;
#endif // end USE_SYMMETRIC
#endif
                    SHUFFLE_WARP_DATA
                    tj = (tj + 1) & (TILE_SIZE - 1);
                }
            }

            // Write results.
#ifdef INCLUDE_FORCES
            atomicAdd(&forceBuffers[atom1], static_cast<unsigned long long>(realToFixedPoint(force.x)));
            atomicAdd(&forceBuffers[atom1+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.y)));
            atomicAdd(&forceBuffers[atom1+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(force.z)));
#ifdef USE_NEIGHBOR_LIST
            unsigned int atom2 = atomIndices[threadIdx.x];
#else
            unsigned int atom2 = y*TILE_SIZE + tgx;
#endif
            if (atom2 < PADDED_NUM_ATOMS) {
                atomicAdd(&forceBuffers[atom2], static_cast<unsigned long long>(realToFixedPoint(shflForce.x)));
                atomicAdd(&forceBuffers[atom2+PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.y)));
                atomicAdd(&forceBuffers[atom2+2*PADDED_NUM_ATOMS], static_cast<unsigned long long>(realToFixedPoint(shflForce.z)));
            }
#endif
        }
        pos++;
    }
    
    // Third loop: single pairs that aren't part of a tile. findInteractingBlocks emits these when a
    // tile has too few interacting atoms to be worth a full warp-tile; here each thread handles one
    // (atom1,atom2) pair directly (no shuffle, no exclusion mask -- these pairs are never excluded).
    // Forces go through saveSingleForce (same 2^32 fixed-point atomic accumulation).

#if USE_NEIGHBOR_LIST
    const unsigned int numPairs = interactionCount[1];
    if (numPairs > maxSinglePairs)
        return; // There wasn't enough memory for the neighbor list.
    for (int i = blockIdx.x*blockDim.x+threadIdx.x; i < numPairs; i += blockDim.x*gridDim.x) {
        int2 pair = singlePairs[i];
        int atom1 = pair.x;
        int atom2 = pair.y;
        real4 posq1 = posq[atom1];
        real4 posq2 = posq[atom2];
        LOAD_ATOM1_PARAMETERS
        LOAD_ATOM2_PARAMETERS_FROM_GLOBAL
        real3 delta = make_real3(posq2.x-posq1.x, posq2.y-posq1.y, posq2.z-posq1.z);
#ifdef USE_PERIODIC
        APPLY_PERIODIC_TO_DELTA(delta)
#endif
        real r2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
        real invR = RSQRT(r2);
        real r = r2*invR;
#ifdef USE_SYMMETRIC
        real dEdR = 0.0f;
#else
        real3 dEdR1 = make_real3(0);
        real3 dEdR2 = make_real3(0);
#endif
        bool hasExclusions = false;
        bool isExcluded = false;
        real tempEnergy = 0.0f;
        const real interactionScale = 1.0f;
        COMPUTE_INTERACTION
        energy += tempEnergy;
#ifdef INCLUDE_FORCES
#ifdef USE_SYMMETRIC
        real3 dEdR1 = delta*dEdR;
        real3 dEdR2 = -dEdR1;
#endif
        saveSingleForce(atom1, -dEdR1, forceBuffers);
        saveSingleForce(atom2, -dEdR2, forceBuffers);
#endif
    }
#endif
#ifdef INCLUDE_ENERGY
    energyBuffer[blockIdx.x*blockDim.x+threadIdx.x] += energy;
#endif
    SAVE_DERIVATIVES
}
