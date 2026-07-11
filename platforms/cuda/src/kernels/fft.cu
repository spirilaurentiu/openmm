/**
 * @file fft.cu
 * @brief CUDA-syntax template for one pass of OpenMM's hand-written batched-1D
 *        mixed-radix FFT, the building block of the 3D FFT used for PME
 *        reciprocal-space electrostatics.
 *
 * @par Provenance and build status (READ FIRST)
 * In this vendored OpenMM the CUDA platform's 3D FFT is class CudaFFT3D
 * (CudaFFT3D.cpp), implemented entirely with cuFFT: its constructor calls
 * cufftPlan3d and its execFFT dispatches cufftExec{R2C,C2R,Z2Z,C2C,D2Z,Z2D}. No
 * CUDA C++ driver reads CudaKernelSources::fft or substitutes any of the macros
 * this file depends on. This kernel is therefore ORPHANED in the CUDA build: it is
 * stringified into CudaKernelSources::fft by the source-embedding step but is never
 * compiled by NVRTC and never launched. It is the CUDA-syntax twin of the live
 * OpenCL kernel fft.cl. The authoritative injector for every macro below is the
 * OpenCL driver method OpenCLFFT3D::createKernel (OpenCLFFT3D.cpp, the non-USE_VKFFT
 * branch), which JIT-compiles fft.cl after populating a replacements map. Read each
 * macro's semantics as "what createKernel would inject"; the CUDA-only macro
 * THREADS_PER_BLOCK has no in-tree injector and is documented as an assumption.
 * Kernel-source comments are stripped before JIT (strip_comments.py), so this
 * documentation carries zero runtime cost.
 *
 * @par Algorithm (contract, not mechanism)
 * A 3D FFT is realized as three passes of batched 1D FFTs, one per axis. Each pass
 * transforms every "row" (a line along the transform axis) with a Stockham
 * self-sorting mixed-radix Cooley-Tukey FFT, and writes the result transposed so
 * the following pass reads its own axis contiguously. The driver permutes the
 * logical grid so the axis under transform is always local Z (innermost); it
 * factors the row length ZSIZE into primes drawn from {2,3,4,5,7} and emits one
 * straight-line radix butterfly block per factor as the COMPUTE_FFT string. This
 * file supplies only the twiddle-table build, the shared-memory staging, and the
 * launch scaffold around that generated body; the butterflies and the final store
 * live entirely inside COMPUTE_FFT. Per 3D transform, createKernel builds six
 * specializations (z/x/y axis x forward/inverse), each a distinct macro set.
 *
 * @par Compile-time transform parameters (each: meaning + how createKernel chooses it)
 * String-substituted before JIT; clangd cannot resolve them. Types are the C type
 * of the substituted literal.
 *
 * - real, real2      : NOT injected here; come from the context's global precision
 *                      defines. `real` is float or double; `real2` the matching
 *                      2-vector (float2/double2). A complex value is real2 with
 *                      re in .x, im in .y. sizeof(real2) is 8 (single) or 16 (double)
 *                      bytes, which sets the shared-memory footprint below.
 *
 * - XSIZE, YSIZE, ZSIZE : int. Dimensions of the axis-permuted grid THIS pass sees.
 *                      The transform runs along ZSIZE; XSIZE*YSIZE independent rows
 *                      of length ZSIZE are batched. createKernel receives them as the
 *                      permuted (xsize,ysize,zsize) arguments and injects them via
 *                      replacements["XSIZE"/"YSIZE"/"ZSIZE"]. ZSIZE governs both the
 *                      twiddle table and the shared-array sizing, so it is the single
 *                      most performance-relevant parameter for occupancy tuning.
 *
 * - SIGN             : int literal +1 (forward) or -1 (inverse); replacements["SIGN"]
 *                      = forward ? "1" : "-1". Sets the twiddle exponent sign in
 *                      w[k] = exp(-i*SIGN*2*M_PI*k/ZSIZE) and the imaginary-rotation
 *                      sign of every butterfly. The transform is UNNORMALIZED: a
 *                      forward followed by an inverse scales data by the point count.
 *
 * - INPUT_TYPE       : `real` or `real2`; element type of @p in. `real` only for the
 *                      forward pass over the innermost real axis of an unpacked
 *                      real-to-complex FFT (i.e. inputIsReal && axis==0 && forward),
 *                      else `real2`.
 *
 * - OUTPUT_TYPE      : `real` or `real2`; element type of @p out. `real` only when
 *                      writing the final real result of an inverse real transform
 *                      (inputIsReal && axis==2 && !forward), else `real2`.
 *
 * - INPUT_IS_REAL    : 0/1 preprocessor guard. 1 iff (inputIsReal && axis==0 &&
 *                      forward): @p in is a plain real array, each sample widened to
 *                      real2{v,0} on load. Mutually exclusive with the next two.
 *
 * - INPUT_IS_PACKED  : 0/1 guard. 1 iff (inputIsReal && axis==0 && !forward): @p in is
 *                      the Hermitian half-complex grid from a forward real transform,
 *                      loaded through loadComplexValue which rebuilds the full spectrum
 *                      by conjugate symmetry.
 *
 * - OUTPUT_IS_PACKED : 0/1 guard. 1 iff (inputIsReal && axis==2 && forward): only the
 *                      non-redundant Hermitian half along local X (rows x < XSIZE/2+1)
 *                      is loaded/transformed; redundant conjugate rows are skipped.
 *                      The narrowed store itself is emitted inside COMPUTE_FFT.
 *
 * - BLOCKS_PER_GROUP : int >= 1; replacements["BLOCKS_PER_GROUP"]. Number of
 *                      independent length-ZSIZE rows packed into one thread block, so
 *                      short rows share a block for occupancy. createKernel sets it to
 *                      1 when loopRequired (row longer than the device max block size,
 *                      or CPU device), else max(1, maxThreads/zsize). It sizes the
 *                      data0/data1 shared arrays (BLOCKS_PER_GROUP*ZSIZE each).
 *
 * - COMPUTE_FFT      : multi-line code string; replacements["COMPUTE_FFT"]. The
 *                      per-radix butterfly passes followed by the transposed (and
 *                      possibly Hermitian-packed / real-narrowed) store to @p out.
 *                      Ping-pongs between shared arrays data0 and data1, reading
 *                      data(stage%2) and writing data(1-stage%2) each pass, and
 *                      multiplies by w[] via multiplyComplex.
 *
 * - M_PI             : double literal for pi; replacements["M_PI"] =
 *                      doubleToString(M_PI). Lets device code avoid math.h.
 *
 * @note Assumed: THREADS_PER_BLOCK. This CUDA-only macro (the number of threads
 *       cooperating on ONE length-ZSIZE row; total block width is
 *       BLOCKS_PER_GROUP*THREADS_PER_BLOCK, and thread i belongs to sub-block
 *       threadIdx.x/THREADS_PER_BLOCK) has NO injector in this tree: the OpenCL
 *       kernel uses get_local_size(0) and its ZSIZE-derived slicing instead, and the
 *       CUDA driver that once supplied THREADS_PER_BLOCK was removed with the cuFFT
 *       switch. By analogy to OpenCL's threads = blocksPerGroup*zsize, its intended
 *       value is ZSIZE (one thread per row element). Treated as an assumption; logged
 *       in findings as untraceable.
 *
 * @note createKernel also injects LOOP_REQUIRED (0/1). The staging scaffold in this
 *       file corresponds to the non-looped path (one thread per row element); the
 *       looped variant, for rows exceeding the block size, is realized through the
 *       generated COMPUTE_FFT and the LOOP_REQUIRED path of the .cl template.
 *
 * @note PACKED_AXIS / PACKED_XSIZE / PACKED_YSIZE / PACKED_ZSIZE are NOT referenced
 *       by this kernel; they belong to the sibling kernel fftR2C.cu and are not part
 *       of this file's contract.
 *
 * @see OpenCLFFT3D::createKernel (macro injector), fft.cl (the live twin),
 *      CudaFFT3D (the cuFFT path that actually runs on CUDA).
 */

/**
 * @brief Complex product of two interleaved real2 operands (.x real, .y imag).
 *
 * Returns (a+bi)(c+di) = (ac-bd) + (ad+bc)i. Applied inside COMPUTE_FFT to twiddle
 * each butterfly output against the shared twiddle table w[]. Pure and register-only;
 * no memory, sync, or divergence effects.
 *
 * @param[in] c1 first complex operand.
 * @param[in] c2 second complex operand.
 * @return the complex product c1*c2.
 */
static __inline__ __device__ real2 multiplyComplex(real2 c1, real2 c2) {
    return make_real2(c1.x*c2.x-c1.y*c2.y, c1.x*c2.y+c1.y*c2.x);
}

/**
 * @brief Reconstruct full-spectrum element F(x,y,z) from a Hermitian half-complex grid.
 *
 * A forward real-to-complex transform stores only the non-redundant half of the
 * spectrum along local Z (z in [0, ZSIZE/2]), i.e. inputZSize = ZSIZE/2+1 planes laid
 * out contiguously as in[x*YSIZE*inputZSize + y*inputZSize + z]. The inverse pass
 * needs the full length-ZSIZE spectrum, so for z > ZSIZE/2 this returns the reflected,
 * conjugated value F(x,y,z) = conj(F(-x mod X, -y mod Y, -z mod Z)); negated indices
 * fold into [0,size) as (size - idx) with idx==0 mapping to itself. Reached only when
 * INPUT_IS_PACKED (inverse pass over the real axis).
 *
 * @par Memory access
 * One global read of @p in per call. For z <= ZSIZE/2 the read is the natural
 * contiguous slot; for z > ZSIZE/2 it is a reflected, generally non-coalesced slot
 * (relevant to profilers: the upper half of every packed inverse row scatters loads).
 *
 * @param[in] in device pointer to the packed half-complex grid; borrowed, read-only,
 *               z-minor with z-stride inputZSize = ZSIZE/2+1. Must be sized for at
 *               least XSIZE*YSIZE*inputZSize real2 elements.
 * @param[in] x  full-grid first index in [0, XSIZE).
 * @param[in] y  full-grid second index in [0, YSIZE).
 * @param[in] z  full-grid third index in [0, ZSIZE); z > ZSIZE/2 triggers the
 *               conjugate-symmetry reflection.
 * @return the complex spectrum value F(x,y,z).
 */
static __inline__ __device__ real2 loadComplexValue(const real2* __restrict__ in, int x, int y, int z) {
    const int inputZSize = ZSIZE/2+1;
    if (z < inputZSize)
        return in[x*YSIZE*inputZSize+y*inputZSize+z];
    int xp = (x == 0 ? 0 : XSIZE-x);
    int yp = (y == 0 ? 0 : YSIZE-y);
    real2 value = in[xp*YSIZE*inputZSize+yp*inputZSize+(ZSIZE-z)];
    return make_real2(value.x, -value.y);
}

/**
 * @brief One 3D-FFT pass: batched 1D mixed-radix FFT of every length-ZSIZE row, stored
 *        transposed for the next axis.
 *
 * Treats the permuted grid as XSIZE*YSIZE independent rows of length ZSIZE and
 * transforms each along local Z. It builds the twiddle table w[] in shared memory,
 * stages one row per sub-block into shared data0, then runs the injected COMPUTE_FFT
 * (radix-2/3/4/5/7 butterflies ping-ponging data0<->data1 and twiddling via
 * multiplyComplex against w[]) whose final statement writes the transposed result to
 * @p out. Evaluates the N-point DFT X[k] = sum_n x[n]*exp(-i*SIGN*2*M_PI*n*k/N) by
 * Cooley-Tukey factorization of N=ZSIZE into radix stages.
 *
 * @par Launch configuration
 * 1D launch. blockDim.x = BLOCKS_PER_GROUP*THREADS_PER_BLOCK (assumed
 * BLOCKS_PER_GROUP*ZSIZE; see the file-level THREADS_PER_BLOCK note); one sub-block of
 * THREADS_PER_BLOCK threads cooperates on one length-ZSIZE row, BLOCKS_PER_GROUP rows
 * per block. gridDim.x is unconstrained: the grid-stride loop over baseIndex covers all
 * XSIZE*YSIZE rows regardless of grid size. In the live OpenCL twin, the driver's
 * executeKernel launches a total work size of XSIZE*YSIZE*ZSIZE (halved when packing)
 * with local size blocksPerGroup*zsize.
 *
 * @par Shared memory
 * Three STATIC shared arrays: w[ZSIZE], data0[BLOCKS_PER_GROUP*ZSIZE],
 * data1[BLOCKS_PER_GROUP*ZSIZE]. Total = (ZSIZE + 2*BLOCKS_PER_GROUP*ZSIZE) *
 * sizeof(real2) bytes per block (sizeof(real2) = 8 single / 16 double). This is the
 * dominant occupancy limiter; BLOCKS_PER_GROUP trades block width against per-block
 * shared use. (The OpenCL twin passes these three buffers as dynamic __local args of
 * the same sizes; here they are declared static, so no dynamic-shared launch argument
 * is needed.) data0 and data1 alias no global memory; the twiddle table is rebuilt per
 * launch, not cached across passes.
 *
 * @par Memory access and transpose (contract)
 * Read: element (x,y,z) of a row comes from the contiguous slot
 * in[x*(YSIZE*ZSIZE) + y*ZSIZE + z] (transform axis z is innermost/coalesced), or via
 * loadComplexValue when INPUT_IS_PACKED. Write (emitted by COMPUTE_FFT): the
 * transformed row lands at out[y*(ZSIZE*XSIZE) + z*XSIZE + x] for the complex case, or
 * out[y*(ZSIZE*(XSIZE/2+1)) + z*(XSIZE/2+1) + x] when OUTPUT_IS_PACKED. This is a
 * deliberate transpose: the just-transformed axis (local Z) moves to an outer stride
 * and local X becomes innermost, so the NEXT pass reads its own axis contiguously. The
 * store is therefore strided in global memory (stride XSIZE or XSIZE/2+1), which is the
 * primary coalescing cost of this kernel. Each output element is written exactly once;
 * no atomics, so the result is bitwise deterministic run-to-run for a fixed launch
 * shape.
 *
 * @par Synchronization and participation
 * Two block-wide __syncthreads() barriers per grid-stride iteration path: after filling
 * w[] and after staging the row, before COMPUTE_FFT reads across the whole row. All
 * threads of the block must reach both barriers, so the OUTPUT_IS_PACKED row-skip and
 * the index<XSIZE*YSIZE tail guard gate only the loads, never the barriers. COMPUTE_FFT
 * contributes one further __syncthreads() per radix pass. No cross-block synchronization
 * and no warp-primitive participation assumptions.
 *
 * @param[in] in  input array, element type INPUT_TYPE (real when INPUT_IS_REAL, else
 *                real2); borrowed, read-only device buffer. Logical layout
 *                in[x*(YSIZE*ZSIZE) + y*ZSIZE + z], or the packed half-complex grid
 *                (z-stride ZSIZE/2+1) when INPUT_IS_PACKED. May not alias @p out.
 * @param[out] out output array, element type OUTPUT_TYPE; borrowed device buffer written
 *                exactly once per element by COMPUTE_FFT in the transposed (and possibly
 *                Hermitian-packed or real-narrowed) layout above. Distinct from @p in;
 *                the driver ping-pongs the two grid buffers between passes.
 *
 * @pre blockDim.x == BLOCKS_PER_GROUP*THREADS_PER_BLOCK.
 * @pre ZSIZE factors entirely into {2,3,4,5,7}; otherwise createKernel refuses to emit
 *      COMPUTE_FFT (enforced host-side, e.g. via findLegalDimension).
 * @pre Total static shared bytes (see @par Shared memory) fit the target SM budget; the
 *      OpenCL driver retries with a smaller block when the device rejects the size.
 * @pre @p in and @p out are distinct, correctly sized device buffers for this pass's
 *      element types and layouts.
 * @post @p out holds the transposed 1D transform of every row; @p in is unmodified.
 *
 * @note Unnormalized: forward-then-inverse over all three axes scales data by
 *       XSIZE*YSIZE*ZSIZE (the point count).
 */

extern "C" __global__ void execFFT(const INPUT_TYPE* __restrict__ in, OUTPUT_TYPE* __restrict__ out) {
    // Non-local: in the live OpenCL twin (OpenCLFFT3D::createKernel) these three are
    // DYNAMIC __local kernel args, sized at launch by the driver's setArg(2/3/4, ...).
    // Here they are STATIC __shared__, so ZSIZE and BLOCKS_PER_GROUP must be compile-time
    // constants and no dynamic-shared-memory launch argument is passed. w = twiddle
    // table; data0/data1 = ping-pong row buffers COMPUTE_FFT alternates between per pass.
    __shared__ real2 w[ZSIZE];
    __shared__ real2 data0[BLOCKS_PER_GROUP*ZSIZE];
    __shared__ real2 data1[BLOCKS_PER_GROUP*ZSIZE];
    // Cooperatively fill the twiddle table w[k] = e^{-i*SIGN*2*pi*k/ZSIZE}
    // (SIGN = +1 forward, -1 inverse). Indexed inside COMPUTE_FFT as w[j*ZSIZE/(r*L)].
    for (int i = threadIdx.x; i < ZSIZE; i += blockDim.x)
        w[i] = make_real2(cos(-(SIGN)*i*2*M_PI/ZSIZE), sin(-(SIGN)*i*2*M_PI/ZSIZE));
    __syncthreads();
    
    // Non-local / ASSUMED: THREADS_PER_BLOCK has no injector in this tree (the OpenCL twin
    // uses get_local_size(0)/ZSIZE slicing; the CUDA driver that once set it was dropped
    // with the cuFFT switch). Inferred value = ZSIZE, so blockDim.x = BLOCKS_PER_GROUP*ZSIZE,
    // matching OpenCL's threads = blocksPerGroup*zsize. Sub-block id selects this thread's row.
    const int block = threadIdx.x/THREADS_PER_BLOCK;
    // Grid-stride loop over all XSIZE*YSIZE rows; each iteration this block handles the
    // BLOCKS_PER_GROUP rows starting at baseIndex (this thread's row = baseIndex+block).
    for (int baseIndex = blockIdx.x*BLOCKS_PER_GROUP; baseIndex < XSIZE*YSIZE; baseIndex += gridDim.x*BLOCKS_PER_GROUP) {
        int index = baseIndex+block;
        // Decode the flat row index into (x,y); the transform runs along the remaining z.
        int x = index/YSIZE;
        int y = index-x*YSIZE;
#if OUTPUT_IS_PACKED
        // Hermitian symmetry: only rows with x in the lower half carry independent data,
        // so skip loading/transforming the redundant upper-half rows entirely.
        if (x < XSIZE/2+1) {
#endif
        // Stage row z-line into shared memory (data0), THREADS_PER_BLOCK threads striding
        // over its ZSIZE elements. The guard drops the padding row of the last partial tile.
        if (index < XSIZE*YSIZE)
            for (int i = threadIdx.x-block*THREADS_PER_BLOCK; i < ZSIZE; i += THREADS_PER_BLOCK)
    #if INPUT_IS_REAL
                // Real input: widen each sample to a complex value with zero imaginary part.
                data0[i+block*ZSIZE] = make_real2(in[x*(YSIZE*ZSIZE)+y*ZSIZE+i], 0);
    #elif INPUT_IS_PACKED
                // Packed half-complex input: reconstruct the full spectrum via Hermitian symmetry.
                data0[i+block*ZSIZE] = loadComplexValue(in, x, y, i);
    #else
                // Plain complex input: direct contiguous copy from the row-major grid.
                data0[i+block*ZSIZE] = in[x*(YSIZE*ZSIZE)+y*ZSIZE+i];
    #endif
#if OUTPUT_IS_PACKED
        }
#endif
        // Ensure the whole row is staged before the butterflies read across it.
        __syncthreads();
        // Non-local: this expands to the straight-line code string generated by the factor
        // loop in OpenCLFFT3D::createKernel (radix 2/3/4/5/7 butterflies, ping-ponging
        // data0<->data1, twiddled via multiplyComplex against w[]). Its final statement is
        // the ONLY global write of this kernel: a transposed store to out at stride XSIZE
        // (or XSIZE/2+1 when OUTPUT_IS_PACKED) -- uncoalesced, unlike the coalesced staging
        // read above, and thus the primary global-bandwidth cost of the pass.
        COMPUTE_FFT
    }
}
