/** ==========================================================================
 * @file fftR2C.cu
 * @brief Hermitian packing/unpacking kernels wrapped around a full complex 3D
 *        FFT so that a real 3D transform (R2C) and its inverse (C2R) can be
 *        computed with a complex transform of HALF the data.
 *
 * ==========================================================================
 * 0. PROVENANCE AND ORPHANED STATUS - READ FIRST
 * ==========================================================================
 *
 * This translation unit is NOT part of the CUDA execution path. In this
 * vendored tree the CUDA 3D FFT is class CudaFFT3D, which is implemented
 * entirely on top of cuFFT (cufftPlan3d / cufftExecR2C / cufftExecC2R and the
 * Z2Z/C2C variants); it never builds a module from these kernels and never
 * launches packForwardData, unpackForwardData, packBackwardData, or
 * unpackBackwardData. The source text of this file is embedded as the string
 * CudaKernelSources::fftR2C, but that symbol has no consumer anywhere in the
 * CUDA platform - it is dead. Treat this file as the CUDA-syntax twin of the
 * OpenCL kernel fftR2C.cl and nothing more.
 *
 * The authoritative, live implementation is class OpenCLFFT3D (basename
 * OpenCLFFT3D.cpp), specifically its non-VkFFT branch. Every compile-time macro
 * below is computed by the OpenCLFFT3D constructor, the kernels are launched by
 * OpenCLFFT3D::execFFT, and the surrounding butterfly-FFT kernels are generated
 * by OpenCLFFT3D::createKernel. All contracts stated here are recovered from
 * those three symbols; there is no CUDA launch site to recover them from.
 *
 * ==========================================================================
 * 1. What these kernels do (the "pack a real transform into a complex one" trick)
 * ==========================================================================
 *
 * A real N-point DFT has conjugate (Hermitian) symmetry: X[k] = conj(X[N-k]).
 * So an N-point real signal has only ~N/2+1 independent complex outputs, and an
 * N/2-point complex FFT carries exactly enough information to reconstruct it.
 * The device is: pick one EVEN axis of the real grid (the packed axis),
 * interleave its adjacent real samples into the real/imag lanes of a complex
 * value (h[n] = x[2n] + i*x[2n+1]), run a full complex 3D FFT on the resulting
 * half-size complex grid, then algebraically split the even/odd sub-transforms
 * back apart with twiddle factors. The forward and inverse pipelines, as
 * sequenced by OpenCLFFT3D::execFFT:
 *
 *   R2C (forward): packForwardData -> complex 3D FFT (3 axis passes) -> unpackForwardData
 *   C2R (inverse): packBackwardData -> inverse complex 3D FFT (3 axis passes) -> unpackBackwardData
 *
 * These four kernels are ONLY the pre/post packing steps; the butterfly FFT
 * itself is the separately generated execFFT kernel (basename fft.cu / fft.cl).
 *
 * ==========================================================================
 * 2. Two DISTINCT reductions - do not conflate them
 * ==========================================================================
 *
 *  1. The PACKED-AXIS reduction (internal, computational): real pairs along the
 *     chosen even axis are folded into complex lanes, halving THAT axis. Which
 *     axis is chosen is PACKED_AXIS; its halved length is the corresponding
 *     PACKED_*SIZE. This is invisible to callers of the FFT.
 *  2. The HERMITIAN OUTPUT reduction (externally visible R2C storage format):
 *     the emitted half-complex spectrum is always stored reduced along Z (the
 *     inner/last axis) to length ZSIZE/2+1, with full XSIZE x YSIZE extent
 *     (outputZSize / inputZSize = ZSIZE/2+1). This is the standard R2C layout,
 *     produced solely by unpackForwardData and consumed solely by
 *     loadComplexValue, and is INDEPENDENT of which axis was packed.
 *
 * ==========================================================================
 * 3. Grid index conventions (row-major, x outer / y middle / z inner)
 * ==========================================================================
 *
 *   Real grid            (XSIZE x YSIZE x ZSIZE):          idx = x*YSIZE*ZSIZE + y*ZSIZE + z
 *   Packed complex grid  (PACKED_X x PACKED_Y x PACKED_Z): idx = x*PACKED_YSIZE*PACKED_ZSIZE + y*PACKED_ZSIZE + z
 *   Half-complex grid    (XSIZE x YSIZE x (ZSIZE/2+1)):    idx = x*YSIZE*outputZSize + y*outputZSize + z
 *
 * ==========================================================================
 * 4. Launch contract (recovered from OpenCLFFT3D::execFFT)
 * ==========================================================================
 *
 * All four kernels are launched with a work-item / thread count of exactly
 * gridSize = XSIZE*YSIZE*ZSIZE/2, which equals the product of the three
 * PACKED_*SIZE values (one axis is halved, so the packed grid holds half the
 * real cells). Each kernel is a 1D grid-stride loop over that gridSize, so the
 * grid and block shape are otherwise unconstrained: any block size that covers
 * up to gridSize is legal, and correctness does not depend on it.
 *
 * Argument order for every kernel is (input buffer, output buffer). The driver
 * ping-pongs two device buffers (the caller's in and out arrays): the pack
 * kernel reads the caller buffer and writes the scratch buffer, the three FFT
 * passes alternate between the two, and the unpack kernel reads the FFT result
 * and writes the caller's destination buffer. All buffers are device global
 * memory owned by the caller (the FFT object borrows them for the duration of
 * execFFT); these kernels neither allocate nor free. Completion is ordered by
 * the platform command queue / stream; results are valid only after that
 * queue/stream is synchronized. No kernel here synchronizes across blocks and
 * none is a cooperative launch.
 *
 * packForwardData and unpackBackwardData use NO shared memory. unpackForwardData
 * and packBackwardData each declare a per-block __shared__ twiddle table `w` of
 * PACKED_<PACKED_AXIS>SIZE complex entries (PACKED_<axis>SIZE * sizeof(real2)
 * bytes, fixed at compile time). Every thread of a block cooperatively fills its
 * block's copy of `w`, then __syncthreads() before any read; consequently all
 * threads of a block must reach that barrier (no divergence may skip it). The
 * table is indexed by the thread's coordinate along the packed axis, whose range
 * is exactly [0, PACKED_<axis>SIZE), so the static size always covers every read.
 *
 * NOTE (CUDA vs OpenCL difference): here `w` is a statically sized __shared__
 * array. The OpenCL twin instead receives `w` as a __local kernel argument
 * (arg index 2) whose byte size the OpenCLFFT3D constructor sets at build time
 * to bufferSize*sizeof(real2), bufferSize = PACKED_<axis>SIZE; the two are the
 * same allocation expressed differently.
 *
 * ==========================================================================
 * 5. Compile-time parameters (injected as preprocessor #defines before JIT)
 * ==========================================================================
 *
 * These macros are NOT #defined in this translation unit; a host driver builds
 * the module with a defines map, so an editor/clangd cannot resolve them. The
 * values below are exactly those computed by the OpenCLFFT3D constructor.
 *
 *  @param XSIZE, YSIZE, ZSIZE  (int) Full real-grid dimensions along x, y, z,
 *        substituted verbatim from the xsize/ysize/zsize constructor arguments.
 *        These are the caller's (e.g. PME) real charge-grid dimensions.
 *
 *  @param PACKED_AXIS  (int, 0|1|2) The axis chosen for the real->complex
 *        packing: 0=x, 1=y, 2=z. The OpenCLFFT3D constructor picks the FIRST
 *        even axis in x,y,z order (xsize%2==0 -> 0; else ysize%2==0 -> 1; else
 *        zsize%2==0 -> 2). If NO axis is even, packing is disabled
 *        (packRealAsComplex=false) and these four kernels are never built or
 *        launched; a plain complex FFT with a real input path is used instead.
 *        Hence whenever this file is compiled, PACKED_AXIS names an axis whose
 *        full size is even, and that evenness is a standing precondition of
 *        packForwardData and unpackBackwardData.
 *
 *  @param PACKED_XSIZE, PACKED_YSIZE, PACKED_ZSIZE  (int) Dimensions of the
 *        packed complex grid. Exactly ONE of them - the one selected by
 *        PACKED_AXIS - equals the full size / 2; the other two equal their full
 *        sizes unchanged. Concretely:
 *          PACKED_AXIS==0: (XSIZE/2, YSIZE,   ZSIZE)
 *          PACKED_AXIS==1: (XSIZE,   YSIZE/2, ZSIZE)
 *          PACKED_AXIS==2: (XSIZE,   YSIZE,   ZSIZE/2)
 *        So X/Y/ZSIZE are the full real extents; PACKED_*SIZE is the reduced
 *        extent along whichever single axis was packed.
 *
 *  @param M_PI  (double literal) Value of pi, injected as a numeric literal
 *        (doubleToString(M_PI)); used to build the twiddle factors.
 *
 * NOTE: real, real2, make_real2 are not macros defined here - they come from the
 * platform's numeric prelude and resolve to float/float2 (single precision) or
 * double/double2 (double precision) per the context's precision mode. A real2 is
 * 8 or 16 bytes accordingly; the __shared__ `w` size scales with it.
 * ==========================================================================
 */

/**
 * @brief R2C forward, step 1 of 3: interleave adjacent real samples along the
 *        packed axis into the real/imag lanes of a half-size complex grid.
 *
 * Forms h = x[even] + i*x[odd] for every packed cell, where even/odd are the two
 * real samples at positions 2k and 2k+1 along PACKED_AXIS. This is the step that
 * folds a real grid into a complex grid half as large so a single complex FFT
 * can carry the real transform. Writes each output cell exactly once; no
 * symmetry is elided here (the Hermitian reduction happens in unpackForwardData).
 *
 * @param[in]  in   Real input grid, XSIZE x YSIZE x ZSIZE, row-major, device
 *                  global memory; borrowed, read-only.
 * @param[out] out  Packed complex grid, PACKED_XSIZE x PACKED_YSIZE x
 *                  PACKED_ZSIZE, device global memory; borrowed. Each cell is
 *                  written once as (x[2k], x[2k+1]) with the pair taken along
 *                  PACKED_AXIS.
 *
 * @pre  The full size of the packed axis is even (guaranteed by the OpenCLFFT3D
 *       PACKED_AXIS selection); otherwise this kernel is not built.
 * @pre  @p in and @p out do not alias (distinct ping-pong buffers).
 * @par Launch configuration:
 *      1D grid-stride loop over gridSize = product of PACKED_*SIZE; grid/block
 *      shape otherwise unconstrained, no shared memory, no synchronization.
 */
extern "C" __global__ void packForwardData(const real* __restrict__ in, real2* __restrict__ out) {
    const int gridSize = PACKED_XSIZE*PACKED_YSIZE*PACKED_ZSIZE;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        int x = index/(PACKED_YSIZE*PACKED_ZSIZE);
        int remainder = index-x*(PACKED_YSIZE*PACKED_ZSIZE);
        int y = remainder/PACKED_ZSIZE;
        int z = remainder-y*PACKED_ZSIZE;
        // Only the branch selected by PACKED_AXIS is compiled. The 2*coord+1 read
        // cannot overrun because the OpenCLFFT3D constructor picks PACKED_AXIS only
        // when that axis's full size is even (else these kernels are not built).
#if PACKED_AXIS == 0
        real2 value = make_real2(in[2*x*YSIZE*ZSIZE+y*ZSIZE+z], in[(2*x+1)*YSIZE*ZSIZE+y*ZSIZE+z]);
#elif PACKED_AXIS == 1
        real2 value = make_real2(in[x*YSIZE*ZSIZE+2*y*ZSIZE+z], in[x*YSIZE*ZSIZE+(2*y+1)*ZSIZE+z]);
#else
        real2 value = make_real2(in[x*YSIZE*ZSIZE+y*ZSIZE+2*z], in[x*YSIZE*ZSIZE+y*ZSIZE+(2*z+1)]);
#endif
        out[index] = value;
    }
}

/**
 * @brief R2C forward, step 3 of 3: split the complex FFT of the packed grid back
 *        into the true half-complex spectrum of the original real grid.
 *
 * After a full complex 3D FFT of the packed grid, each output blends the DFTs of
 * the even- and odd-indexed real samples along the packed axis. This kernel
 * un-blends them: for packed point z1 = H(x,y,z) and its packed-grid Hermitian
 * partner z2 = H(xp,yp,zp) (reflection across all three packed axes with the
 * 0->0 fixed point), it applies the split-radix recombination
 *   X = ( (z1 + conj(z2)) - i*w*(z1 - conj(z2)) ) / 2,
 * with twiddle w = (sin theta, cos theta), theta = 2*pi*k/N, N = full length of
 * the packed axis, k = the packed-axis coordinate. It emits both the spectral
 * point and, for the folded rows, its conjugate mirror, so the output grid is
 * fully populated including the packed-axis Nyquist plane. This kernel alone
 * produces the externally visible Z-reduced R2C layout.
 *
 * @param[in]  in   Packed complex grid, PACKED_XSIZE x PACKED_YSIZE x
 *                  PACKED_ZSIZE, device global memory - the complex FFT of
 *                  packForwardData's result; borrowed, read-only.
 * @param[out] out  Half-complex spectrum, XSIZE x YSIZE x (ZSIZE/2+1), Z-reduced
 *                  (the externally visible R2C layout), device global memory;
 *                  borrowed.
 *
 * @pre  @p in and @p out do not alias.
 * @par Shared memory:
 *      Per-block __shared__ real2 w[PACKED_<PACKED_AXIS>SIZE], filled
 *      cooperatively as w[i] = (sin(2*pi*i/N), cos(2*pi*i/N)); all threads of the
 *      block must reach the __syncthreads() that follows before any read of `w`.
 * @par Launch configuration:
 *      1D grid-stride loop over gridSize = product of PACKED_*SIZE; grid/block
 *      shape otherwise unconstrained.
 * @note Write pattern: for the direct point only z < ZSIZE/2+1 is stored; the
 *       redundant upper-Z half is elided. A second write emits the
 *       conjugate-symmetric partner at the full-grid reflection (XSIZE-x,
 *       YSIZE-y, ZSIZE-z), with a dedicated branch writing the packed-axis
 *       Nyquist row (index PACKED_<axis>SIZE) when the packed coordinate is 0.
 */
extern "C" __global__ void unpackForwardData(const real2* __restrict__ in, real2* __restrict__ out) {
    // Compute the phase factors (twiddles). w[i] = (sin, cos) of 2*pi*i/N, where
    // N is the FULL length of the packed axis (XSIZE/YSIZE/ZSIZE per PACKED_AXIS).
    // CUDA-only shape: `w` is per-block __shared__, statically sized here. The
    // live OpenCL twin (OpenCLFFT3D constructor / fftR2C.cl) instead passes `w`
    // as a __local kernel argument sized at build time - same table, and it is
    // rebuilt independently by every block, so shrinking the block does not
    // reduce this per-block traffic.
#if PACKED_AXIS == 0
    __shared__ real2 w[PACKED_XSIZE];
    for (int i = threadIdx.x; i < PACKED_XSIZE; i += blockDim.x)
        w[i] = make_real2(sin(i*2*M_PI/XSIZE), cos(i*2*M_PI/XSIZE));
#elif PACKED_AXIS == 1
    __shared__ real2 w[PACKED_YSIZE];
    for (int i = threadIdx.x; i < PACKED_YSIZE; i += blockDim.x)
        w[i] = make_real2(sin(i*2*M_PI/YSIZE), cos(i*2*M_PI/YSIZE));
#else
    __shared__ real2 w[PACKED_ZSIZE];
    for (int i = threadIdx.x; i < PACKED_ZSIZE; i += blockDim.x)
        w[i] = make_real2(sin(i*2*M_PI/ZSIZE), cos(i*2*M_PI/ZSIZE));
#endif
    __syncthreads();

    // Transform the data.
    
    const int gridSize = PACKED_XSIZE*PACKED_YSIZE*PACKED_ZSIZE;
    const int outputZSize = ZSIZE/2+1;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        int x = index/(PACKED_YSIZE*PACKED_ZSIZE);
        int remainder = index-x*(PACKED_YSIZE*PACKED_ZSIZE);
        int y = remainder/PACKED_ZSIZE;
        int z = remainder-y*PACKED_ZSIZE;
        int xp = (x == 0 ? 0 : PACKED_XSIZE-x);
        int yp = (y == 0 ? 0 : PACKED_YSIZE-y);
        int zp = (z == 0 ? 0 : PACKED_ZSIZE-z);
        real2 z1 = in[x*PACKED_YSIZE*PACKED_ZSIZE+y*PACKED_ZSIZE+z];
        real2 z2 = in[xp*PACKED_YSIZE*PACKED_ZSIZE+yp*PACKED_ZSIZE+zp];
#if PACKED_AXIS == 0
        real2 wfac = w[x];
#elif PACKED_AXIS == 1
        real2 wfac = w[y];
#else
        real2 wfac = w[z];
#endif
        real2 output = make_real2((z1.x+z2.x - wfac.x*(z1.x-z2.x) + wfac.y*(z1.y+z2.y))/2, (z1.y-z2.y - wfac.y*(z1.x-z2.x) - wfac.x*(z1.y+z2.y))/2);
        if (z < outputZSize)
            out[x*YSIZE*outputZSize+y*outputZSize+z] = output;
        xp = (x == 0 ? 0 : XSIZE-x);
        yp = (y == 0 ? 0 : YSIZE-y);
        zp = (z == 0 ? 0 : ZSIZE-z);
        // Nyquist plane of the packed axis: PACKED_<axis>SIZE equals full/2, the
        // one full-grid row along the packed axis that the direct write above can
        // never reach (its packed coordinate folds to 0). loadComplexValue on the
        // C2R side reads this row back, so it must be populated here.
        if (zp < outputZSize) {
#if PACKED_AXIS == 0
            if (x == 0)
                out[PACKED_XSIZE*YSIZE*outputZSize+yp*outputZSize+zp] = make_real2((z1.x-z1.y+z2.x-z2.y)/2, (-z1.x-z1.y+z2.x+z2.y)/2);
#elif PACKED_AXIS == 1
            if (y == 0)
                out[xp*YSIZE*outputZSize+PACKED_YSIZE*outputZSize+zp] = make_real2((z1.x-z1.y+z2.x-z2.y)/2, (-z1.x-z1.y+z2.x+z2.y)/2);
#else
            if (z == 0)
                out[xp*YSIZE*outputZSize+yp*outputZSize+PACKED_ZSIZE] = make_real2((z1.x-z1.y+z2.x-z2.y)/2, (-z1.x-z1.y+z2.x+z2.y)/2);
#endif
            else
                out[xp*YSIZE*outputZSize+yp*outputZSize+zp] = make_real2(output.x, -output.y);
        }
    }
}

/**
 * @brief Fetch spectrum value (x,y,z) from the Z-reduced half-complex grid,
 *        reconstructing the elided upper-Z half on the fly via Hermitian symmetry.
 *
 * The half-complex grid stores only z in [0, ZSIZE/2] (inputZSize = ZSIZE/2+1).
 * For a requested z at or below that it reads directly; for z in the missing
 * upper half it returns conj(X[XSIZE-x, YSIZE-y, ZSIZE-z]) (reflected in-range
 * read with negated imaginary part). Lets the C2R path address the spectrum as
 * if it were full size. Device-side, called from packBackwardData; contains no
 * synchronization and no shared-memory access, so it imposes no divergence
 * constraint on its caller.
 *
 * @param[in] in     Half-complex grid, XSIZE x YSIZE x (ZSIZE/2+1), device
 *                   global memory; borrowed, read-only.
 * @param[in] x,y,z  Logical coordinates into the FULL (unreduced) spectrum.
 * @return The complex spectral value, conjugated when reconstructed from the
 *         reflected point.
 * @note Reflection uses the 0->0 fixed point on x and y; z is reflected as
 *       ZSIZE-z only in the branch where z >= inputZSize.
 */
static __inline__ __device__ real2 loadComplexValue(const real2* __restrict__ in, int x, int y, int z) {
    const int inputZSize = ZSIZE/2+1;
    if (z < inputZSize)
        return in[x*YSIZE*inputZSize+y*inputZSize+z];
    // z >= inputZSize is the upper-Z half that unpackForwardData deliberately
    // never stored (its `z < outputZSize` guard). Reconstruct it here from the
    // stored conjugate mirror. The elision is always along Z (inputZSize =
    // ZSIZE/2+1) regardless of which axis PACKED_AXIS packed.
    int xp = (x == 0 ? 0 : XSIZE-x);
    int yp = (y == 0 ? 0 : YSIZE-y);
    real2 value = in[xp*YSIZE*inputZSize+yp*inputZSize+(ZSIZE-z)];
    return make_real2(value.x, -value.y);
}

/**
 * @brief C2R inverse, step 1 of 3: fold the half-complex spectrum into a
 *        half-size complex grid whose inverse complex FFT yields the packed real
 *        pairs. Algebraic inverse of unpackForwardData.
 *
 * Reads z1 = X(x,y,z) and its packed-axis-reflected partner z2 (packed
 * coordinate reflected to PACKED_<axis>SIZE - coord, other axes mirrored), both
 * via loadComplexValue so the elided upper-Z half is reconstructed
 * transparently. Splits into even = (z1 + conj z2)/2 and odd = (z1 - conj z2)/2
 * rotated by twiddle wfac, then writes (even.x - odd.y, even.y + odd.x) once per
 * packed cell. The result feeds the inverse complex 3D FFT.
 *
 * @param[in]  in   Half-complex spectrum, XSIZE x YSIZE x (ZSIZE/2+1) (the R2C
 *                  layout), device global memory; borrowed, read-only.
 * @param[out] out  Packed complex grid, PACKED_XSIZE x PACKED_YSIZE x
 *                  PACKED_ZSIZE, device global memory; borrowed.
 *
 * @pre  @p in and @p out do not alias.
 * @par Shared memory:
 *      Per-block __shared__ real2 w[PACKED_<PACKED_AXIS>SIZE], filled as
 *      w[i] = (cos(2*pi*i/N), sin(2*pi*i/N)) - the inverse-rotation conjugate of
 *      the (sin, cos) table in unpackForwardData; all threads of the block must
 *      reach the following __syncthreads() before any read of `w`.
 * @par Launch configuration:
 *      1D grid-stride loop over gridSize = product of PACKED_*SIZE; grid/block
 *      shape otherwise unconstrained.
 * @note The /2 factors here are undone by the x2 scaling in unpackBackwardData.
 */
extern "C" __global__ void packBackwardData(const real2* __restrict__ in, real2* __restrict__ out) {
    // Compute the phase factors (twiddles). w[i] = (cos, sin) of 2*pi*i/N, the
    // inverse-rotation conjugate of the (sin, cos) table used on the forward side.

#if PACKED_AXIS == 0
    __shared__ real2 w[PACKED_XSIZE];
    for (int i = threadIdx.x; i < PACKED_XSIZE; i += blockDim.x)
        w[i] = make_real2(cos(i*2*M_PI/XSIZE), sin(i*2*M_PI/XSIZE));
#elif PACKED_AXIS == 1
    __shared__ real2 w[PACKED_YSIZE];
    for (int i = threadIdx.x; i < PACKED_YSIZE; i += blockDim.x)
        w[i] = make_real2(cos(i*2*M_PI/YSIZE), sin(i*2*M_PI/YSIZE));
#else
    __shared__ real2 w[PACKED_ZSIZE];
    for (int i = threadIdx.x; i < PACKED_ZSIZE; i += blockDim.x)
        w[i] = make_real2(cos(i*2*M_PI/ZSIZE), sin(i*2*M_PI/ZSIZE));
#endif
    __syncthreads();

    // Transform the data.
    
    const int gridSize = PACKED_XSIZE*PACKED_YSIZE*PACKED_ZSIZE;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        int x = index/(PACKED_YSIZE*PACKED_ZSIZE);
        int remainder = index-x*(PACKED_YSIZE*PACKED_ZSIZE);
        int y = remainder/PACKED_ZSIZE;
        int z = remainder-y*PACKED_ZSIZE;
        int xp = (x == 0 ? 0 : PACKED_XSIZE-x);
        int yp = (y == 0 ? 0 : PACKED_YSIZE-y);
        int zp = (z == 0 ? 0 : PACKED_ZSIZE-z);
        real2 z1 = loadComplexValue(in, x, y, z);
#if PACKED_AXIS == 0
        real2 wfac = w[x];
        real2 z2 = loadComplexValue(in, PACKED_XSIZE-x, yp, zp);
#elif PACKED_AXIS == 1
        real2 wfac = w[y];
        real2 z2 = loadComplexValue(in, xp, PACKED_YSIZE-y, zp);
#else
        real2 wfac = w[z];
        real2 z2 = loadComplexValue(in, xp, yp, PACKED_ZSIZE-z);
#endif
        real2 even = make_real2((z1.x+z2.x)/2, (z1.y-z2.y)/2);
        real2 odd = make_real2((z1.x-z2.x)/2, (z1.y+z2.y)/2);
        odd = make_real2(odd.x*wfac.x-odd.y*wfac.y, odd.y*wfac.x+odd.x*wfac.y);
        out[x*PACKED_YSIZE*PACKED_ZSIZE+y*PACKED_ZSIZE+z] = make_real2(even.x-odd.y, even.y+odd.x);
    }
}

/**
 * @brief C2R inverse, step 3 of 3: unpack the inverse complex FFT of the packed
 *        grid into a full-size real grid. Inverse of packForwardData.
 *
 * Each packed complex value carries two adjacent real samples along the packed
 * axis in its real/imag lanes; scatters value.x -> the even sample (2k) and
 * value.y -> the odd sample (2k+1) along PACKED_AXIS, each written exactly once.
 * The x2 factor restores the normalization removed by the /2 folds in
 * packBackwardData / unpackForwardData. Pure de-interleave: no symmetry logic
 * (the Hermitian reconstruction already happened in loadComplexValue /
 * packBackwardData).
 *
 * @param[in]  in   Packed complex grid, PACKED_XSIZE x PACKED_YSIZE x
 *                  PACKED_ZSIZE - output of the inverse complex FFT, device
 *                  global memory; borrowed, read-only.
 * @param[out] out  Full real grid, XSIZE x YSIZE x ZSIZE, row-major, device
 *                  global memory; borrowed.
 *
 * @pre  The full size of the packed axis is even (guaranteed by the OpenCLFFT3D
 *       PACKED_AXIS selection).
 * @pre  @p in and @p out do not alias.
 * @par Launch configuration:
 *      1D grid-stride loop over gridSize = product of PACKED_*SIZE; grid/block
 *      shape otherwise unconstrained, no shared memory, no synchronization.
 */
extern "C" __global__ void unpackBackwardData(const real2* __restrict__ in, real* __restrict__ out) {
    const int gridSize = PACKED_XSIZE*PACKED_YSIZE*PACKED_ZSIZE;
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < gridSize; index += blockDim.x*gridDim.x) {
        int x = index/(PACKED_YSIZE*PACKED_ZSIZE);
        int remainder = index-x*(PACKED_YSIZE*PACKED_ZSIZE);
        int y = remainder/PACKED_ZSIZE;
        int z = remainder-y*PACKED_ZSIZE;
        real2 value = 2*in[index];
#if PACKED_AXIS == 0
        out[2*x*YSIZE*ZSIZE+y*ZSIZE+z] = value.x;
        out[(2*x+1)*YSIZE*ZSIZE+y*ZSIZE+z] = value.y;
#elif PACKED_AXIS == 1
        out[x*YSIZE*ZSIZE+2*y*ZSIZE+z] = value.x;
        out[x*YSIZE*ZSIZE+(2*y+1)*ZSIZE+z] = value.y;
#else
        out[x*YSIZE*ZSIZE+y*ZSIZE+2*z] = value.x;
        out[x*YSIZE*ZSIZE+y*ZSIZE+(2*z+1)] = value.y;
#endif
    }
}
