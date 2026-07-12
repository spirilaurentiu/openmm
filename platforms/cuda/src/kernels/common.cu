/**
 * @file common.cu
 * @brief CUDA backing of OpenMM's portable "common compute" abstraction layer.
 *
 * OpenMM's `common` framework lets a *single* kernel source file be compiled
 * unchanged for three GPU/CPU backends: CUDA (this file), OpenCL
 * (`platforms/opencl/src/kernels/common.cl`), and the CPU/C++ reference
 * (`common.cc`). Portability is achieved by writing kernels against a fixed set
 * of UPPERCASE macros (KERNEL, DEVICE, LOCAL, GLOBAL_ID, SYNC_THREADS,
 * ATOMIC_ADD, ...) and portable typedefs (mm_long, real, mixed). Each backend
 * supplies one small prologue file that expands those macros to the constructs
 * native to its language. This file is that prologue for CUDA; every macro below
 * maps a portable name onto a CUDA/PTX construct.
 *
 * ------------------------------------------------------------------------------
 * How this file reaches the compiler (see CudaContext::createModule in
 * platforms/cuda/src/CudaContext.cpp, ~line 472):
 *
 * Every kernel module NVRTC compiles at runtime is assembled, in order, as:
 *   1. `#define`s from CudaContext::compilationDefines (precision-, PBC-, and
 *      warp-intrinsic-dependent; see "externally injected macros" below).
 *   2. `typedef ... real/real2..; typedef ... mixed..; typedef unsigned int tileflags;`
 *      emitted directly by createModule.
 *   3. THIS FILE, injected verbatim as `CudaKernelSources::common`.
 *   4. Per-module `#define`s passed to createModule.
 *   5. The actual kernel source.
 *
 * Consequence: `real` (used by realToFixedPoint below) and the precision macros
 * are ALREADY in scope by the time this file is textually inserted — this file
 * does not define `real`; step 2 does.
 *
 * ------------------------------------------------------------------------------
 * Comments have ZERO runtime cost: openmm/cmake_modules/strip_comments.py strips
 * all comments before the source is string-encoded into CudaKernelSources, so
 * NVRTC never sees them. Document freely.
 *
 * ------------------------------------------------------------------------------
 * Macros DEFINED here vs. macros INJECTED by the C++ layer.
 *
 * DEFINED HERE (this file is the source of truth for the portability layer):
 *   KERNEL, DEVICE, LOCAL, LOCAL_ARG, GLOBAL, RESTRICT, LOCAL_ID, LOCAL_SIZE,
 *   GLOBAL_ID, GLOBAL_SIZE, GROUP_ID, NUM_GROUPS, SYNC_THREADS, MEM_FENCE,
 *   ATOMIC_ADD, FLT_MAX, SUPPORTS_64_BIT_ATOMICS, SUPPORTS_DOUBLE_PRECISION,
 *   mm_long, mm_ulong, realToFixedPoint.
 *
 * INJECTED by CudaContext.cpp (NOT defined here; this file relies on them being
 * present above it). Traced to CudaContext::CudaContext / createModule:
 *   - `real` / `real2..4`, `mixed` / `mixed2..4`, `tileflags`: emitted as
 *     `typedef`s by createModule (double vs float chosen by precision mode).
 *   - USE_DOUBLE_PRECISION (= "1"): set only when useDoublePrecision
 *     (CudaContext.cpp ~line 227). USE_MIXED_PRECISION (= "1") only in mixed mode
 *     (~line 239). In single precision NEITHER is defined. Downstream kernels
 *     `#ifdef` on these; this prologue itself does not, but they gate the typedef
 *     of `real` that realToFixedPoint consumes.
 *   - make_real2/3/4, make_mixed2/3/4 (~lines 228-255): expand to the
 *     make_double2/3/4 or make_float2/3/4 constructors per precision.
 *   - SQRT/RSQRT/EXP/... math-function names (~lines 279-294): double vs float
 *     variants.
 *   - SYNC_WARPS, SHFL, BALLOT (~lines 215-222): warp intrinsics, sync vs
 *     non-sync variants chosen by CUDA driver version (>= 9000 => *_sync).
 *   - APPLY_PERIODIC_TO_DELTA / _POS / _POS_WITH_CENTER (~lines 306-355):
 *     triclinic vs rectangular PBC wrapping.
 *
 * ------------------------------------------------------------------------------
 * Backend contrast (why the same macro differs across files):
 *   CUDA here            OpenCL (common.cl)          role
 *   KERNEL=extern"C"     KERNEL=__kernel             entry point; extern "C"
 *     __global__                                     defeats C++ name mangling so
 *                                                    getKernel() can look it up.
 *   DEVICE=__device__    DEVICE=(empty)              non-entry function.
 *   GLOBAL=(empty)       GLOBAL=__global             CUDA global pointers need no
 *                                                    qualifier; OpenCL requires one.
 *   LOCAL_ARG=(empty)    LOCAL_ARG=__local           shared-mem *parameter*
 *                                                    qualifier (see LOCAL group).
 *   LOCAL_ID=threadIdx.x LOCAL_ID=get_local_id(0)    indexing model.
 *
 * The original one-line header comment ("CUDA definitions for the macros and
 * functions needed for the common compute framework") was accurate but minimal;
 * it is preserved in spirit and expanded here. Nothing in it was misleading.
 */

/**
 * @name Function qualifiers
 * @brief Map portable qualifier names onto CUDA function-space specifiers.
 *
 * KERNEL  -> `extern "C" __global__`: a grid entry point launchable from the
 *            host. `extern "C"` suppresses C++ name mangling so the mangled-name-
 *            agnostic host lookup (CudaContext::getKernel) finds the symbol by its
 *            plain source name. (OpenCL: `__kernel`; no extern "C" needed.)
 * DEVICE  -> `__device__`: a helper callable only from device code. (OpenCL:
 *            empty — non-kernel functions are implicitly device functions.)
 * @{
 */
#define KERNEL extern "C" __global__
#define DEVICE __device__
/** @} */

/**
 * @name Memory-space qualifiers
 * @brief Map portable address-space names onto CUDA storage/pointer qualifiers.
 *
 * LOCAL     -> `__shared__`: block-shared, on-chip memory shared by all threads
 *              in a thread block (the "local" work-group in OpenCL terminology).
 *              NOTE the naming inversion vs. OpenCL, where "local" == block-shared;
 *              CUDA's own "local memory" (per-thread spill) is unrelated and NOT
 *              what this maps to.
 * LOCAL_ARG -> (empty): qualifier applied to a *shared-memory pointer parameter*.
 *              CUDA declares `__shared__` arrays at function/file scope and passes
 *              plain pointers, so no parameter qualifier is needed. OpenCL passes
 *              `__local` pointers as kernel arguments, so there LOCAL_ARG=__local.
 *              This macro is what makes the same parameter list compile on both.
 * GLOBAL    -> (empty): pointer into global device memory. CUDA global pointers
 *              carry no qualifier; OpenCL requires `__global`, hence GLOBAL there.
 * RESTRICT  -> `__restrict__`: no-alias promise enabling load/store reordering and
 *              caching; identical intent to C99 `restrict` (OpenCL: `restrict`).
 * @{
 */
#define LOCAL __shared__
#define LOCAL_ARG
#define GLOBAL
#define RESTRICT __restrict__
/** @} */

/**
 * @name Thread/work-item indexing
 * @brief Portable NDRange indexing over CUDA's threadIdx/blockIdx/blockDim/gridDim.
 *
 * The common framework uses OpenCL's flat, 1-D NDRange model (get_local_id(0),
 * get_global_id(0), ...). These macros reconstruct that model from CUDA's
 * hierarchical built-ins; all kernels in this framework are launched 1-D.
 *
 * LOCAL_ID    -> `threadIdx.x`               : thread index within its block.
 * LOCAL_SIZE  -> `blockDim.x`                : threads per block (work-group size).
 * GLOBAL_ID   -> `blockIdx.x*blockDim.x+threadIdx.x` : unique global thread index.
 * GLOBAL_SIZE -> `blockDim.x*gridDim.x`      : total threads in the grid; used as
 *                the stride in grid-stride loops (`for i in [GLOBAL_ID; N;
 *                GLOBAL_SIZE)`), the framework's standard data-parallel idiom.
 * GROUP_ID    -> `blockIdx.x`                : block (work-group) index.
 * NUM_GROUPS  -> `gridDim.x`                 : number of blocks in the grid.
 *
 * @pre 1-D launch geometry (only the .x dimension is consulted).
 * @{
 */
#define LOCAL_ID threadIdx.x
#define LOCAL_SIZE blockDim.x
#define GLOBAL_ID (blockIdx.x*blockDim.x+threadIdx.x)
#define GLOBAL_SIZE (blockDim.x*gridDim.x)
#define GROUP_ID blockIdx.x
#define NUM_GROUPS gridDim.x
/** @} */

/**
 * @name Synchronization and memory ordering
 * @brief Portable barrier / fence primitives.
 *
 * SYNC_THREADS -> `__syncthreads();`: block-wide barrier. Every thread in the
 *                 block must reach it (do not place in divergent control flow) and
 *                 it also orders shared-memory accesses across the barrier.
 * MEM_FENCE    -> `__threadfence_block();`: block-scope memory fence — makes this
 *                 thread's prior writes visible to other threads in the block,
 *                 WITHOUT the barrier/convergence of SYNC_THREADS. NOTE: narrower
 *                 than the OpenCL mapping, where MEM_FENCE fences BOTH local and
 *                 global memory (`mem_fence(CLK_LOCAL_MEM_FENCE+CLK_GLOBAL_MEM_FENCE)`);
 *                 __threadfence_block orders shared+global writes but only at
 *                 block scope. Adequate for intra-block producer/consumer patterns;
 *                 cross-block ordering needs __threadfence()/__threadfence_system,
 *                 which this macro does NOT provide — a point to check when auditing
 *                 kernels that assume grid-wide visibility.
 * @{
 */
#define SYNC_THREADS __syncthreads();
#define MEM_FENCE __threadfence_block();
/** @} */

/**
 * @brief Portable atomic add; maps to CUDA hardware `atomicAdd`.
 * @param dest  pointer to the accumulator (global or shared memory).
 * @param value increment; overload selected by its type.
 *
 * On CUDA this is a single hardware atomic. Its dominant use in this framework is
 * accumulating the per-atom force buffer, which is typed `long long`
 * (CudaContext.cpp ~line 257, `force.initialize<long long>(..., paddedNumAtoms*3)`)
 * and holds fixed-point values produced by realToFixedPoint(). Integer atomicAdd
 * is associative and commutative and therefore order-independent and bit-exact —
 * see realToFixedPoint for why that matters. Contrast: OpenCL's ATOMIC_ADD maps to
 * `atom_add`, and on devices lacking cl_khr_int64_base_atomics common.cl provides a
 * software 64-bit emulation (split lo/hi 32-bit atomics with carry); CUDA needs no
 * such emulation (see SUPPORTS_64_BIT_ATOMICS).
 */
#define ATOMIC_ADD(dest, value) atomicAdd(dest, value)

/**
 * @brief Largest finite single-precision float, as a sentinel constant.
 *
 * Defined here because NVRTC compiles the kernel string without the host `<cfloat>`
 * header. Kernels use it as a "+infinity" initializer for min-reductions, distance
 * cutoffs, etc. NOTE: this is always the FLOAT max even in double-precision builds
 * (where `real`==double); it is a named large constant, not `real`-typed — fine for
 * sentinel use, but do not read it as "max representable real".
 */
#define FLT_MAX 3.40282347e+38f

/**
 * @name Portable 64-bit integer typedefs
 * @brief mm_long / mm_ulong = guaranteed-64-bit signed/unsigned integers.
 *
 * On CUDA `long long` / `unsigned long long` are 64-bit on every supported target,
 * whereas plain `long` is platform-dependent (32-bit under the Windows/LLP64 ABI),
 * so the framework standardizes on the `long long` family. These are the storage
 * type for fixed-point force/energy accumulators. Contrast: OpenCL (common.cl)
 * typedefs mm_long to `long`, which is already exactly 64-bit in OpenCL C.
 * @{
 */
using mm_long = long long;
using mm_ulong = unsigned long long;
/** @} */

/**
 * @name Backend capability flags
 * @brief Advertise features the CUDA backend always has.
 *
 * SUPPORTS_64_BIT_ATOMICS (=1): native 64-bit integer atomics — enables the
 *   fixed-point force accumulation path unconditionally. (OpenCL gates the
 *   equivalent on the cl_khr_int64_base_atomics extension and otherwise emulates.)
 * SUPPORTS_DOUBLE_PRECISION (=1): device supports `double`, so double- and
 *   mixed-precision kernels may be built. Both are compile-time constants because
 *   every CUDA-capable device OpenMM targets provides these; kernels `#if` on them
 *   to select code paths.
 * @{
 */
#define SUPPORTS_64_BIT_ATOMICS 1
#define SUPPORTS_DOUBLE_PRECISION 1
/** @} */

/**
 * @brief Convert a floating-point quantity to Q32 fixed-point for deterministic
 *        atomic accumulation.
 *
 * @param value  a `real`-typed scalar — a force component (kJ/mol/nm) in the
 *               dominant use, or any quantity summed atomically across threads.
 *               `real` is float or double per the build's precision mode (typedef
 *               injected by createModule above this file).
 * @return the value scaled by 2^32 and truncated to a signed 64-bit integer
 *         (`long long`), i.e. a fixed-point number with 32 fractional bits.
 *
 * @par Why fixed point.
 * Per-atom forces are the sum of contributions from many threads (bonded, angle,
 * torsion, nonbonded tiles, ...) accumulated via ATOMIC_ADD in arbitrary,
 * run-to-run-varying thread order. Floating-point addition is NOT associative:
 * `(a+b)+c != a+(b+c)` under rounding, so a floating atomicAdd yields results that
 * depend on scheduling and are not bit-reproducible. Integer addition IS
 * associative and commutative, so summing fixed-point integers with atomicAdd is
 * order-independent and bit-exact — identical totals regardless of how the GPU
 * interleaves threads. This determinism is the reason the force buffer is typed
 * `long long` (CudaContext.cpp ~line 257) rather than `real`.
 *
 * @par The 0x100000000 (= 2^32) scale factor.
 * Multiplying by 2^32 places 32 bits below the binary point (Q<32>.<32> layout in
 * the 64-bit word): the integer part occupies the high 32 bits, the fraction the
 * low 32. This fixes the accumulation grid (ULP == 2^-32 in force units), which is
 * what makes summation exact. A signed 64-bit accumulator then holds forces with
 * magnitude up to ~2^31 before overflow — ample for physical per-atom forces.
 * Downstream code divides the accumulated integer back by 2^32 (multiplies by
 * 2^-32) to recover a floating-point force. The identical constant and scheme are
 * used by the OpenCL backend's realToFixedPoint, keeping the two backends
 * bitwise-comparable.
 *
 * @note `static_cast<long long>` truncates toward zero; the fractional bits below
 *       2^-32 are discarded per contribution. `inline` + `__device__`: header-style
 *       device helper, no linkage cost.
 */
__device__ inline auto realToFixedPoint(real value) -> long long {
    return static_cast<long long>(value * 0x100000000);
}
