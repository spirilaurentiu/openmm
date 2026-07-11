/**
 * @file parallel.cu
 * @brief Device-side per-atom force reduction for OpenMM's multi-GPU CUDA platform.
 *
 * @par Role
 * This source is JIT-compiled at runtime (NVRTC) from the string table
 * CudaKernelSources::parallel and launched, not with the `<<< >>>` syntax, but
 * via cuLaunchKernel through CudaContext::executeKernel. The host driver is the
 * class CudaParallelCalcForcesAndEnergyKernel; the sole kernel here, sumForces,
 * is launched from its finishComputation method. All of the kernel's contract
 * (buffer origins, launch shape, gather completion) is established there rather
 * than by injected preprocessor macros.
 *
 * @par Multi-GPU force-combining strategy (host context)
 * A parallel run spawns N CUDA contexts, one per participating device. Every
 * context evaluates forces on the SAME full atom set over a load-balanced share
 * of the nonbonded work plus its other contributions, so each context's own
 * force buffer (its getForce array) is a PARTIAL per-atom force. The physical
 * total on each atom is the element-wise sum of the N partial buffers. The
 * primary context (index 0) keeps its partial forces in place; the N-1 secondary
 * contexts' partial buffers are gathered onto device 0 into one contiguous array
 * (the host CudaArray contextForces). The gather is performed by the secondary
 * contexts inside FinishComputationTask::execute:
 *   - peer path (peerAccessSupported): each secondary context issues cuMemcpyAsync
 *     from its getForce array into contextForces at offset (contextIndex-1)*bufferSize,
 *     and the primary stream is made to wait on each secondary's completion event
 *     before the launch;
 *   - staging path (no peer access): each secondary downloads its getForce array
 *     into a slice of the portable pinned host buffer pinnedForceBuffer, and
 *     finishComputation uploads the whole pinned buffer into contextForces before
 *     the launch.
 * sumForces then folds all N-1 gathered buffers into context 0's force buffer in
 * place, leaving the complete total force on the primary device for the integrator.
 *
 * @par Force representation and determinism
 * Forces are 64-bit fixed-point integers (long long): a component in kJ/mol/nm is
 * stored as round(value * 2^32), OpenMM's deterministic accumulation format. The
 * scale is identical on every device, so partial forces combine by plain integer
 * addition with no rescaling. Integer addition is associative and commutative, so
 * the total is bitwise identical regardless of device count, gather order, or the
 * launch's block/grid shape. Accumulator overflow is not checked (wraps mod 2^64).
 *
 * @par Buffer layout
 * A full force buffer is bufferSize = 3 * paddedNumAtoms int64 elements in
 * component-major order (all x, then all y, then all z). sumForces is agnostic to
 * this split: it treats each buffer as a flat length-bufferSize vector and reduces
 * element-wise, so the x/y/z structure matters only to the caller.
 *
 * @par Compile-time parameters
 * None. The module is built from CudaKernelSources::parallel with an empty defines
 * map, so unlike typical OpenMM CUDA kernels no NUM_ATOMS / PADDED_NUM_ATOMS /
 * device-count macros are injected. The two quantities such macros would supply,
 * the buffer length and the number of contributing devices, arrive instead as the
 * runtime scalar arguments bufferSize and numBuffers.
 */

/**
 * @brief Add every secondary context's partial force buffer into the primary device's force buffer, in place.
 *
 * @par Effect
 * For every element index in [0, @p bufferSize), computes
 *   force[index] = force[index] + sum over d in [0, numBuffers) of buffer[index + d*bufferSize]
 * i.e. context 0's own partial force (already in @p force) plus the partial forces
 * of the @p numBuffers secondary contexts concatenated in @p buffer. Each output
 * element is written exactly once, from a fixed-order integer sum; the write is a
 * deterministic function of the inputs independent of launch configuration. See
 * the file header for the surrounding multi-GPU strategy.
 *
 * @par Work decomposition
 * One logical thread per output element index; a grid-stride loop lets a
 * grid smaller than @p bufferSize cover the whole range. Threads own disjoint
 * output indices with no cross-thread dependence.
 *
 * @par Launch configuration
 * Launched via CudaContext::executeKernel(sumKernel, args, bufferSize) on the
 * primary context's current stream: 1D grid of blocks of the context's default
 * ThreadBlockSize, grid capped at the context's numThreadBlocks (hence the
 * grid-stride loop is required for correctness, not merely an optimization). No
 * dynamic shared memory. Result is valid on device 0 once that stream completes;
 * the call is asynchronous with respect to the host.
 *
 * @par Synchronization scope
 * None across blocks or threads: no __syncthreads, no shared memory, no atomics.
 * Correctness relies solely on the disjoint-index write pattern.
 *
 * @param[in,out] force     Primary device (context 0) force buffer, @p bufferSize
 *                          int64 fixed-point elements. On entry holds context 0's
 *                          partial forces; on exit holds the total force across all
 *                          contexts. Resides on and is borrowed from device 0
 *                          (context 0's getForce device pointer).
 * @param[in]     buffer    Gathered partial forces of the @p numBuffers secondary
 *                          contexts, concatenated: block d (d in [0, numBuffers))
 *                          occupies buffer[d*bufferSize .. (d+1)*bufferSize) and
 *                          holds the full force buffer of secondary context d+1.
 *                          Resides on device 0 (the contextForces device pointer),
 *                          length numBuffers*bufferSize int64. Must not alias
 *                          @p force (distinct allocations; asserted by __restrict__).
 * @param[in]     bufferSize Number of int64 elements in one force buffer
 *                          (3 * paddedNumAtoms): simultaneously the per-context
 *                          stride within @p buffer and the number of output
 *                          elements produced.
 * @param[in]     numBuffers Number of secondary contexts contributing to @p buffer
 *                          (total device count minus one). When 0 (single device)
 *                          the inner sum is empty and @p force is left unchanged.
 *
 * @pre @p force and @p buffer reside on the launching (primary) device; @p buffer
 *      holds at least numBuffers*bufferSize and @p force at least bufferSize valid
 *      int64 elements.
 * @pre All secondary partial buffers are already gathered into @p buffer and
 *      visible on device 0 before launch. finishComputation guarantees this by
 *      making the primary stream wait on each secondary's peer-copy event (peer
 *      path) or by uploading the pinned staging buffer (staging path) first.
 * @pre @p bufferSize and @p numBuffers are non-negative and consistent with the
 *      allocations above.
 *
 * @post Each element of @p force equals the sum of the corresponding element of
 *       the primary partial buffer and all secondary partial buffers.
 *
 * @note Purely memory-bound streaming: one read of @p force, @p numBuffers reads
 *       of @p buffer, and one write of @p force per output element.
 */
extern "C" __global__ void sumForces(long long* __restrict__ force, long long* __restrict__ buffer, int bufferSize, int numBuffers) {
    int totalSize = bufferSize*numBuffers;
    for (int index = blockDim.x*blockIdx.x+threadIdx.x; index < bufferSize; index += blockDim.x*gridDim.x) {
        long long sum = force[index];
        // Stride-bufferSize walk visits the same component index in each secondary
        // context's block: block d holds context d+1's partial force, placed at
        // offset d*bufferSize by finishComputation (peer cuMemcpyAsync or pinned-host
        // staging). This layout is fixed on the host; the strided read cannot be made
        // contiguous without changing those gather offsets too.
        for (int i = index; i < totalSize; i += bufferSize)
            sum += buffer[i];
        force[index] = sum;
    }
}
