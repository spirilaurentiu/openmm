/**
 * @file utilities.cu
 * @brief General-purpose CUDA device kernels for the OpenMM CUDA platform:
 *        flat-buffer zeroing, block-partial energy reduction, and per-atom
 *        charge scatter into the packed position array.
 *
 * @par JIT compilation and injected types
 * These kernels are not built by nvcc. OpenMM string-encodes this file into
 * CudaKernelSources.cpp and compiles it at runtime with NVRTC. Ahead of this
 * text, CudaContext::createModule prepends a generated typedef preamble
 * (in openmm/platforms/cuda/src/CudaContext.cpp) that defines the element
 * types named below. Comments are stripped before encoding
 * (openmm/cmake_modules/strip_comments.py), so this documentation has no runtime
 * cost.
 *
 * @par Injected element types (compile-time parameters)
 * - real, real4: position/charge precision. double and double4 when the platform
 *   precision property is "double" (useDoublePrecision), otherwise float and
 *   float4 (both "single" and "mixed"). Selected in the CudaContext ctor
 *   and emitted by CudaContext::createModule.
 * - mixed, mixed2..4: energy/accumulator precision. double when precision is
 *   "double" OR "mixed" (useDoublePrecision || useMixedPrecision), otherwise
 *   float. Emitted by CudaContext::createModule. This keeps a double energy
 *   accumulator over single-precision positions in mixed mode.
 * - int4, make_int4: CUDA built-in vector type and constructor; not injected.
 *
 * @par Launch and stream conventions shared by every kernel here
 * All launch through CudaContext::executeKernel as a 1-D
 * grid of 1-D blocks on the context's current stream. executeKernel caps the
 * grid at numThreadBlocks (in CudaContext::executeKernel), so the grid can be smaller than
 * the data; each kernel therefore covers its range with a grid-stride loop and
 * is correct for any grid size. No kernel here synchronizes across blocks. All
 * launches are asynchronous with respect to the host unless a caller synchronizes.
 */

extern "C" {

/**
 * @brief Zero @p size contiguous 32-bit words of a device buffer.
 *
 * Every word in [0, @p size) is set to a zero bit pattern, which reads back as
 * 0, 0.0f, 0.0, or 0LL for any 32- or 64-bit element the buffer may hold. Serves
 * the clear{One..Six}Buffers kernels, which reset force, energy, and parameter
 * accumulators between force evaluations.
 *
 * @param[out] buffer Device buffer, borrowed; overwritten with zero, never read.
 *                    Accessed 16 bytes at a time, so the base address must be
 *                    16-byte aligned (satisfied by the CUDA allocator). Must not
 *                    alias any other buffer cleared in the same launch.
 * @param[in]  size   Length in 32-bit WORDS, not bytes and not elements. The host
 *                    forms it as byteSize/4 (CudaContext::clearBuffer;
 *                    CudaContext::addAutoclearBuffer).
 *
 * @pre All threads of the grid call this; participation is grid-wide and no block
 *      or grid synchronization occurs or is required.
 * @post Each of the @p size words is written exactly once; the result is
 *       independent of grid shape and deterministic.
 * @note Any grid/block shape is valid; the grid is capped at numThreadBlocks and
 *       the full range is covered by a grid-stride sweep. Callers use 128
 *       threads/block.
 */
__device__ void clearSingleBuffer(int* __restrict__ buffer, int size) {
    int index = blockDim.x*blockIdx.x+threadIdx.x;
    int4* buffer4 = (int4*) buffer;
    int sizeDiv4 = size/4;
    while (index < sizeDiv4) {
        buffer4[index] = make_int4(0);
        index += blockDim.x*gridDim.x;
    }
    if (blockDim.x*blockIdx.x+threadIdx.x == 0)
        for (int i = sizeDiv4*4; i < size; i++)
            buffer[i] = 0;
}

/**
 * @brief Kernel entry point that zeroes one device buffer.
 *
 * @param[out] buffer Device buffer to zero; see clearSingleBuffer for the
 *                    alignment and write contract.
 * @param[in]  size   Length in 32-bit words (byteSize/4).
 *
 * @par Launch configuration
 * 1-D grid, 128 threads/block; grid = min(ceil(size/128), numThreadBlocks)
 * (CudaContext::clearBuffer then CudaContext::executeKernel).
 * @par Shared memory
 * None.
 * @note Launched directly by CudaContext::clearBuffer and, for a lone leftover
 *       autoclear buffer, by CudaContext::clearAutoclearBuffers.
 *       Asynchronous on the current stream.
 */
__global__ void clearBuffer(int* __restrict__ buffer, int size) {
    // size counts 32-bit WORDS (host CudaContext::clearBuffer passes byteSize/4), not bytes or elements - an easy caller mistake.
    clearSingleBuffer(buffer, size);
}

/**
 * @brief Kernel entry point that zeroes two independent device buffers in one
 *        launch.
 *
 * The clear{Two..Six}Buffers family coalesces several independent buffer clears,
 * each with its own length, into a single launch to amortize launch overhead. A
 * thread does the work for every buffer its index reaches; the buffers are
 * disjoint, so the order among them is irrelevant.
 *
 * @param[out] buffer1 First device buffer to zero; see clearSingleBuffer.
 * @param[in]  size1   Length of @p buffer1 in 32-bit words.
 * @param[out] buffer2 Second device buffer to zero; disjoint from @p buffer1.
 * @param[in]  size2   Length of @p buffer2 in 32-bit words.
 *
 * @par Launch configuration
 * 1-D grid, 128 threads/block; grid sized to max(size1,size2) words
 * (CudaContext::clearAutoclearBuffers).
 * @par Shared memory
 * None.
 * @note Emitted by the autoclear pass when exactly two buffers remain after
 *       grouping in sixes. Asynchronous on the current stream.
 */
__global__ void clearTwoBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
}

/**
 * @brief Kernel entry point that zeroes three independent device buffers in one
 *        launch. Shared semantics are documented on clearTwoBuffers.
 *
 * @param[out] buffer1 First device buffer to zero; see clearSingleBuffer.
 * @param[in]  size1   Length of @p buffer1 in 32-bit words.
 * @param[out] buffer2 Second device buffer to zero; disjoint from the others.
 * @param[in]  size2   Length of @p buffer2 in 32-bit words.
 * @param[out] buffer3 Third device buffer to zero; disjoint from the others.
 * @param[in]  size3   Length of @p buffer3 in 32-bit words.
 *
 * @par Launch configuration
 * 1-D grid, 128 threads/block; grid sized to max(size1,size2,size3) words
 * (CudaContext::clearAutoclearBuffers).
 * @par Shared memory
 * None.
 */
__global__ void clearThreeBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
}

/**
 * @brief Kernel entry point that zeroes four independent device buffers in one
 *        launch. Shared semantics are documented on clearTwoBuffers.
 *
 * @param[out] buffer1 First device buffer to zero; see clearSingleBuffer.
 * @param[in]  size1   Length of @p buffer1 in 32-bit words.
 * @param[out] buffer2 Second device buffer to zero; disjoint from the others.
 * @param[in]  size2   Length of @p buffer2 in 32-bit words.
 * @param[out] buffer3 Third device buffer to zero; disjoint from the others.
 * @param[in]  size3   Length of @p buffer3 in 32-bit words.
 * @param[out] buffer4 Fourth device buffer to zero; disjoint from the others.
 * @param[in]  size4   Length of @p buffer4 in 32-bit words.
 *
 * @par Launch configuration
 * 1-D grid, 128 threads/block; grid sized to max(size1..size4) words
 * (CudaContext::clearAutoclearBuffers).
 * @par Shared memory
 * None.
 */
__global__ void clearFourBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3, int* __restrict__ buffer4, int size4) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
}

/**
 * @brief Kernel entry point that zeroes five independent device buffers in one
 *        launch. Shared semantics are documented on clearTwoBuffers.
 *
 * @param[out] buffer1 First device buffer to zero; see clearSingleBuffer.
 * @param[in]  size1   Length of @p buffer1 in 32-bit words.
 * @param[out] buffer2 Second device buffer to zero; disjoint from the others.
 * @param[in]  size2   Length of @p buffer2 in 32-bit words.
 * @param[out] buffer3 Third device buffer to zero; disjoint from the others.
 * @param[in]  size3   Length of @p buffer3 in 32-bit words.
 * @param[out] buffer4 Fourth device buffer to zero; disjoint from the others.
 * @param[in]  size4   Length of @p buffer4 in 32-bit words.
 * @param[out] buffer5 Fifth device buffer to zero; disjoint from the others.
 * @param[in]  size5   Length of @p buffer5 in 32-bit words.
 *
 * @par Launch configuration
 * 1-D grid, 128 threads/block; grid sized to max(size1..size5) words
 * (CudaContext::clearAutoclearBuffers).
 * @par Shared memory
 * None.
 */
__global__ void clearFiveBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3, int* __restrict__ buffer4, int size4, int* __restrict__ buffer5, int size5) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
    clearSingleBuffer(buffer5, size5);
}

/**
 * @brief Kernel entry point that zeroes six independent device buffers in one
 *        launch. Shared semantics are documented on clearTwoBuffers.
 *
 * Six is the widest variant and the primary path of the autoclear pass, which
 * consumes registered buffers in groups of six and falls through to the narrower
 * variants for the trailing 1-5.
 *
 * @param[out] buffer1 First device buffer to zero; see clearSingleBuffer.
 * @param[in]  size1   Length of @p buffer1 in 32-bit words.
 * @param[out] buffer2 Second device buffer to zero; disjoint from the others.
 * @param[in]  size2   Length of @p buffer2 in 32-bit words.
 * @param[out] buffer3 Third device buffer to zero; disjoint from the others.
 * @param[in]  size3   Length of @p buffer3 in 32-bit words.
 * @param[out] buffer4 Fourth device buffer to zero; disjoint from the others.
 * @param[in]  size4   Length of @p buffer4 in 32-bit words.
 * @param[out] buffer5 Fifth device buffer to zero; disjoint from the others.
 * @param[in]  size5   Length of @p buffer5 in 32-bit words.
 * @param[out] buffer6 Sixth device buffer to zero; disjoint from the others.
 * @param[in]  size6   Length of @p buffer6 in 32-bit words.
 *
 * @par Launch configuration
 * 1-D grid, 128 threads/block; grid sized to max(size1..size6) words
 * (CudaContext::clearAutoclearBuffers).
 * @par Shared memory
 * None.
 */
__global__ void clearSixBuffers(int* __restrict__ buffer1, int size1, int* __restrict__ buffer2, int size2, int* __restrict__ buffer3, int size3, int* __restrict__ buffer4, int size4, int* __restrict__ buffer5, int size5, int* __restrict__ buffer6, int size6) {
    clearSingleBuffer(buffer1, size1);
    clearSingleBuffer(buffer2, size2);
    clearSingleBuffer(buffer3, size3);
    clearSingleBuffer(buffer4, size4);
    clearSingleBuffer(buffer5, size5);
    clearSingleBuffer(buffer6, size6);
}

/**
 * @brief Reduce the per-lane energy buffer to one partial sum per block.
 *
 * Force kernels accumulate potential energy into @p energyBuffer, one slot per
 * energy lane. Each block sums a strided slice of that buffer and writes its
 * single partial to result[blockIdx.x]. This is a PARTIAL reduction: the final
 * sum over the per-block partials is completed on the host
 * (CudaContext::reduceEnergy).
 *
 * @param[in]  energyBuffer  Device array of @p bufferSize energy lanes, mixed
 *                           precision, read-only, borrowed. Host passes
 *                           energyBuffer.getSize() as @p bufferSize.
 * @param[out] result        Device output, mixed precision, one entry per block;
 *                           result[blockIdx.x] receives this block's partial.
 *                           Must hold at least gridDim.x entries. Borrowed.
 * @param[in]  bufferSize    Number of valid lanes in @p energyBuffer.
 * @param[in]  workGroupSize Reduction width; equals blockDim.x. Host passes 512
 *                           (CudaContext::reduceEnergy).
 *
 * @par Launch configuration
 * 1-D grid, blockDim.x = @p workGroupSize = 512; grid = energySum.getSize()
 * blocks, which equals the device multiprocessor count and is not clamped by
 * numThreadBlocks (the reduceEnergy launch in CudaContext::reduceEnergy, energySum
 * sized to the multiprocessor count in CudaContext::initialize, and numThreadBlocks
 * set in the CudaContext constructor). Each block fills exactly one
 * result slot, so every result entry is written.
 * @par Shared memory
 * Dynamic: @p workGroupSize * sizeof(mixed) bytes, passed as
 * workGroupSize * energyBuffer.getElementSize() (CudaContext::reduceEnergy).
 *
 * @pre blockDim.x == @p workGroupSize: the dynamic shared buffer is sized to
 *      workGroupSize and the in-block combine spans [0, workGroupSize); a
 *      mismatch drops contributions or accesses shared memory out of bounds.
 * @post result[blockIdx.x] holds the sum of this block's assigned lanes. The
 *       in-block combine reaches a barrier that all threads of the block must
 *       execute, so callers must not let threads of the block exit early. No
 *       cross-block synchronization occurs; the cross-block sum is the host's job.
 * @note Summation order within a block is fixed for a given launch shape, so each
 *       partial is deterministic; mixed = double (except on a single-precision
 *       platform) preserves energy accuracy over float positions and forces.
 */
__global__ void reduceEnergy(const mixed* __restrict__ energyBuffer, mixed* __restrict__ result, int bufferSize, int workGroupSize) {
    extern __shared__ mixed tempBuffer[];
    const unsigned int thread = threadIdx.x;
    mixed sum = 0;
    // Grid is only energySum-size blocks (== multiprocessor count, per host CudaContext::reduceEnergy), far fewer than bufferSize lanes, so each thread folds many lanes here.
    for (unsigned int index = blockDim.x*blockIdx.x+threadIdx.x; index < bufferSize; index += blockDim.x*gridDim.x)
        sum += energyBuffer[index];
    // Shared tempBuffer is sized to workGroupSize by the host; correctness relies on blockDim.x == workGroupSize (host passes both = 512).
    tempBuffer[thread] = sum;
    for (int i = 1; i < workGroupSize; i *= 2) {
        __syncthreads();
        if (thread%(i*2) == 0 && thread+i < workGroupSize)
            tempBuffer[thread] += tempBuffer[thread+i];
    }
    if (thread == 0)
        // Per-block PARTIAL only; host CudaContext::reduceEnergy sums all blocks. grid == result size (unclamped), so every result slot is written.
        result[blockIdx.x] = tempBuffer[0];
}

/**
 * @brief Scatter per-atom charges into the w lane of the packed position array.
 *
 * OpenMM packs each atom as real4 (x, y, z, q) in @p posq. This kernel rewrites
 * only the charge lane (.w), leaving positions untouched; it applies charge
 * changes that occur without atoms moving. @p charges is indexed in the original
 * atom order while @p posq is in the engine's internal (reordered) order, so
 * @p atomOrder maps each internal slot back to its source charge.
 *
 * @param[in]  charges   Device array of charges in original atom order, real
 *                       precision, read-only, borrowed. Host uploads
 *                       chargeBuffer of numAtoms entries (CudaContext::setCharges).
 * @param[in,out] posq   Device position+charge array in internal order, real4.
 *                       Only .w is modified; x, y, z are preserved. Borrowed.
 * @param[in]  atomOrder Device permutation; atomOrder[i] indexes @p charges for
 *                       internal slot i (host passes atomIndexDevice in
 *                       CudaContext::setCharges).
 * @param[in]  numAtoms  Count of real atoms; internal slots [0, numAtoms) are
 *                       updated, padded slots beyond it are left unchanged.
 *
 * @pre For every i in [0, @p numAtoms), atomOrder[i] lies in [0, @p numAtoms):
 *      @p charges holds only numAtoms entries.
 * @par Launch configuration
 * 1-D grid, 64 threads/block (default CudaContext::ThreadBlockSize);
 * grid = min(ceil(numAtoms/64), numThreadBlocks) (CudaContext::setCharges then CudaContext::executeKernel).
 * @par Shared memory
 * None.
 * @post Each updated slot's charge is written exactly once; work is independent
 *       per atom and no synchronization occurs. Asynchronous on the current stream.
 */
__global__ void setCharges(real* __restrict__ charges, real4* __restrict__ posq, int* __restrict__ atomOrder, int numAtoms) {
    for (int i = blockDim.x*blockIdx.x+threadIdx.x; i < numAtoms; i += blockDim.x*gridDim.x)
        // atomOrder is CudaContext::atomIndexDevice (internal->original permutation); charges (chargeBuffer) holds only numAtoms entries, so atomOrder[i] must be < numAtoms.
        posq[i].w = charges[atomOrder[i]];
}
}
