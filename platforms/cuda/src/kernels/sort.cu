/**
 * @file sort.cu
 * @brief OpenMM's GPU distribution ("bucket") sort, JIT-specialized per element type.
 *
 * Device side of OpenMM's general-purpose GPU sort. The host driver is class CudaSort
 * (CudaSort.cpp): its constructor JIT-compiles this source with NVRTC after
 * string-substituting the compile-time parameters below, and CudaSort::sort() launches the
 * kernels in sequence through CudaContext::executeKernel, which issues a 1D cuLaunchKernel
 * with gridSize = min(ceil(threads/blockSize), context.getNumThreadBlocks()). These kernels
 * are NOT invoked with the <<< >>> syntax and are not visible to callers as C symbols; the
 * only entry point is CudaSort::sort(). Robosample's dominant user is
 * CudaNonbondedUtilities (neighbor-list construction: sorts the atom-block list by an
 * interaction-locality key, uniform=false); OpenMM core also uses it for
 * CommonCalcNonbondedForce (PME atom-grid index) and CommonCalcConstantPotentialForce.
 *
 * @par Sort contract (guaranteed on return of CudaSort::sort)
 *      The array is replaced by a permutation of itself in nondecreasing key order
 *      (key = SORT_KEY applied to each element). Every input element appears exactly once;
 *      no element is created or dropped. Among elements with equal keys the relative order
 *      is unspecified on the bucket path (see @ref determinism) and stable (by original
 *      index) on the sortShortList2 path. The sort is a pure device-side operation on the
 *      caller's array; results are valid after the stream reaches the last kernel.
 *
 * @par Phases (bucket path — taken when the input does not fit one block; !isShortList)
 *      A single-pass distribution sort: keys are mapped to numBuckets contiguous, key-range-
 *      monotone buckets by value (not by radix digit), then each bucket is bitonic-sorted
 *      independently. The driver sizes buckets so the average holds ~sortKernelSize/2
 *      elements, small enough to sort in one block's shared memory. Concatenating sorted
 *      buckets in bucket order yields the globally sorted array because the bucketing is
 *      key-monotone. The five kernels form a pipeline; each depends on the previous:
 *        1. computeRange           — histogram-support: global min/max key (UNIFORM only)
 *                                     and zero the bucket counters.
 *        2. assignElementsToBuckets / assignElementsToBuckets2
 *                                  — bucketing: assign each element a bucket and reserve its
 *                                     intra-bucket slot via atomic counter increment.
 *        3. computeBucketPositions — prefix scan: inclusive scan of bucket sizes into
 *                                     cumulative end offsets.
 *        4. copyDataToBuckets      — scatter: place each element at bucketStart+slot.
 *        5. sortBuckets            — per-bucket bitonic sort into the final array.
 *
 * @par Short-list path (input fits in one launch — CudaSort::isShortList true)
 *      One of two kernels runs; the bucket phases do not. sortShortList2 (rank-by-counting,
 *      out-of-place, stable) is used when dataLength <= ThreadBlockSize*numThreadBlocks;
 *      otherwise sortShortList (single-block in-shared-memory bitonic). The driver picks
 *      exactly one path per sort() call; the two families never cooperate.
 *
 * @par Dispatch boundary and runtime-selected element type (see CudaSort::createSort)
 *      This file is not templated in C++; the element type is chosen at runtime by the
 *      SortTrait subclass passed to CudaContext::createSort / the CudaSort constructor,
 *      whose accessor strings are pasted in by CudaSort::CudaSort via
 *      context.replaceStrings() before NVRTC parses the file (comments are stripped first by
 *      strip_comments.py, so this documentation has zero JIT cost). The invariant that the
 *      trait's DATA_TYPE/KEY_TYPE/SORT_KEY actually match the real element type of the array
 *      later handed to sort() is not checkable by the compiler and is a @pre of
 *      CudaSort::createSort; sort() verifies only element SIZE (getElementSize vs
 *      getDataSize) and total length, not type layout, so a matching-size but wrong-layout
 *      trait is undefined behavior. The concrete traits that exist in this build are:
 *        - CudaNonbondedUtilities::BlockSortTrait  : DATA_TYPE=unsigned int, KEY_TYPE=
 *          unsigned int, SORT_KEY="value", getDataSize()==4, uniform=false. Robosample
 *          neighbor list.
 *        - CommonCalcNonbondedForce::SortTrait and CommonCalcConstantPotentialForce::
 *          SortTrait : DATA_TYPE=int2, KEY_TYPE=int, SORT_KEY="value.y", getDataSize()==8,
 *          uniform defaulted true. Sort (index,key) pairs by the key field.
 *        - the FloatTrait example in CudaSort.h / TestCudaSort.cpp : DATA_TYPE=float,
 *          KEY_TYPE=float, SORT_KEY="value".
 *      Any element type whose values are relocatable by plain copy and whose KEY_TYPE is an
 *      ordered scalar (min/max/</> and a float-cast are applied to it) is admissible; there
 *      is no runtime tag switch and no unsupported-type error path — the type is fixed at
 *      JIT time, one module per trait.
 *
 * @par Compile-time parameter meanings (string macros, not C++ template params)
 *      DATA_TYPE  the element stored, moved, and sorted (getDataType).
 *      KEY_TYPE   the ordered scalar compared; the result type of SORT_KEY (getKeyType).
 *                 Equals DATA_TYPE for self-keyed elements, a field type for structs.
 *      SORT_KEY   key-EXTRACTION expression (getSortKey), spliced as the body of getValue()
 *                 with a local `value` of type DATA_TYPE in scope (e.g. "value", "value.y").
 *      MIN_KEY    KEY_TYPE constant <= every key; the identity seed for the max reduction
 *                 in computeRange (getMinKey).
 *      MAX_KEY    KEY_TYPE constant >= every key; the identity seed for the min reduction in
 *                 computeRange and the key of MAX_VALUE (getMaxKey).
 *      MAX_VALUE  DATA_TYPE sentinel whose key equals MAX_KEY (getMaxValue); pads a bucket's
 *                 shared buffer so the fixed-width bitonic network sorts padding to the high
 *                 end and never writes it back.
 *      UNIFORM    "1" if the caller declared the key distribution ~uniform, else "0". It
 *                 #if-guards the range reduction (only the uniform map needs global min/max)
 *                 and, via kernel NAME selection in the driver, chooses the uniform vs
 *                 non-uniform bucketing kernel — only the matching one is launched.
 *
 * @par determinism
 *      No floating-point atomics touch element data; the only atomics are integer
 *      increments of bucket counters (assignElementsToBuckets*), whose race decides which
 *      slot equal-keyed elements land in. The final key order is therefore deterministic,
 *      but the relative order of equal-keyed elements on the bucket path varies run to run.
 *      For BlockSortTrait the key is the element itself, so no two distinct elements share a
 *      key and the whole result is deterministic.
 *
 * @note The bucket-assignment arithmetic casts keys to float regardless of KEY_TYPE (see
 *       assignElementsToBuckets*); this affects only which bucket a key lands in, never
 *       final correctness, since sortBuckets re-sorts on the true KEY_TYPE and bucket
 *       boundaries stay monotone in float.
 */

/**
 * @brief Extract the sort key from an element.
 *
 * The single point where element-to-key projection is defined; every comparison and every
 * bucket assignment routes through it, so key semantics are fixed by the injected SORT_KEY
 * expression alone. Pure and side-effect-free; safe to call from any thread with no
 * participation or synchronization requirement.
 *
 * @param[in] value  The element whose key is wanted; passed by value, not mutated.
 * @return Its KEY_TYPE sort key.
 */
__device__ KEY_TYPE getValue(DATA_TYPE value) {
    // SORT_KEY is the trait's getSortKey() text spliced in by CudaSort::replaceStrings; the
    // local 'value' (type DATA_TYPE) is the only name it may name (e.g. "value" or "value.y").
    return SORT_KEY;
}

extern "C" {

/**
 * @brief Short-list path: bitonic-sort an entire array in one block via shared memory.
 *
 * Phase role: standalone whole-array sort, no bucketing. On return @p data holds all
 * @p length elements in nondecreasing key order. The bitonic network is length-generalized:
 * @p length need not be a power of two. Not stable.
 *
 * @par Launch configuration
 *      Exactly one block (CudaSort::sort passes threads == blockSize == sortKernelSize, so
 *      the grid collapses to a single block). blockDim.x == sortKernelSize; the array is
 *      streamed through the block in strides of blockDim.x, so @p length may exceed
 *      blockDim.x. Chosen only on the short-list path when the array is too large for
 *      sortShortList2 yet still fits shared memory.
 * @par Shared memory
 *      Dynamic; the driver passes dataLength*getDataSize() bytes, i.e. one DATA_TYPE per
 *      element. @pre available dynamic shared memory >= length*sizeof(DATA_TYPE).
 * @par Synchronization
 *      Block-wide __syncthreads() after the load and after every compare-exchange stage; all
 *      threads of the block must reach every stage, so the launch must be non-divergent at
 *      block granularity. No cross-block synchronization (there is only one block).
 *
 * @param[in,out] data   Device array of DATA_TYPE, length @p length; sorted in place. Borrowed
 *                       (the caller retains ownership); every slot is read then overwritten.
 * @param[in]     length Number of valid elements; need not be a power of two.
 *
 * @note No MAX_VALUE padding: out-of-range bitonic partners are skipped by the ixj < length
 *       guard rather than padded.
 */
__global__ void sortShortList(DATA_TYPE* __restrict__ data, unsigned int length) {
    // Load the data into local memory.
    
    extern __shared__ DATA_TYPE dataBuffer[];
    for (int index = threadIdx.x; index < length; index += blockDim.x)
        dataBuffer[index] = data[index];
    __syncthreads();

    // Perform a bitonic sort in local memory.

    for (unsigned int k = 2; k < 2*length; k *= 2) {
        for (unsigned int j = k/2; j > 0; j /= 2) {
            for (unsigned int i = threadIdx.x; i < length; i += blockDim.x) {
                int ixj = i^j;
                if (ixj > i && ixj < length) {
                    DATA_TYPE value1 = dataBuffer[i];
                    DATA_TYPE value2 = dataBuffer[ixj];
                    bool ascending = ((i&k) == 0);
                    for (unsigned int mask = k*2; mask < 2*length; mask *= 2)
                        ascending = ((i&mask) == 0 ? !ascending : ascending);
                    KEY_TYPE lowKey  = (ascending ? getValue(value1) : getValue(value2));
                    KEY_TYPE highKey = (ascending ? getValue(value2) : getValue(value1));
                    if (lowKey > highKey) {
                        dataBuffer[i] = value2;
                        dataBuffer[ixj] = value1;
                    }
                }
            }
            __syncthreads();
        }
    }

    // Write the data back to global memory.

    for (int index = threadIdx.x; index < length; index += blockDim.x)
        data[index] = dataBuffer[index];
}

/**
 * @brief Short-list path: stable rank-by-counting sort, one thread per element.
 *
 * Phase role: standalone whole-array sort, out-of-place. Each element's final position is
 * its rank = the number of elements that sort strictly before it, ties broken by original
 * index; the result is therefore a stable, collision-free permutation. Cost is O(n) per
 * thread / O(n^2) total but fully parallel; the driver prefers it when
 * dataLength <= ThreadBlockSize*numThreadBlocks (one element per resident thread).
 *
 * @par Launch configuration
 *      Grid sized to cover every element (CudaSort::sort passes threads == dataLength,
 *      blockSize defaulting to ThreadBlockSize == 64). @pre the grid covers all elements;
 *      threads with globalId >= length load index 0 harmlessly and skip their store.
 * @par Shared memory
 *      Static, a fixed 64-element DATA_TYPE tile; no dynamic shared memory is passed.
 * @par Synchronization
 *      Block-wide __syncthreads() around each tile load; all threads of the block must
 *      participate in every tile. No cross-block synchronization.
 *
 * @param[in]  dataIn  Source array of DATA_TYPE, read-only, length @p length. Borrowed.
 * @param[out] dataOut Destination array, length @p length, distinct from @p dataIn; each
 *                     element is scattered to its rank exactly once. The driver passes the
 *                     `buckets` scratch here and copies it back over the caller's array.
 * @param[in]  length  Number of valid elements.
 *
 * @warning The shared tile is fixed at 64 elements while tiles are read in widths of
 *          blockDim.x; correctness requires blockDim.x <= 64. Satisfied today only because
 *          ThreadBlockSize == 64. A larger launch block would overflow the tile.
 */
__global__ void sortShortList2(const DATA_TYPE* __restrict__ dataIn, DATA_TYPE* __restrict__ dataOut, unsigned int length) {
    // 64 == CudaContext::ThreadBlockSize, the only block width CudaSort::sort launches this
    // kernel with; the tile loads below index dataBuffer[threadIdx.x], so blockDim.x must stay <= 64.
    __shared__ DATA_TYPE dataBuffer[64];
    int globalId = blockDim.x*blockIdx.x+threadIdx.x;
    DATA_TYPE value = dataIn[globalId < length ? globalId : 0];
    KEY_TYPE key = getValue(value);
    int count = 0;
    for (int blockStart = 0; blockStart < length; blockStart += blockDim.x) {
        int numInBlock = min(blockDim.x, length-blockStart);
        __syncthreads();
        if (threadIdx.x < numInBlock)
            dataBuffer[threadIdx.x] = dataIn[blockStart+threadIdx.x];
        __syncthreads();
        for (int i = 0; i < numInBlock; i++) {
            KEY_TYPE otherKey = getValue(dataBuffer[i]);
            if (otherKey < key || (otherKey == key && blockStart+i < globalId))
                count++;
        }
    }
    if (globalId < length)
        dataOut[count] = value;
}

/**
 * @brief Bucket phase 1: reduce global min/max key (UNIFORM only) and zero bucket counters.
 *
 * Phase role: histogram support. On return, when UNIFORM==1, @p range holds the exact
 * {min,max} key over the whole array (the domain the uniform bucketing map spans); when
 * UNIFORM==0 the reduction is compiled out and @p range is left untouched (the non-uniform
 * bucketing kernel estimates its own map from a sample). In both configurations every bucket
 * counter in @p bucketOffset is zeroed, establishing the precondition for the atomic
 * accumulation in phase 2.
 *
 * @par Launch configuration
 *      Exactly one block (CudaSort::sort passes threads == blockSize == rangeKernelSize).
 *      blockDim.x == rangeKernelSize, the largest power of two <= the device max block size
 *      (then clamped to length); @pre blockDim.x is a power of two, required by the tree
 *      reduction.
 * @par Shared memory
 *      Dynamic, used only when UNIFORM==1: the driver passes 2*rangeKernelSize*getKeySize()
 *      bytes, split into two KEY_TYPE halves (running min, running max). @pre when UNIFORM==1,
 *      available dynamic shared memory >= 2*blockDim.x*sizeof(KEY_TYPE). Zero is passed and
 *      required when UNIFORM==0.
 * @par Synchronization
 *      Block-wide __syncthreads() between reduction steps; all threads participate. Single
 *      block, so no cross-block sync.
 *
 * @param[in]  data         Source array, read-only, length @p length. Borrowed.
 * @param[in]  length       Number of elements.
 * @param[out] range        Two-element device buffer {min key, max key}; written only when
 *                          UNIFORM==1, otherwise left unchanged.
 * @param[in]  numBuckets   Number of bucket counters to clear.
 * @param[out] bucketOffset numBuckets-element counter array; every entry set to 0.
 *
 * @note MAX_KEY seeds the running minimum and MIN_KEY the running maximum (reduction
 *       identities).
 */
__global__ void computeRange(const DATA_TYPE* __restrict__ data, unsigned int length, KEY_TYPE* __restrict__ range,
        unsigned int numBuckets, unsigned int* __restrict__ bucketOffset) {
#if UNIFORM
    extern __shared__ KEY_TYPE minBuffer[];
    KEY_TYPE* maxBuffer = minBuffer+blockDim.x;
    KEY_TYPE minimum = MAX_KEY;
    KEY_TYPE maximum = MIN_KEY;

    // Each thread calculates the range of a subset of values.

    for (unsigned int index = threadIdx.x; index < length; index += blockDim.x) {
        KEY_TYPE value = getValue(data[index]);
        minimum = min(minimum, value);
        maximum = max(maximum, value);
    }

    // Now reduce them.

    minBuffer[threadIdx.x] = minimum;
    maxBuffer[threadIdx.x] = maximum;
    __syncthreads();
    for (unsigned int step = 1; step < blockDim.x; step *= 2) {
        if (threadIdx.x+step < blockDim.x && threadIdx.x%(2*step) == 0) {
            minBuffer[threadIdx.x] = min(minBuffer[threadIdx.x], minBuffer[threadIdx.x+step]);
            maxBuffer[threadIdx.x] = max(maxBuffer[threadIdx.x], maxBuffer[threadIdx.x+step]);
        }
        __syncthreads();
    }
    minimum = minBuffer[0];
    maximum = maxBuffer[0];
    if (threadIdx.x == 0) {
        range[0] = minimum;
        range[1] = maximum;
    }
#endif

    // Clear the bucket counters in preparation for the next kernel.

    for (unsigned int index = threadIdx.x; index < numBuckets; index += blockDim.x)
        bucketOffset[index] = 0;
}

/**
 * @brief Bucket phase 2 (uniform variant): assign each element to a bucket by a linear map.
 *
 * Phase role: bucketing. Selected by kernel NAME when UNIFORM==1. Buckets tile the [min,max]
 * key range with uniform width; each element maps to floor((key-min)/width), clamped to the
 * last bucket, so the map is key-monotone. On return each element has recorded its bucket
 * (@p bucketOfElement) and a unique reserved slot within that bucket (@p offsetInBucket),
 * and @p bucketOffset holds the final per-bucket counts. Guarantee: over all elements the
 * (bucket, slot) pairs are a contiguous, collision-free numbering per bucket, which is what
 * makes the phase-4 scatter race-free.
 *
 * @par Launch configuration
 *      Grid-strided over all elements; CudaSort::sort passes threads == dataLength,
 *      blockSize == 128. No shared memory, no cross-thread synchronization.
 *
 * @param[in]     data           Source array, read-only, length @p length. Borrowed.
 * @param[in]     length         Number of elements.
 * @param[in]     numBuckets     Bucket count.
 * @param[in]     range          Two-element {min,max} key produced by computeRange.
 * @param[in,out] bucketOffset   numBuckets-element counters; atomically incremented, ending
 *                               as per-bucket sizes. @pre pre-zeroed by computeRange.
 * @param[out]    bucketOfElement Per-element chosen bucket, indexed in input order, length @p length.
 * @param[out]    offsetInBucket Per-element reserved slot within its bucket, indexed in input
 *                               order, length @p length.
 *
 * @note Keys are cast to float for the mapping arithmetic regardless of KEY_TYPE; single
 *       precision affects only bucket placement, not final order (phase 5 re-sorts).
 */
__global__ void assignElementsToBuckets(const DATA_TYPE* __restrict__ data, unsigned int length, unsigned int numBuckets, const KEY_TYPE* __restrict__ range,
        unsigned int* __restrict__ bucketOffset, unsigned int* __restrict__ bucketOfElement, unsigned int* __restrict__ offsetInBucket) {
    // The whole key->bucket map runs in float regardless of KEY_TYPE. Misplacement only
    // unbalances buckets; sortBuckets re-sorts each bucket on the true KEY_TYPE, so order stays exact.
    float minValue = (float) (range[0]);
    float maxValue = (float) (range[1]);
    float bucketWidth = (maxValue-minValue)/numBuckets;
    for (unsigned int index = blockDim.x*blockIdx.x+threadIdx.x; index < length; index += blockDim.x*gridDim.x) {
        float key = (float) getValue(data[index]);
        unsigned int bucketIndex = min((unsigned int) ((key-minValue)/bucketWidth), numBuckets-1);
        offsetInBucket[index] = atomicAdd(&bucketOffset[bucketIndex], 1);
        bucketOfElement[index] = bucketIndex;
    }
}

/**
 * @brief Bucket phase 2 (non-uniform variant): assign elements via an adaptive, sample-
 *        estimated piecewise-linear key->bucket map.
 *
 * Phase role: bucketing for skewed key distributions. Selected by kernel NAME when
 * UNIFORM==0 (the Robosample neighbor-list case). A 64-element sample of the array is
 * sorted, and its order statistics define a key-monotone 8-segment approximate-CDF map that
 * hands each segment an equal share of the bucket range, so expected occupancy is
 * equalized and no small set of huge buckets can overflow the shared-memory bucket sort.
 * Output contract is identical to the uniform variant: each element records a bucket and a
 * unique reserved slot, and @p bucketOffset ends holding per-bucket counts, giving a
 * collision-free per-bucket numbering.
 *
 * @par Launch configuration
 *      CudaSort::sort passes threads == dataLength, blockSize == 128, so @pre blockDim.x
 *      >= 64 (the sample load and sort use the first 64 threads). Grid-strided over all
 *      elements.
 * @par Shared memory
 *      Static only: a 64-key sample buffer plus three 9-entry float segment tables; no
 *      dynamic shared memory is passed.
 * @par Synchronization
 *      Block-wide __syncthreads() after the sample load, after every sample-sort stage, and
 *      after the thread-0 segment build; all threads of the block must participate.
 *
 * @param[in]     data           Source array, read-only, length @p length. Borrowed.
 * @param[in]     length         Number of elements.
 * @param[in]     numBuckets     Bucket count.
 * @param[in]     range          Ignored by this variant; present only to share the launch
 *                               signature with the uniform kernel.
 * @param[in,out] bucketOffset   numBuckets-element counters; atomically incremented, ending
 *                               as per-bucket sizes. @pre pre-zeroed by computeRange.
 * @param[out]    bucketOfElement Per-element chosen bucket, input-order indexed, length @p length.
 * @param[out]    offsetInBucket Per-element reserved slot within its bucket, input-order
 *                               indexed, length @p length.
 *
 * @note Keys are handled in float throughout; a degenerate (zero-width) segment collapses to
 *       its base bucket index. The map is monotone in the sample keys, preserving global sort
 *       correctness independent of the estimate's quality.
 */
__global__ void assignElementsToBuckets2(const DATA_TYPE* __restrict__ data, unsigned int length, unsigned int numBuckets, const KEY_TYPE* __restrict__ range,
        unsigned int* __restrict__ bucketOffset, unsigned int* __restrict__ bucketOfElement, unsigned int* __restrict__ offsetInBucket) {
    // Only one of assignElementsToBuckets / assignElementsToBuckets2 is ever launched: CudaSort
    // resolves the kernel by name from the uniform flag. Robosample's BlockSortTrait passes
    // uniform=false, so the neighbor-list sort always runs this non-uniform variant.

    // Load 64 datapoints and sort them to get an estimate of the data distribution.

    __shared__ KEY_TYPE elements[64];
    if (threadIdx.x < 64) {
        int index = (int) (threadIdx.x*length/64.0);
        elements[threadIdx.x] = getValue(data[index]);
    }
    __syncthreads();
    for (unsigned int k = 2; k <= 64; k *= 2) {
        for (unsigned int j = k/2; j > 0; j /= 2) {
            if (threadIdx.x < 64) {
                int ixj = threadIdx.x^j;
                if (ixj > threadIdx.x) {
                    KEY_TYPE value1 = elements[threadIdx.x];
                    KEY_TYPE value2 = elements[ixj];
                    bool ascending = (threadIdx.x&k) == 0;
                    KEY_TYPE lowKey = (ascending ? value1 : value2);
                    KEY_TYPE highKey = (ascending ? value2 : value1);
                    if (lowKey > highKey) {
                        elements[threadIdx.x] = value2;
                        elements[ixj] = value1;
                    }
                }
            }
            __syncthreads();
        }
    }

    // Create a function composed of linear segments mapping data values to bucket indices.

    __shared__ float segmentLowerBound[9];
    __shared__ float segmentBaseIndex[9];
    __shared__ float segmentIndexScale[9];
    if (threadIdx.x == 0) {
        segmentLowerBound[0] = elements[0]-0.2f*(elements[5]-elements[0]);
        segmentLowerBound[1] = elements[5];
        segmentLowerBound[2] = elements[10];
        segmentLowerBound[3] = elements[20];
        segmentLowerBound[4] = elements[30];
        segmentLowerBound[5] = elements[40];
        segmentLowerBound[6] = elements[50];
        segmentLowerBound[7] = elements[60];
        segmentLowerBound[8] = elements[63]+0.2f*(elements[63]-elements[58]);
        segmentBaseIndex[0] = numBuckets/16;
        segmentBaseIndex[1] = 3*numBuckets/16;
        segmentBaseIndex[2] = 5*numBuckets/16;
        segmentBaseIndex[3] = 7*numBuckets/16;
        segmentBaseIndex[4] = 9*numBuckets/16;
        segmentBaseIndex[5] = 11*numBuckets/16;
        segmentBaseIndex[6] = 13*numBuckets/16;
        segmentBaseIndex[7] = 15*numBuckets/16;
        segmentBaseIndex[8] = numBuckets;
        for (int i = 0; i < 8; i++)
            if (segmentLowerBound[i+1] == segmentLowerBound[i])
                segmentIndexScale[i] = 0;
            else
                segmentIndexScale[i] = (segmentBaseIndex[i+1]-segmentBaseIndex[i])/(segmentLowerBound[i+1]-segmentLowerBound[i]);
    }
    __syncthreads();

    // Assign elements to buckets.

    for (unsigned int index = blockDim.x*blockIdx.x+threadIdx.x; index < length; index += blockDim.x*gridDim.x) {
        float key = (float) getValue(data[index]);
        int segment;
        for (segment = 0; segment < 7 && key > segmentLowerBound[segment+1]; segment++)
            ;
        unsigned int bucketIndex = segmentBaseIndex[segment]+(key-segmentLowerBound[segment])*segmentIndexScale[segment];
        bucketIndex = min(max(0, bucketIndex), numBuckets-1);
        // This atomic race is the sort's ONLY run-to-run nondeterminism: it fixes which slot
        // equal-keyed elements take. Final key order is unaffected; for BlockSortTrait the key is the
        // element itself (all keys distinct), so the neighbor-list result is fully deterministic.
        offsetInBucket[index] = atomicAdd(&bucketOffset[bucketIndex], 1);
        bucketOfElement[index] = bucketIndex;
    }
}

/**
 * @brief Bucket phase 3: inclusive prefix scan of bucket sizes into cumulative end offsets.
 *
 * Phase role: prefix scan. Converts @p bucketOffset in place from per-bucket counts (phase 2)
 * to an INCLUSIVE scan. On return bucketOffset[i] is the one-past-the-end offset of bucket i,
 * equivalently the start offset of bucket i+1. Downstream consumers therefore read bucket i
 * as the half-open range [ (i==0 ? 0 : bucketOffset[i-1]), bucketOffset[i] ) — the start of a
 * bucket is the PREDECESSOR's entry, not its own (see copyDataToBuckets and sortBuckets). The
 * last entry equals the total element count.
 *
 * @par Launch configuration
 *      Exactly one block (CudaSort::sort passes threads == blockSize == positionsKernelSize).
 *      numBuckets may exceed blockDim.x, so buckets are scanned in tiles of blockDim.x with a
 *      running carry across tiles. @pre blockDim.x is a power of two (required by the scan).
 * @par Shared memory
 *      Dynamic; the driver passes positionsKernelSize*sizeof(int) bytes, one unsigned int per
 *      thread. @pre available dynamic shared memory >= blockDim.x*sizeof(unsigned int).
 * @par Synchronization
 *      Block-wide __syncthreads() around each tile and each scan step; all threads
 *      participate. Single block, no cross-block sync.
 *
 * @param[in]     numBuckets   Number of buckets.
 * @param[in,out] bucketOffset numBuckets-element array; per-bucket counts in, inclusive-scan
 *                             end offsets out.
 */
__global__ void computeBucketPositions(unsigned int numBuckets, unsigned int* __restrict__ bucketOffset) {
    extern __shared__ unsigned int posBuffer[];
    unsigned int globalOffset = 0;
    for (unsigned int startBucket = 0; startBucket < numBuckets; startBucket += blockDim.x) {
        // Load the bucket sizes into local memory.

        unsigned int globalIndex = startBucket+threadIdx.x;
        __syncthreads();
        posBuffer[threadIdx.x] = (globalIndex < numBuckets ? bucketOffset[globalIndex] : 0);
        __syncthreads();

        // Perform a parallel prefix sum.

        for (unsigned int step = 1; step < blockDim.x; step *= 2) {
            unsigned int add = (threadIdx.x >= step ? posBuffer[threadIdx.x-step] : 0);
            __syncthreads();
            posBuffer[threadIdx.x] += add;
            __syncthreads();
        }

        // Write the results back to global memory.

        // INCLUSIVE scan: entry i becomes bucket i's END offset (== bucket i+1's start).
        // copyDataToBuckets and sortBuckets recover bucket i's start as bucketOffset[i-1] (0 for i==0).
        if (globalIndex < numBuckets)
            bucketOffset[globalIndex] = posBuffer[threadIdx.x]+globalOffset;
        globalOffset += posBuffer[blockDim.x-1];
    }
}

/**
 * @brief Bucket phase 4: scatter each element into its bucket's contiguous region.
 *
 * Phase role: scatter. Each element is written to buckets[start + slot], where start is its
 * bucket's begin offset ( bucketOffset[b-1], or 0 for bucket 0, per the inclusive scan) and
 * slot is the reserved offsetInBucket from phase 2. Guarantee: on return @p buckets is
 * bucket-major — buckets are contiguous and in ascending bucket order — but elements within a
 * bucket remain unsorted (phase 5 sorts them). Every (bucket, slot) pair is unique by
 * construction, so no two threads write the same slot; the scatter is race-free without
 * atomics.
 *
 * @par Launch configuration
 *      Grid-strided over all elements; CudaSort::sort passes threads == dataLength, blockSize
 *      defaulting to ThreadBlockSize == 64. No shared memory, no synchronization.
 *
 * @param[in]  data            Source array, read-only, length @p length. Borrowed.
 * @param[out] buckets         Bucket-contiguous scratch, length @p length; each slot written once.
 * @param[in]  length          Number of elements.
 * @param[in]  bucketOffset    Inclusive end offsets from phase 3 (numBuckets entries).
 * @param[in]  bucketOfElement Per-element bucket from phase 2, length @p length.
 * @param[in]  offsetInBucket  Per-element reserved slot from phase 2, length @p length.
 */
__global__ void copyDataToBuckets(const DATA_TYPE* __restrict__ data, DATA_TYPE* __restrict__ buckets, unsigned int length, const unsigned int* __restrict__ bucketOffset, const unsigned int* __restrict__ bucketOfElement, const unsigned int* __restrict__ offsetInBucket) {
    for (unsigned int index = blockDim.x*blockIdx.x+threadIdx.x; index < length; index += blockDim.x*gridDim.x) {
        DATA_TYPE element = data[index];
        unsigned int bucketIndex = bucketOfElement[index];
        // bucketOffset is computeBucketPositions' INCLUSIVE scan, so the predecessor entry
        // [bucketIndex-1] is exactly this bucket's start offset.
        unsigned int offset = (bucketIndex == 0 ? 0 : bucketOffset[bucketIndex-1]);
        buckets[offset+offsetInBucket[index]] = element;
    }
}

/**
 * @brief Bucket phase 5: bitonic-sort each bucket independently and emit the sorted array.
 *
 * Phase role: final per-bucket sort. One block handles one bucket at a time; blocks
 * grid-stride over buckets. Bucket i spans [ (i==0 ? 0 : bucketOffset[i-1]), bucketOffset[i] )
 * in the phase-4 layout. On return the corresponding range of @p data is sorted by key;
 * since buckets are key-monotone, the concatenation is the globally sorted array. Two regimes
 * with identical output contract:
 *   - bucket length <= blockDim.x (the common, driver-targeted case): sorted in shared memory,
 *     the unused high slots padded with MAX_VALUE (key == MAX_KEY) so padding sorts to the top
 *     and is never written back.
 *   - bucket length > blockDim.x (overflow safety net): sorted in place in global memory with
 *     the length-generalized network, covering arbitrarily large buckets when the adaptive
 *     bucketing still overflows one past blockDim.x.
 *
 * @par Launch configuration
 *      CudaSort::sort sizes the grid to round dataLength up to a multiple of sortKernelSize
 *      and passes blockSize == sortKernelSize; the grid is capped at numThreadBlocks, so
 *      blocks grid-stride over buckets. @pre blockDim.x is a power of two (bitonic network).
 * @par Shared memory
 *      Dynamic; the driver passes sortKernelSize*getDataSize() bytes, one DATA_TYPE per
 *      thread. @pre available dynamic shared memory >= blockDim.x*sizeof(DATA_TYPE).
 * @par Synchronization
 *      Block-wide __syncthreads() between stages; the global-memory regime adds
 *      __threadfence_block() for cross-thread visibility of its in-place writes. All threads
 *      of the block must participate in every stage. No cross-block synchronization: buckets
 *      are disjoint, so distinct blocks never touch the same @p data range.
 *
 * @param[out] data         Final sorted output array; this kernel writes each bucket's range.
 * @param[in]  buckets      Bucket-contiguous source from phase 4, read-only. Borrowed.
 * @param[in]  numBuckets   Number of buckets.
 * @param[in]  bucketOffset Inclusive end offsets (bucket boundaries) from phase 3, numBuckets
 *                          entries.
 *
 * @pre MAX_VALUE is a valid DATA_TYPE whose key equals the global maximum key (MAX_KEY);
 *      the small-bucket regime relies on padding sorting strictly above real data.
 */
__global__ void sortBuckets(DATA_TYPE* __restrict__ data, const DATA_TYPE* __restrict__ buckets, unsigned int numBuckets, const unsigned int* __restrict__ bucketOffset) {
    extern __shared__ DATA_TYPE dataBuffer[];
    for (unsigned int index = blockIdx.x; index < numBuckets; index += gridDim.x) {
        unsigned int startIndex = (index == 0 ? 0 : bucketOffset[index-1]);
        unsigned int endIndex = bucketOffset[index];
        unsigned int length = endIndex-startIndex;
        if (length <= blockDim.x) {
            // Load the data into local memory.

            if (threadIdx.x < length)
                dataBuffer[threadIdx.x] = buckets[startIndex+threadIdx.x];
            else
                dataBuffer[threadIdx.x] = MAX_VALUE;
            __syncthreads();

            // Perform a bitonic sort in local memory.

            for (unsigned int k = 2; k <= blockDim.x; k *= 2) {
                for (unsigned int j = k/2; j > 0; j /= 2) {
                    int ixj = threadIdx.x^j;
                    if (ixj > threadIdx.x) {
                        DATA_TYPE value1 = dataBuffer[threadIdx.x];
                        DATA_TYPE value2 = dataBuffer[ixj];
                        bool ascending = (threadIdx.x&k) == 0;
                        KEY_TYPE lowKey = (ascending ? getValue(value1) : getValue(value2));
                        KEY_TYPE highKey = (ascending ? getValue(value2) : getValue(value1));
                        if (lowKey > highKey) {
                            dataBuffer[threadIdx.x] = value2;
                            dataBuffer[ixj] = value1;
                        }
                    }
                    __syncthreads();
                }
            }

            // Write the data to the sorted array.

            if (threadIdx.x < length)
                data[startIndex+threadIdx.x] = dataBuffer[threadIdx.x];
        }
        else {
            // Cold path: CudaSort targets ~blockDim.x/2 elements per bucket, so length>blockDim.x is a
            // rare overflow. This in-global-memory bitonic sort is a correctness safety net, not hot.
            // Copy the bucket data over to the output array.

            for (unsigned int i = threadIdx.x; i < length; i += blockDim.x)
                data[startIndex+i] = buckets[startIndex+i];
            __threadfence_block();
            __syncthreads();

            // Perform a bitonic sort in global memory.

            for (unsigned int k = 2; k < 2*length; k *= 2) {
                for (unsigned int j = k/2; j > 0; j /= 2) {
                    for (unsigned int i = threadIdx.x; i < length; i += blockDim.x) {
                        int ixj = i^j;
                        if (ixj > i && ixj < length) {
                            DATA_TYPE value1 = data[startIndex+i];
                            DATA_TYPE value2 = data[startIndex+ixj];
                            bool ascending = ((i&k) == 0);
                            for (unsigned int mask = k*2; mask < 2*length; mask *= 2)
                                ascending = ((i&mask) == 0 ? !ascending : ascending);
                            KEY_TYPE lowKey  = (ascending ? getValue(value1) : getValue(value2));
                            KEY_TYPE highKey = (ascending ? getValue(value2) : getValue(value1));
                            if (lowKey > highKey) {
                                data[startIndex+i] = value2;
                                data[startIndex+ixj] = value1;
                            }
                        }
                    }
                    __threadfence_block();
                    __syncthreads();
                }
            }
        }
    }
}

}