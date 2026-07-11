/**
 * This file contains CUDA definitions for the macros and functions needed for the
 * common compute framework.
 */

#define KERNEL extern "C" __global__
#define DEVICE __device__
#define LOCAL __shared__
#define LOCAL_ARG
#define GLOBAL
#define RESTRICT __restrict__
#define LOCAL_ID threadIdx.x
#define LOCAL_SIZE blockDim.x
#define GLOBAL_ID (blockIdx.x*blockDim.x+threadIdx.x)
#define GLOBAL_SIZE (blockDim.x*gridDim.x)
#define GROUP_ID blockIdx.x
#define NUM_GROUPS gridDim.x
#define SYNC_THREADS __syncthreads();
#define MEM_FENCE __threadfence_block();
#define ATOMIC_ADD(dest, value) atomicAdd(dest, value)
#define FLT_MAX 3.40282347e+38f

using mm_long = long long;
using mm_ulong = unsigned long long;

#define SUPPORTS_64_BIT_ATOMICS 1
#define SUPPORTS_DOUBLE_PRECISION 1

__device__ inline auto realToFixedPoint(real value) -> long long {
    return static_cast<long long>(value * 0x100000000);
}
