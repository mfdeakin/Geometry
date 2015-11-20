
#ifndef _CUDADEF_H_
#define _CUDADEF_H_

#ifdef __CUDACC__
#include <cuda.h>
#define NDEBUG
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

#endif
