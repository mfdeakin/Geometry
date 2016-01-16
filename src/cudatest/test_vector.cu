
#include "cudadef.h"

#include "vector.hpp"

#include <random>

#include <gtest/gtest.h>

constexpr static const int dim = 3;
using fptype = float;

using V = Geometry::Vector<dim, fptype>;

TEST(Vector, CudaCopy) {
	constexpr const fptype minValue = -(1 << 30);
	constexpr const fptype maxValue = 1 << 30;
	std::random_device rd;
	std::mt19937_64 rng(rd());
	std::uniform_real_distribution<fptype> pdf(minValue, maxValue);
	constexpr const int numTests = 1000;
	for(int t = 0; t < numTests; t++) {
		V vSrc;
		for(int i = 0; i < dim; i++) {
			fptype coord = pdf(rng);
			vSrc.set(i, coord);
		}
		{
			std::shared_ptr<V::VectorData> ptr1 = vSrc.cudaCopy();
			V vDest1;
			vDest1.cudaRetrieve(ptr1);
			ASSERT_EQ(vDest1, vSrc);
			V vDest2;
			vDest2.cudaRetrieve(ptr1.get());
			ASSERT_EQ(vDest2, vSrc);
		}
		{
			V::VectorData *ptr2 = NULL;
			cudaMalloc(&ptr2, sizeof(*ptr2));
			vSrc.cudaCopy(ptr2);
			V vDest3;
			vDest3.cudaRetrieve(ptr2);
			ASSERT_EQ(vDest3, vSrc);
			cudaFree(ptr2);
		}
		{
			V::VectorData *ptr3 = NULL;
			cudaMalloc(&ptr3, sizeof(*ptr3));
			std::shared_ptr<V::VectorData> sptr3(ptr3, cudaFree);
			vSrc.cudaCopy(sptr3);
			V vDest4;
			vDest4.cudaRetrieve(sptr3);
			ASSERT_EQ(vDest4, vSrc);
			V vDest5;
			vDest5.cudaRetrieve(sptr3.get());
			ASSERT_EQ(vDest5, vSrc);
		}
	}
}

__global__ void CudaAdd(V::VectorData *lhs, V::VectorData *rhs,
                        V::VectorData *dest, int numVectors) {
	const int idx = threadIdx.x + blockIdx.x * blockDim.x;
	const int numThreads = blockDim.x * gridDim.x;
	const int vpt = numVectors / numThreads + 1;
	const int threadIdx = idx * vpt;
	for(int i = 0; i < vpt && threadIdx + i < numVectors; i++) {
		V vLHS(lhs[threadIdx + i]);
		V vRHS(rhs[threadIdx + i]);
		V sum = vLHS + vRHS;
		sum.copy(&dest[threadIdx + i]);
	}
}

TEST(Vector, CudaAdd) {
	constexpr const fptype minValue = -(1 << 30);
	constexpr const fptype maxValue = 1 << 30;
	std::random_device rd;
	std::mt19937_64 rng(rd());
	std::uniform_real_distribution<fptype> pdf(minValue, maxValue);
	constexpr const int numVecs = 100;
	constexpr const int numTests = 100;
	for(int t = 0; t < numTests; t++) {
		V lhs[numVecs];
		V rhs[numVecs];
		V::VectorData *lhsMem, *rhsMem, *destMem;
		cudaMalloc(&lhsMem, sizeof(V::VectorData[numVecs]));
		cudaMalloc(&rhsMem, sizeof(V::VectorData[numVecs]));
		cudaMalloc(&destMem, sizeof(V::VectorData[numVecs]));
		for(int i = 0; i < numVecs; i++) {
			for(int j = 0; j < dim; j++) {
				lhs[i].set(j, pdf(rng));
				rhs[i].set(j, pdf(rng));
			}
			lhs[i].cudaCopy(&lhsMem[i]);
			rhs[i].cudaCopy(&rhsMem[i]);
		}
		constexpr const int numBlocks = 32;
		constexpr const int numThreads = 32;
		CudaAdd<<<numBlocks, numThreads>>>(lhsMem, rhsMem, destMem, numVecs);
		for(int i = 0; i < numVecs; i++) {
			V sum = lhs[i] + rhs[i];
			V cmp;
			cmp.cudaRetrieve(&destMem[i]);
			ASSERT_EQ(sum, cmp);
		}
	}
}
