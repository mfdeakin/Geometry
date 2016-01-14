
#include "cudadef.h"

#include "vector.hpp"

#include <random>

#include <gtest/gtest.h>

constexpr static const int dim = 3;
using fptype = float;

using V = Geometry::Vector<dim, fptype>;

TEST(CudaCopy, Vector) {
	constexpr const fptype minValue = -(1 << 30);
	constexpr const fptype maxValue = 1 << 30;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
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
