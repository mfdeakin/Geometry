
#include "cudadef.h"

#include "origin.hpp"
#include "vector.hpp"
#include "point.hpp"

#include <random>

#include <gtest/gtest.h>

constexpr static const int dim = 3;
using fptype = float;

using A = Array<fptype, dim>;
using O = Geometry::Origin<dim, fptype>;
using V = Geometry::Vector<dim, fptype>;
using P = Geometry::Point<dim, fptype>;

TEST(Point, CudaCopy) {
	constexpr const fptype minValue = -(1 << 30);
	constexpr const fptype maxValue = 1 << 30;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
	std::random_device rd;
	std::mt19937_64 rng(rd());
	std::uniform_real_distribution<fptype> pdf(minValue, maxValue);
	constexpr const int numTests = 1000;
	for(int t = 0; t < numTests; t++) {
		V offset;
		A origin;
		for(int i = 0; i < dim; i++) {
			offset.set(i, pdf(rng));
			origin[i] = pdf(rng);
		}
		P pSrc(O(origin), offset);
		{
			std::shared_ptr<P::PointData> ptr1 = pSrc.cudaCopy();
			P pDest1;
			pDest1.cudaRetrieve(ptr1);
			ASSERT_EQ(pDest1, pSrc);
			P pDest2;
			pDest2.cudaRetrieve(ptr1.get());
			ASSERT_EQ(pDest2, pSrc);
		}
		{
			P::PointData *ptr2 = NULL;
			cudaMalloc(&ptr2, sizeof(*ptr2));
			pSrc.cudaCopy(ptr2);
			P pDest3;
			pDest3.cudaRetrieve(ptr2);
			ASSERT_EQ(pDest3, pSrc);
			cudaFree(ptr2);
		}
		{
			P::PointData *ptr3 = NULL;
			cudaMalloc(&ptr3, sizeof(*ptr3));
			std::shared_ptr<P::PointData> sptr3(ptr3, cudaFree);
			pSrc.cudaCopy(sptr3);
			P pDest4;
			pDest4.cudaRetrieve(sptr3);
			ASSERT_EQ(pDest4, pSrc);
			P pDest5;
			pDest5.cudaRetrieve(sptr3.get());
			ASSERT_EQ(pDest5, pSrc);
		}
	}
}
