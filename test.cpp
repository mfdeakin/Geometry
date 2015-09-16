
#include "origin.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "bsgtree.hpp"
#include "accurate_math.hpp"
#include <array>
#include <math.h>
#include <stdio.h>
#include <gtest/gtest.h>

TEST(Quadric, Classify) {
  constexpr int dim = 3;
  constexpr const unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;
  typedef float fptype;
  struct teststruct {
    AccurateMath::QuadType expected;
    fptype coeffs[numCoeffs];
  } tests[] = {{AccurateMath::QUADT_HYPERBOLOID_TWO,
                {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0,
                 23.0, 29.0}},
               {AccurateMath::QUADT_ELLIPSOID_RE,
                {1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0,
                 0.0, 0.0}}};
  Geometry::Origin<dim, fptype> o;
  Geometry::Quadric<dim, fptype> q(o);
  for(auto t : tests) {
    for(unsigned i = 0; i < numCoeffs; i++) {
      q.coeff(i) = t.coeffs[i];
    }
    auto quadtype = AccurateMath::classifyQuadric(q);
    EXPECT_EQ(t.expected, quadtype);
  }
}

int main(int argc, char **argv) {
  constexpr int dim = 3;
  typedef float fptype;
  Geometry::Origin<dim, fptype> o;
  Geometry::Vector<dim, fptype> ptOffset;
  Geometry::Point<dim, fptype> pt(o, ptOffset);
  Geometry::Vector<dim, fptype> dir;
  Geometry::Line<dim, fptype> line(pt, dir);
  Geometry::Quadric<dim, fptype> q(o);
  auto orthogs = dir.calcOrthogonals();
  q.coeff(0, 0) = 2.0;
  q.coeff(1, 1) = 3.0;
  q.coeff(2, 2) = 5.0;
  q.coeff(3, 3) = 7.0;

  q.coeff(0, 1) = 11.0;
  q.coeff(0, 2) = 13.0;
  q.coeff(0, 3) = 17.0;

  q.coeff(1, 2) = 19.0;
  q.coeff(1, 3) = 23.0;

  q.coeff(2, 3) = 29.0;
  auto quadtype = AccurateMath::classifyQuadric(q);
  assert(quadtype < AccurateMath::QUADT_ERRORINVALID);
  printf("Quadric Type: %d, %s\n", quadtype,
         AccurateMath::QuadTypeNames[quadtype]);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
