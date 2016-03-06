
#include <gtest/gtest.h>

#include "quadric.hpp"
#include "origin.hpp"
#include "quadric_classify.hpp"
#include "accurate_math.hpp"

TEST(Quadric, SimpleClassify) {
  constexpr int dim = 3;
  constexpr const unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;
  typedef float fptype;
  struct teststruct {
    QuadricClassify::QuadType expected;
    fptype coeffs[numCoeffs];
  } tests[] = {
      {QuadricClassify::QUADT_HYPERBOLOID_ONE,
       {1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0}},
      {QuadricClassify::QUADT_HYPERBOLOID_TWO,
       {1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {QuadricClassify::QUADT_PARABOLOID_ELL,
       {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}},
      {QuadricClassify::QUADT_PARABOLOID_HYP,
       {1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}},

      {QuadricClassify::QUADT_ELLIPSOID_RE,
       {1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {QuadricClassify::QUADT_ELLIPSOID_IM,
       {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {QuadricClassify::QUADT_CONE_IM,
       {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {QuadricClassify::QUADT_CONE_RE,
       {1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {QuadricClassify::QUADT_CYLINDER_PAR,
       {0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0}},
      {QuadricClassify::QUADT_CYLINDER_HYP,
       {1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      /* This one should be imaginary,
       * but the identification isn't done yet */
      {QuadricClassify::QUADT_CYLINDER_ELL,
       {1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      /* This one should be real,
       * but the identification isn't done yet */
      {QuadricClassify::QUADT_CYLINDER_ELL,
       {1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      /* This one should be imaginary,
       * but the identification isn't done yet */
      {QuadricClassify::QUADT_PARALLELPLANES,
       {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      /* This one should be real,
       * but the identification isn't done yet */
      {QuadricClassify::QUADT_PARALLELPLANES,
       {1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {QuadricClassify::QUADT_INTERSECTPLANES_IM,
       {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {QuadricClassify::QUADT_INTERSECTPLANES_RE,
       {1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {QuadricClassify::QUADT_COINCIDENTPLANES,
       {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      /* This one should be a real plane */
      {QuadricClassify::QUADT_PLANE,
       {0, 0, 0, 1, 0, 0, 0, 0, 0, 1}},
  };
  Geometry::Origin<dim, fptype> o;
  Geometry::Quadric<dim, fptype> q(o);
  for(auto t : tests) {
    for(unsigned i = 0; i < numCoeffs; i++) {
      q.setCoeff(i, t.coeffs[i]);
    }
    auto quadtype = QuadricClassify::classifyQuadric(q);
    EXPECT_EQ(t.expected, quadtype);
  }
}

TEST(Quadric, HyperboloidClassify) {
  constexpr int dim = 3;
  constexpr const unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;
  typedef float fptype;
  struct teststruct {
    QuadricClassify::QuadType expected;
    fptype coeffs[numCoeffs];
  } tests[] = {
      {QuadricClassify::QUADT_HYPERBOLOID_ONE,
       {1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0}},
      {QuadricClassify::QUADT_HYPERBOLOID_TWO,
       {1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {QuadricClassify::QUADT_HYPERBOLOID_TWO,
       {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0,
        29.0}},
      {QuadricClassify::QUADT_HYPERBOLOID_TWO,
       {2.04395023038941230182, 3.1249825792834792,
        5.00483590573495728934729, 7.1293812379432,
        11.423897932984, 13.8579460935, 17.12983798173429,
        19.59378466987, 23.23947828357, 29.239478239748}},
      {QuadricClassify::QUADT_HYPERBOLOID_ONE,
       {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0,
        -29.0}},
      {QuadricClassify::QUADT_HYPERBOLOID_ONE,
       {2.04395023038941230182, 3.1249825792834792,
        5.00483590573495728934729, 7.1293812379432,
        11.423897932984, 13.8579460935, 17.12983798173429,
        19.59378466987, 23.23947828357, -29.239478239748}},
  };
  Geometry::Origin<dim, fptype> o;
  Geometry::Quadric<dim, fptype> q(o);
  for(auto t : tests) {
    for(unsigned i = 0; i < numCoeffs; i++) {
      q.setCoeff(i, t.coeffs[i]);
    }
    auto quadtype = QuadricClassify::classifyQuadric(q);
    EXPECT_EQ(t.expected, quadtype);
  }
}

TEST(Quadric, EllipsoidClassify) {
  constexpr const int dim = 3;
  // The following test fails for single precision floats
  typedef double fptype;
  constexpr const unsigned numCoeffs =
      Geometry::Quadric<dim, fptype>::numCoeffs;
  struct teststruct {
    QuadricClassify::QuadType expected;
    int detSign;
    int r3, r4;
    int eigenSameSign;
    fptype coeffs[numCoeffs];
  } tests[] = {
      QuadricClassify::QUADT_ELLIPSOID_RE,
      -1,
      3,
      4,
      1,
      {12.4697265625, 93139.4, 3567848359.0322265625,
       891985376.72802734375, 0.0, 0.0,
       2.0 * -6.23486328125, 0.0, 2.0 * -46569.705078125,
       2.0 * -1783924179.51611328125}};
  Geometry::Origin<dim, fptype> o;
  Geometry::Quadric<dim, fptype> q(o);
  for(auto t : tests) {
    for(unsigned i = 0; i < numCoeffs; i++) {
      q.setCoeff(i, t.coeffs[i]);
    }
    int detSign = QuadricClassify::classifyCalcDetSign(q);
    EXPECT_EQ(t.detSign, detSign);
    Array<int, 2> ranks =
        QuadricClassify::classifyCalcRank(q);
    EXPECT_EQ(t.r3, ranks[0]);
    EXPECT_EQ(t.r4, ranks[1]);
    int eigenSign =
        QuadricClassify::classifyCalcEigenSign(q);
    EXPECT_EQ(t.eigenSameSign, eigenSign);
    auto quadtype = QuadricClassify::classifyQuadric(q);
    EXPECT_EQ(t.expected, quadtype);
  }
}
