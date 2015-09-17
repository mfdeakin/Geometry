
#include <gtest/gtest.h>

#include "quadrics.hpp"
#include "origin.hpp"
#include "accurate_math.hpp"

TEST(Quadric, SimpleClassify) {
  constexpr int dim = 3;
  constexpr const unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;
  typedef float fptype;
  struct teststruct {
    AccurateMath::QuadType expected;
    fptype coeffs[numCoeffs];
  } tests[] = {
      {AccurateMath::QUADT_HYPERBOLOID_ONE,
       {1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0}},
      {AccurateMath::QUADT_HYPERBOLOID_TWO,
       {1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {AccurateMath::QUADT_PARABOLOID_ELL,
       {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}},
      {AccurateMath::QUADT_PARABOLOID_HYP,
       {1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}},

      {AccurateMath::QUADT_ELLIPSOID_RE,
       {1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {AccurateMath::QUADT_ELLIPSOID_IM,
       {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {AccurateMath::QUADT_CONE_IM,
       {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {AccurateMath::QUADT_CONE_RE,
       {1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {AccurateMath::QUADT_CYLINDER_PAR,
       {0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0}},
      {AccurateMath::QUADT_CYLINDER_HYP,
       {1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      /* This one should be imaginary,
       * but the identification isn't done yet */
      {AccurateMath::QUADT_CYLINDER_ELL,
       {1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      /* This one should be real,
       * but the identification isn't done yet */
      {AccurateMath::QUADT_CYLINDER_ELL,
       {1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      /* This one should be imaginary,
       * but the identification isn't done yet */
      {AccurateMath::QUADT_PARALLELPLANES,
       {1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      /* This one should be real,
       * but the identification isn't done yet */
      {AccurateMath::QUADT_PARALLELPLANES,
       {1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {AccurateMath::QUADT_INTERSECTPLANES_IM,
       {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {AccurateMath::QUADT_INTERSECTPLANES_RE,
       {1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},

      {AccurateMath::QUADT_COINCIDENTPLANES,
       {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
  };
  Geometry::Origin<dim, fptype> o;
  Geometry::Quadric<dim, fptype> q(o);
  for(auto t : tests) {
    for(unsigned i = 0; i < numCoeffs; i++) {
      q.coeff(i) = t.coeffs[i];
    }
    auto quadtype = AccurateMath::classifyQuadric(q);
    printf(
        "Expected Quadric Type: %s\n"
        "Returned Quadric Type: %s\n",
        AccurateMath::QuadTypeNames[t.expected],
        AccurateMath::QuadTypeNames[quadtype]);
    EXPECT_EQ(t.expected, quadtype);
  }
}

TEST(Quadric, HyperboloidClassify) {
  constexpr int dim = 3;
  constexpr const unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;
  typedef float fptype;
  struct teststruct {
    AccurateMath::QuadType expected;
    fptype coeffs[numCoeffs];
  } tests[] = {
      {AccurateMath::QUADT_HYPERBOLOID_ONE,
       {1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0}},
      {AccurateMath::QUADT_HYPERBOLOID_TWO,
       {1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
      {AccurateMath::QUADT_HYPERBOLOID_TWO,
       {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0,
        29.0}},
      {AccurateMath::QUADT_HYPERBOLOID_TWO,
       {2.04395023038941230182, 3.1249825792834792,
        5.00483590573495728934729, 7.1293812379432,
        11.423897932984, 13.8579460935, 17.12983798173429,
        19.59378466987, 23.23947828357, 29.239478239748}},
      {AccurateMath::QUADT_HYPERBOLOID_ONE,
       {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0,
        -29.0}},
      {AccurateMath::QUADT_HYPERBOLOID_ONE,
       {2.04395023038941230182, 3.1249825792834792,
        5.00483590573495728934729, 7.1293812379432,
        11.423897932984, 13.8579460935, 17.12983798173429,
        19.59378466987, 23.23947828357, -29.239478239748}},
  };
  Geometry::Origin<dim, fptype> o;
  Geometry::Quadric<dim, fptype> q(o);
  for(auto t : tests) {
    for(unsigned i = 0; i < numCoeffs; i++) {
      q.coeff(i) = t.coeffs[i];
    }
    auto quadtype = AccurateMath::classifyQuadric(q);
    printf(
        "Expected Quadric Type: %s\n"
        "Returned Quadric Type: %s\n",
        AccurateMath::QuadTypeNames[t.expected],
        AccurateMath::QuadTypeNames[quadtype]);
    EXPECT_EQ(t.expected, quadtype);
  }
}
