
#ifndef _ACCURATE_MATH_HPP_
#define _ACCURATE_MATH_HPP_

#include <array>
#include <cmath>
#include <typeinfo>

#include <limits.h>
#include <mpfr.h>

#include "geometry.hpp"
#include "quadrics.hpp"

#include "genericfp.hpp"

namespace AccurateMath {

template <unsigned size, typename fptype>
fptype kahanSum(const fptype(&summands)[size]) {
  fptype sum = 0.0;
  fptype extra = 0.0;
  for(unsigned i = 0; i < size; i++) {
    fptype mod = summands[i] - extra;
    fptype tmp = sum + mod;
    extra = (tmp - sum) - mod;
    sum = tmp;
  }
  return sum;
}

template <typename fptype>
fptype kahanSum(const fptype *summands, unsigned size) {
  fptype sum = 0.0;
  fptype extra = 0.0;
  for(unsigned i = 0; i < size; i++) {
    fptype mod = summands[i] - extra;
    fptype tmp = sum + mod;
    extra = (tmp - sum) - mod;
    sum = tmp;
  }
  return sum;
}

template <typename fptype,
          fptype (*sumFunc)(const fptype *summands,
                            unsigned i, fptype curSum)>
fptype kahanSum(const fptype *summands, unsigned size) {
  fptype sum = 0.0;
  fptype extra = 0.0;
  for(unsigned i = 0; i < size; i++) {
    fptype mod = summands[i] - extra;
    fptype tmp = sumFunc(summands, i, sum);
    extra = (tmp - sum) - mod;
    sum = tmp;
  }
  return sum;
}

template <typename fptype, unsigned size,
          fptype (*sumFunc)(const fptype *summands,
                            unsigned i, fptype curSum)>
fptype kahanSum(const fptype *summands) {
  fptype sum = 0.0;
  fptype extra = 0.0;
  for(unsigned i = 0; i < size; i++) {
    fptype mod = summands[i] - extra;
    fptype tmp = sumFunc(summands, i, sum);
    extra = (tmp - sum) - mod;
    sum = tmp;
  }
  return sum;
}
template <typename fptype>
fptype twoNormSum(const fptype *summands, unsigned i,
                  fptype curSum) {
  return std::fma(summands[i], summands[i], curSum);
}

template <typename fptype>
std::array<fptype, 2> twoSum(fptype a, fptype b) {
  if(a < b) {
    fptype tmp = b;
    b = a;
    a = tmp;
  }
  fptype x = a + b;
  fptype c = x - a;
  fptype y = b - c;
  std::array<fptype, 2> sum = {{x, y}};
  return sum;
}

template <typename fptype>
std::array<fptype, 2> twoProd(fptype lhs, fptype rhs) {
  fptype prod = lhs * rhs;
  fptype err = std::fma(lhs, rhs, -prod);
  std::array<fptype, 2> products = {{prod, err}};
  return products;
}

template <typename fptype>
std::array<fptype, 3> threeFMA(fptype a, fptype b,
                               fptype c) {
  fptype r1 = std::fma(a, b, c);
  std::array<fptype, 2> mult = twoProd(a, b);
  std::array<fptype, 2> sum1 = twoSum(c, mult[1]);
  std::array<fptype, 2> sum2 = twoSum(mult[0], sum1[0]);
  fptype gamma = (sum2[0] - r1) + sum2[1];
  std::array<fptype, 2> sum3 = twoSum(gamma, sum1[1]);
  std::array<fptype, 3> ret = {{r1, sum3[1], sum3[2]}};
  return ret;
}

template <typename fptype>
fptype compensatedDotProd(const fptype *vec1,
                          const fptype *vec2,
                          unsigned dim) {
  std::array<fptype, 2> prod = twoProd(vec1[0], vec2[0]);
  fptype s = prod[0];
  fptype c = prod[1];
  for(unsigned i = 1; i < dim; i++) {
    std::array<fptype, 3> temp =
        threeFMA(vec1[i], vec2[i], s);
    s = temp[0];
    c = c + (temp[1] + temp[2]);
  }
  return s + c;
}

enum QuadType {
  QUADT_COINCIDENTPLANES,
  QUADT_INTERSECTPLANES_IM,
  QUADT_INTERSECTPLANES_RE,
  QUADT_PARALLELPLANES_IM,
  QUADT_PARALLELPLANES_RE,
  QUADT_ELLIPSOID_IM,
  QUADT_ELLIPSOID_RE,
  QUADT_CONE_IM,
  QUADT_CONE_RE,
  QUADT_CYLINDER_ELL_IM,
  QUADT_CYLINDER_ELL_RE,
  QUADT_CYLINDER_HYP,
  QUADT_CYLINDER_PAR,
  QUADT_PARABOLOID_ELL,
  QUADT_PARABOLOID_HYP,
  QUADT_HYPERBOLOID_ONE,
  QUADT_HYPERBOLOID_TWO,
  QUADT_ERROR = -1
};

template <typename fptype>
int classifyCalcDet(
    const Geometry::Quadric<3, fptype> &quad, mpfr_t &det) {
  int err = 0;
  constexpr const int numDetTerms = 17;
  constexpr const int detCoeffs[] = {1,  -1, -1, -1, 2, -1,
                                     1,  2,  -2, -1, 1, 2,
                                     -2, 2,  2,  -1, 1};
  constexpr const int numDetProds = 4;
  constexpr const int detProds[numDetTerms][numDetProds] = {
      {0, 1, 2, 3},
      {2, 3, 4, 4},
      {1, 3, 5, 5},
      {1, 2, 6, 6},
      {3, 4, 5, 7},
      {0, 3, 7, 7},
      {6, 6, 7, 7},
      {2, 4, 6, 8},
      {5, 6, 7, 8},
      {0, 2, 8, 8},
      {5, 5, 8, 8},
      {1, 5, 6, 9},
      {4, 6, 7, 9},
      {4, 5, 8, 9},
      {0, 7, 8, 9},
      {0, 1, 9, 9},
      {4, 4, 9, 9}};
  constexpr const int precision =
      fpconvert<fptype>::precision;
  constexpr const int detTermPrec = precision * numDetProds;
  mpfr_t detTerm, detSum, tmpSum, extra, modAdd;
  mpfr_init2(detTerm, detTermPrec);
  mpfr_init2(detSum, 2 * detTermPrec);
  mpfr_init2(tmpSum, 2 * detTermPrec);
  mpfr_init2(extra, 2 * detTermPrec);
  mpfr_init2(modAdd, 2 * detTermPrec);
  mpfr_set_si(detSum, 0, MPFR_RNDN);
  for(int i = 0; i < numDetTerms; i++) {
    constexpr const int detTermPrec =
        numDetProds * precision;
    err = mpfr_set_si(detTerm, detCoeffs[i], MPFR_RNDN);
    if(err) return err;
    for(int j = 1; j < numDetProds; j++) {
      int coeffIdx = detProds[i][j];
      err = mpfr_mul_d(detTerm, detTerm,
                       quad.currentCoeffs[coeffIdx],
                       MPFR_RNDN);
      if(err) return err;
    }
    err = mpfr_sub(modAdd, detTerm, extra, MPFR_RNDN);
    if(err) return err;
    err = mpfr_add(tmpSum, modAdd, detSum, MPFR_RNDN);
    if(err) return err;
    err = mpfr_sub(extra, tmpSum, detSum, MPFR_RNDN);
    if(err) return err;
    err = mpfr_sub(extra, extra, modAdd, MPFR_RNDN);
    if(err) return err;
    mpfr_set(detSum, tmpSum, MPFR_RNDN);
  }
  mpfr_clear(detTerm);
  mpfr_clear(detSum);
  mpfr_clear(tmpSum);
  mpfr_clear(extra);
  mpfr_clear(modAdd);

  mpfr_set(det, detSum, MPFR_RNDN);
  return 0;
}

/* Warning: Requires fptype to be recognized by genericfp */
template <typename fptype>
QuadType classifyQuadric(
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}
};

#endif
