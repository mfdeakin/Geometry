
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

/* Warning: Requires fptype to be recognized by genericfp */
template <typename fptype>
QuadType classifyQuadric(
    const Geometry::Quadric<3, fptype> &quad) {
  constexpr const int precision =
      fpconvert<fptype>::precision;
  mpfr_t qc[quad.numCoeffs];
  int err = 0;
  for(int i = 0; i < quad.numCoeffs; i++) {
    mpfr_init2(qc[i], precision);
    err += mpfr_set_flt(qc[i], quad.currentCoeffs[i],
                        MPFR_RNDN);
  }


  for(int i = 0; i < quad.numCoeffs; i++) mpfr_clear(qc[i]);
  if(err) return QUADT_ERROR;
  return QUADT_ERROR;
}
};

#endif
