
#ifndef _ACCURATE_MATH_HPP_
#define _ACCURATE_MATH_HPP_

#include <array>
#include <cmath>
#include <typeinfo>

#include <string.h>
#include <limits.h>

#include "genericfp.hpp"

namespace AccurateMath {

template <unsigned size, typename fptype>
static fptype kahanSum(const fptype(&summands)[size]) {
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
static fptype kahanSum(const fptype *summands,
                       unsigned size) {
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
static fptype kahanSum(const fptype *summands,
                       unsigned size) {
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
static fptype kahanSum(const fptype *summands) {
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
static fptype twoNormSum(const fptype *summands, unsigned i,
                         fptype curSum) {
  return std::fma(summands[i], summands[i], curSum);
}

template <typename fptype>
static std::array<fptype, 2> twoSum(fptype a, fptype b) {
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
static std::array<fptype, 2> twoProd(fptype lhs,
                                     fptype rhs) {
  fptype prod = lhs * rhs;
  fptype err = std::fma(lhs, rhs, -prod);
  std::array<fptype, 2> products = {{prod, err}};
  return products;
}

template <typename fptype>
static std::array<fptype, 3> threeFMA(fptype a, fptype b,
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
static fptype compensatedDotProd(const fptype *vec1,
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

template <typename fptype>
static fptype solveLinear(fptype linCoeff,
                          fptype constant) {
  return -constant / linCoeff;
}

template <typename fptype>
static fptype kahanDiscriminant(fptype sqCoeff,
                                fptype linCoeff,
                                fptype constant) {
  std::array<fptype, 2> prod = twoProd(sqCoeff, constant);
  std::array<fptype, 2> sqr =
      twoProd(linCoeff, linCoeff / (fptype)4.0);
  fptype disc = sqr[1] - prod[1];
  disc += sqr[0] - prod[0];
  return disc;
}

template <typename fptype>
static std::array<fptype, 2> kahanQuadratic(
    fptype sqCoeff, fptype linCoeff, fptype constant) {
  if(sqCoeff == 0.0) {
    return {solveLinear(linCoeff, constant), NAN};
  }
  fptype disc =
      kahanDiscriminant(sqCoeff, linCoeff, constant);
  if(disc < 0) return std::array<fptype, 2>({{NAN, NAN}});
  fptype fracPart =
      -std::copysign(
          std::abs(linCoeff / 2.0) + std::sqrt(disc),
          linCoeff);
  std::array<fptype, 2> roots = {
      {fracPart / sqCoeff, constant / fracPart}};
  return roots;
}
};

#endif
