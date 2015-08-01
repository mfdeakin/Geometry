
#ifndef _ACCURATE_MATH_HPP_
#define _ACCURATE_MATH_HPP_

#include <array>
#include <cmath>
#include <limits.h>

#include "genericfp.hpp"

template <unsigned size, typename fptype>
fptype kahanSum(const fptype (&summands)[size]) {
  fptype ret = 0.0;
  fptype c = 0.0;
  for(unsigned i = 0; i < size; i++) {
    fptype mod = summands[i] - c;
    fptype tmp = ret + mod;
    c = (tmp - ret) - mod;
    ret = tmp;
  }
  return ret;
}

template <typename fptype>
fptype kahanSum(const fptype *summands, unsigned size) {
  fptype ret = 0.0;
  fptype c = 0.0;
  for(unsigned i = 0; i < size; i++) {
    fptype mod = summands[i] - c;
    fptype tmp = ret + mod;
    c = (tmp - ret) - mod;
    ret = tmp;
  }
  return ret;
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
std::array<fptype, 3> threeFMA(fptype a, fptype b, fptype c) {
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
    std::array<fptype, 3> temp = threeFMA(vec1[i], vec2[i], s);
    s = temp[0];
    c = c + (temp[1] + temp[2]);
  }
  return s + c;
}

#endif
