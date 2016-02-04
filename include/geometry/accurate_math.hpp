
#ifndef _ACCURATE_MATH_HPP_
#define _ACCURATE_MATH_HPP_

#include <typeinfo>

#include <string.h>
#include <limits.h>

#include "cudadef.h"

#include "array.hpp"
#include "mathfuncs.hpp"
#include "genericfp.hpp"

namespace AccurateMath {

/* kahanSum
 * Implements Kahan Summation over a static array of
 * floating point values
 * Returns the sum
 */
template <unsigned size, typename fptype>
CUDA_CALLABLE static fptype kahanSum(
    const fptype(&summands)[size]) {
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

/* kahanSum
 * Implements Kahan Summation over a dynamic array of
 * floating point values.
 * Requires a valid pointer to the array
 * Returns the sum
 */
template <typename fptype>
CUDA_CALLABLE static fptype kahanSum(const fptype *summands,
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
CUDA_CALLABLE static fptype kahanSum(const fptype *summands,
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
CUDA_CALLABLE static fptype kahanSum(
    const fptype *summands) {
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
CUDA_CALLABLE static fptype twoNormSum(
    const fptype *summands, unsigned i, fptype curSum) {
  return MathFuncs::MathFuncs<fptype>::fma(
      summands[i], summands[i], curSum);
}

/* twoSum
 * An exact sum function
 * Returns an array of the floating point components,
 * with sum[0] > sum[1].
 * The exact value of the sum is sum[0] + sum[1].
 * "Robust Adaptive Floating-Point Geometric Predicates"
 * by Shewchuck
 */
template <typename fptype>
CUDA_CALLABLE static Array<fptype, 2> twoSum(fptype a,
                                             fptype b) {
  if(a < b) {
    fptype tmp = b;
    b = a;
    a = tmp;
  }
  fptype x = a + b;
  fptype c = x - a;
  fptype y = b - c;
  Array<fptype, 2> sum = {{x, y}};
  return sum;
}

/* twoProd
 * An exact product function
 * Returns an array of the floating point components,
 * with products[0] > products[1]
 * The exact value of the product is product[0] + product[1]
 * From "Some Functions Computable with a Fused-mac", by
 * Boldo, Muller
 */
template <typename fptype>
CUDA_CALLABLE static Array<fptype, 2> twoProd(fptype lhs,
                                              fptype rhs) {
  fptype prod = lhs * rhs;
  fptype err =
      MathFuncs::MathFuncs<fptype>::fma(lhs, rhs, -prod);
  Array<fptype, 2> products = {{prod, err}};
  return products;
}

/* threeFMA
 * An exact FMA function
 * Returns an array of the floating point components,
 * with fma[0] > fma[1] > fma[2]
 * The exact value of the FMA is fma[0] + fma[1] + fma[2]
 * From "Some Functions Computable with a Fused-mac", by
 * Boldo, Muller
 */
template <typename fptype>
CUDA_CALLABLE static Array<fptype, 3> threeFMA(fptype a,
                                               fptype b,
                                               fptype c) {
  fptype r1 = MathFuncs::MathFuncs<fptype>::fma(a, b, c);
  Array<fptype, 2> mult = twoProd(a, b);
  Array<fptype, 2> sum1 = twoSum(c, mult[1]);
  Array<fptype, 2> sum2 = twoSum(mult[0], sum1[0]);
  fptype gamma = (sum2[0] - r1) + sum2[1];
  Array<fptype, 2> sum3 = twoSum(gamma, sum1[1]);
  Array<fptype, 3> ret = {{r1, sum3[1], sum3[2]}};
  return ret;
}

/* compensatedDotProd
 * Computes the dot product of two vectors of equal
 * dimension
 * Requires valid arrays of floats as inputs
 * Uses threeFMA in combination with Kahan summation to
 * improve accuracy
 * From "Accurate dot products with FMA" by
 * Graillat, Langlois, and Louvet
 */
template <typename fptype>
CUDA_CALLABLE static fptype compensatedDotProd(
    const fptype *vec1, const fptype *vec2, unsigned dim) {
  Array<fptype, 2> prod = twoProd(vec1[0], vec2[0]);
  fptype s = prod[0];
  fptype c = prod[1];
  for(unsigned i = 1; i < dim; i++) {
    Array<fptype, 3> temp = threeFMA(vec1[i], vec2[i], s);
    s = temp[0];
    c = c + (temp[1] + temp[2]);
  }
  return s + c;
}

/* Solves the linear system linCoeff * t + constant = 0 */
template <typename fptype>
CUDA_CALLABLE static fptype solveLinear(fptype linCoeff,
                                        fptype constant) {
  return -constant / linCoeff;
}

/* Computes the discriminant at double the intial precision
 */
template <typename fptype>
CUDA_CALLABLE static fptype kahanDiscriminant(
    fptype sqCoeff, fptype linCoeff, fptype constant) {
  /* Computes the discriminant */
  Array<fptype, 2> prod = twoProd(sqCoeff, constant);
  Array<fptype, 2> sqr =
      twoProd(linCoeff, linCoeff / (fptype)4.0);
  fptype disc = sqr[1] - prod[1];
  disc += sqr[0] - prod[0];
  return disc;
}

/* Computes the roots of the quadratic equation
 * Uses double precision to avoid a catastrophic
 * loss of precision in the discriminant
 * Computes the root with a larger magnitude normally,
 * and the smaller root in a way that avoids
 * the catastrophic loss of precision.
 */
template <typename fptype>
CUDA_CALLABLE static Array<fptype, 2> kahanQuadratic(
    fptype sqCoeff, fptype linCoeff, fptype constant) {
  if(sqCoeff == 0.0) {
    return {solveLinear(linCoeff, constant), NAN};
  }
  fptype disc =
      kahanDiscriminant(sqCoeff, linCoeff, constant);
  if(disc < 0)
    return Array<fptype, 2>({{fptype(NAN), fptype(NAN)}});
  fptype fracPart = -MathFuncs::MathFuncs<fptype>::copysign(
      MathFuncs::MathFuncs<fptype>::fabs(linCoeff / 2.0) +
          MathFuncs::MathFuncs<fptype>::sqrt(disc),
      linCoeff);
  Array<fptype, 2> roots = {
      {fracPart / sqCoeff, constant / fracPart}};
  return roots;
}

/* Uses Newton's method to optimize a guess,
 * which may be roots computed with the standard quadratic
 * equation.
 * This is not subject to the same errors that the
 * regular root computation is,
 * and will usually bring the root to a more accurate value
 */
template <typename fptype>
CUDA_CALLABLE static fptype newtonsQuadratic(
    fptype sqCoeff, fptype linCoeff, fptype constant,
    fptype guess,
    fptype maxRelErr =
        4 * GenericFP::fpconvert<fptype>::epsilon) {
  fptype estimate = guess;
  fptype eval = MathFuncs::MathFuncs<fptype>::fma(
      MathFuncs::MathFuncs<fptype>::fma(sqCoeff, estimate,
                                        linCoeff),
      estimate, constant);
  while(MathFuncs::MathFuncs<fptype>::fabs(eval) >
        maxRelErr * eval) {
    fptype delta = MathFuncs::MathFuncs<fptype>::fma(
        2 * sqCoeff, estimate, linCoeff);
    if(delta == 0.0 || eval / delta == 0.0) break;
    estimate -= eval / delta;
    eval = MathFuncs::MathFuncs<fptype>::fma(
        MathFuncs::MathFuncs<fptype>::fma(sqCoeff, estimate,
                                          linCoeff),
        estimate, constant);
  }
  return estimate;
}
};

#endif
