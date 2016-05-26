
#ifndef _MATHFUNCS_HPP_
#define _MATHFUNCS_HPP_

#include <cmath>
#include "mpreal.hpp"

#include "genericfp.hpp"

#include "cudadef.h"

#ifdef __CUDA_ARCH__
#define CUDASTL(expression) expression
#else
#define CUDASTL(expression) std::expression
#endif

/* A stupid solution to an annoying limitation of C++
 * But a better one than giving up all the fast hardware
 * math functionality.
 * A better solution would be to write this in a code
 * generation language
 */
namespace MathFuncs {

template <typename fptype>
struct MathFuncs;

template <typename fptype, bool isStd>
struct IsStdMathFuncs;

template <typename fptype, bool isMPFR>
struct IsMPFRMathFuncs;

template <typename fptype>
struct IsStdMathFuncs<fptype, true> {
	CUDA_CALLABLE static fptype frexp(fptype v, int *exp) {
		return CUDASTL(frexp(v, exp));
	}
	CUDA_CALLABLE static fptype ldexp(fptype x, int exp) {
		return CUDASTL(ldexp(x, exp));
	}
  CUDA_CALLABLE static int getPrec(fptype v) {
    return GenericFP::fpconvert<fptype>::precision;
  }
  CUDA_CALLABLE static fptype abs(fptype v) {
    return CUDASTL(abs(v));
  }
  CUDA_CALLABLE static fptype fabs(fptype v) {
    return CUDASTL(fabs(v));
  }
  CUDA_CALLABLE static fptype min(fptype v1, fptype v2) {
#ifdef __CUDA_ARCH__
    return min(v1, v2);
#else
    return std::min(v1, v2);
#endif
  }
  CUDA_CALLABLE static fptype max(fptype v1, fptype v2) {
#ifdef __CUDA_ARCH__
    return max(v1, v2);
#else
    return std::max(v1, v2);
#endif
  }
  CUDA_CALLABLE static fptype fma(fptype m1, fptype m2,
                                  fptype a) {
#ifdef __CUDA_ARCH__
    return fma(m1, m2, a);
#else
    return std::fma(m1, m2, a);
#endif
  }
  CUDA_CALLABLE static bool signbit(fptype val) {
#ifdef __CUDA_ARCH__
    return signbit(val);
#else
    return std::signbit(val);
#endif
  }
  CUDA_CALLABLE static fptype copysign(fptype val,
                                       fptype sign) {
#ifdef __CUDA_ARCH__
    return copysign(val, sign);
#else
    return std::copysign(val, sign);
#endif
  }
  CUDA_CALLABLE static fptype sqr(fptype val) {
    return val * val;
  }
  CUDA_CALLABLE static fptype sqrt(fptype val) {
#ifdef __CUDA_ARCH__
    return sqrt(val);
#else
    return std::sqrt(val);
#endif
  }
  CUDA_CALLABLE static fptype isinf(fptype val) {
#ifdef __CUDA_ARCH__
    return isinf(val);
#else
    return std::isinf(val);
#endif
  }
  CUDA_CALLABLE static fptype isnan(fptype val) {
#ifdef __CUDA_ARCH__
    return isnan(val);
#else
    return std::isnan(val);
#endif
  }
  CUDA_CALLABLE static fptype sin(fptype val) {
#ifdef __CUDA_ARCH__
    return sin(val);
#else
    return std::sin(val);
#endif
  }

  CUDA_CALLABLE static fptype fastMult2Pow(fptype val,
                                           int twoPower) {
    union {
      fptype floatVal;
      GenericFP::fpconvert<fptype> structVal;
    } mult;
    mult.structVal.sign = 0;
    mult.structVal.mantissa = 0;
    mult.structVal.exponent =
        twoPower + GenericFP::fpconvert<fptype>::centralExp;
    return val * mult.floatVal;
  }
};

template <typename fptype>
struct IsStdMathFuncs<fptype, false> {};

template <typename fptype>
struct IsMPFRMathFuncs<fptype, true> {
  static int getPrec(fptype v) { return v.get_prec(); }
  static fptype abs(fptype v) { return mpfr::abs(v); }
  static fptype fabs(fptype v) { return mpfr::fabs(v); }
  static fptype min(fptype v1, fptype v2) {
    return mpfr::min(v1, v2);
  }
  static fptype max(fptype v1, fptype v2) {
    return mpfr::max(v1, v2);
  }
  static fptype fma(fptype m1, fptype m2, fptype a) {
    return mpfr::fma(m1, m2, a);
  }
  static bool signbit(fptype val) {
    return mpfr::signbit(val);
  }
  static fptype copysign(fptype val, fptype sign) {
    return mpfr::copysign(val, sign);
  }
  static fptype sqr(fptype val) { return mpfr::sqr(val); }
  static fptype sqrt(fptype val) { return mpfr::sqrt(val); }
  static fptype isinf(fptype val) {
    return mpfr::isinf(val);
  }
  static fptype isnan(fptype val) {
    return mpfr::isnan(val);
  }
  static fptype sin(fptype val) { return mpfr::sin(val); }
  static fptype fastMult2Pow(fptype val, int twoPower) {
    if(twoPower >= 0)
      return val << twoPower;
    else
      return val >> twoPower;
  }
};

template <typename fptype>
struct IsMPFRMathFuncs<fptype, false> {};

template <typename fptype>
struct MathFuncs
    : public IsStdMathFuncs<
          fptype,
          std::is_same<fptype, float>::value ||
              std::is_same<fptype, double>::value ||
              std::is_same<fptype, long double>::value ||
              std::is_same<fptype, long double>::value>,
      IsMPFRMathFuncs<
          fptype,
          std::is_same<fptype, mpfr::mpreal>::value> {};
};

#endif
