
#ifndef _MATHFUNCS_HPP_
#define _MATHFUNCS_HPP_

#include <cmath>
#include "mpreal.hpp"

#include "genericfp.hpp"

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
  static int getPrec(fptype v) {
    return GenericFP::fpconvert<fptype>::precision;
  }
  static fptype abs(fptype v) { return std::abs(v); }
  static fptype fabs(fptype v) { return std::fabs(v); }
  static fptype min(fptype v1, fptype v2) {
    return std::min(v1, v2);
  }
  static fptype max(fptype v1, fptype v2) {
    return std::max(v1, v2);
  }
  static fptype fma(fptype m1, fptype m2, fptype a) {
    return std::fma(m1, m2, a);
  }
  static bool signbit(fptype val) {
    return std::signbit(val);
  }
  static fptype copysign(fptype val, fptype sign) {
    return std::copysign(val, sign);
  }
  static fptype sqr(fptype val) { return val * val; }
  static fptype sqrt(fptype val) { return std::sqrt(val); }
  static fptype isinf(fptype val) {
    return std::isinf(val);
  }
  static fptype isnan(fptype val) {
    return std::isnan(val);
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
