
#ifndef _GENFP_H_
#define _GENFP_H_

#include <cfloat>

#include <assert.h>

namespace GenericFP {

/* The strucure for little endian architectures */
template <unsigned e, unsigned p>
struct fp {
  unsigned long mantissa : p;
  unsigned long exponent : e;
  unsigned sign : 1;
  static constexpr const unsigned eBits = e, pBits = p,
                                  precision = p + 1;
  static constexpr const unsigned long centralExp =
      (1 << (e - 1)) - 1;
} __attribute__((packed));

/* The bitfield lengths specified by IEEE 754 */
constexpr const unsigned gf16ExpBits = 4;
constexpr const unsigned gf16ManBits = 11;
constexpr const unsigned gf32ExpBits = 8;
constexpr const unsigned gf32ManBits = 23;
constexpr const unsigned gf44ExpBits = 11;
constexpr const unsigned gf44ManBits = 32;
constexpr const unsigned gf64ExpBits = 11;
constexpr const unsigned gf64ManBits = 52;
constexpr const unsigned gf80ExpBits = 15;
constexpr const unsigned gf80ManBits = 64;

/* Easy names for these structures */
typedef fp<gf16ExpBits, gf16ManBits> fp16;
typedef fp<gf32ExpBits, gf32ManBits> fp32;
typedef fp<gf44ExpBits, gf44ManBits> fp44;
typedef fp<gf64ExpBits, gf64ManBits> fp64;
typedef fp<gf80ExpBits, gf80ManBits> fp79;

template <typename fptype>
struct fpconvert;

template <>
struct fpconvert<float> : public fp32 {
  static constexpr const float epsilon = FLT_EPSILON;
  static constexpr const char *fpname = "Single Precision";
};

template <>
struct fpconvert<double> : public fp64 {
  static constexpr const double epsilon = DBL_EPSILON;
  static constexpr const char *fpname = "Double Precision";
};

template <>
struct fpconvert<long double> : public fp79 {
  static constexpr const long double epsilon = LDBL_EPSILON;
  static constexpr const char *fpname =
      "Extended Double Precision";
};

template <typename fptype>
struct fpconvert<fptype> gfFPStruct(fptype value);
template <typename fptype>
fptype gfFPFloat(struct fpconvert<fptype> value);

template <typename T>
bool gfExpAllSet(T f);
template <typename T>
bool gfIsNaN(T f);
template <typename T>
bool gfIsInf(T f);
template <typename T>
bool gfManAllSet(T f);

template <typename T>
bool getMantissaBit(T f, unsigned bitPos);
template <typename T>
bool getExponentBit(T f, unsigned bitPos);
template <typename TDest, typename TSrc>
TDest gfRound(TSrc f);

template <typename T>
float gfFloat(T orig);
template <typename T>
double gfDouble(T orig);

template <unsigned pBits, unsigned eBits>
void gfFPToBinString(
    fp<pBits, eBits> in,
    char(&out)[(pBits + eBits) / 4 + 1 + 6]);

template <unsigned pBits, unsigned eBits>
void gfFPToBinString(
    fp<pBits, eBits> in,
    char(&out)[(pBits + eBits) / 4 + 1 + 6]) {
  snprintf(out, pBits + eBits + 5, "2^%x * 1.%x");
}

template <typename fptype>
struct fpconvert<fptype> gfFPStruct(fptype value) {
  union {
    struct fpconvert<fptype> ret;
    fptype value;
  } convert;
  convert.value = value;
  return convert.ret;
}

template <typename fptype>
fptype gfFPFloat(struct fpconvert<fptype> value) {
  union {
    struct fpconvert<fptype> value;
    fptype ret;
  } convert;
  convert.value = value;
  return convert.ret;
}

template <typename T>
bool gfExpAllSet(T f) {
  unsigned long long cmp = (1 << f.eBits) - 1;
  return cmp == f.exponent;
}

template <typename T>
bool gfManAllSet(T f) {
  unsigned long long cmp = (1 << f.pBits) - 1;
  return cmp == f.mantissa;
}

template <typename T>
bool gfIsNaN(T f) {
  return (gfExpAllSet<T>(f) == true) && (f.mantissa != 0);
}

template <typename T>
bool gfIsInf(T f) {
  return (gfExpAllSet<T>(f) == true) && (f.mantissa == 0);
}

template <typename T>
bool gfGetMantissaBit(T f, unsigned bitPos) {
  assert(bitPos < f.pBits);
  unsigned long selector = 1 << bitPos;
  unsigned long bit = f.mantissa & selector;
  bit >>= bitPos;
  return bit;
}

template <typename T>
bool gfGetExponentBit(T f, unsigned bitPos) {
  assert(bitPos < f.eBits);
  unsigned long selector = 1 << bitPos;
  unsigned long bit = f.exponent & selector;
  bit >>= bitPos;
  return bit;
}

template <typename TDest, typename TSrc>
TDest gfRoundNearest(TSrc src) {
  TDest dest;
  dest.sign = src.sign;
  dest.mantissa = 0;
  dest.exponent = 0;
  /* Compute the exponents corresponding to 1.0 */
  unsigned long srcCenter = (1 << (src.eBits - 1)) - 1;
  unsigned long destCenter = (1 << (dest.eBits - 1)) - 1;
  unsigned long centerDiff = srcCenter - destCenter;
  if(dest.eBits < src.eBits &&
     src.exponent >= 2 * destCenter + centerDiff) {
    /* Round it to infinity */
    dest.exponent = ~dest.exponent;
  }
  /* Verify it doesn't need to be denormalized */
  else if(dest.eBits < src.eBits &&
          src.exponent <= centerDiff) {
    unsigned long numShifts = centerDiff - src.exponent;
    if(numShifts < 8 * dest.pBits) {
      /* Translate the mantissa to a denormalized number */
      dest.mantissa = src.mantissa >> 1;
      dest.mantissa |= 1 << (dest.pBits - 1);
      dest.mantissa >>= numShifts;
    }
    /* Otherwise it's just 0 */
  } else {
    /* Plausibly not going to infinity :) */
    dest.exponent = src.exponent - centerDiff;
    if(dest.pBits >= src.pBits) {
      /* And we are done */
      dest.mantissa = src.mantissa;
    } else {
      unsigned roundingBit = src.pBits - dest.pBits;
      dest.mantissa = src.mantissa >> roundingBit;
      unsigned long truncated =
          src.mantissa & ((1 << roundingBit) - 1);
      /* Check the first truncated bit to see if we
       * need to consider rounding up */
      if((truncated & (1 << (roundingBit - 1))) > 0) {
        unsigned long trailing;
        trailing =
            truncated & ((1 << (roundingBit - 1)) - 1);
        /* Round up if trailing is nonzero or if
         * it is zero whatever direction makes the
         * 0'th bit of the mantissa 0 */
        if(trailing > 0 || ((dest.mantissa & 1) == 1)) {
          /* Round up. */
          dest.mantissa++;
          if(dest.mantissa == 0) {
            dest.exponent++;
          }
        }
      }
    }
  }
  return dest;
}
}

#endif
