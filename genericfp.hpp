
#ifndef _GENFP_H_
#define _GENFP_H_

#include <assert.h>

/* The strucure for little endian architectures */
template<unsigned e, unsigned p> struct fp {
  unsigned long mantissa : p;
  unsigned long exponent : e;
  unsigned sign : 1;
  static const unsigned eBits = e, pBits = p;
} __attribute__((packed));

/* The bitfield lengths specified by IEEE 754 */
const unsigned gf16ExpBits = 4;
const unsigned gf16ManBits = 11;
const unsigned gf32ExpBits = 8;
const unsigned gf32ManBits = 23;
const unsigned gf43ExpBits = 11;
const unsigned gf43ManBits = 31;
const unsigned gf64ExpBits = 11;
const unsigned gf64ManBits = 52;
const unsigned gf79ExpBits = 15;
const unsigned gf79ManBits = 63;

/* Easy names for these structures */
typedef fp<gf16ExpBits, gf16ManBits> fp16;
typedef fp<gf32ExpBits, gf32ManBits> fp32;
typedef fp<gf43ExpBits, gf43ManBits> fp43;
typedef fp<gf64ExpBits, gf64ManBits> fp64;
typedef fp<gf79ExpBits, gf79ManBits> fp79;

template <typename fptype>
struct fpconvert;

template <>
struct fpconvert<float> : public fp32 {};
template <>
struct fpconvert<double> : public fp64 {};

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

template<unsigned pBits, unsigned eBits>
void gfFPToBinString(fp<pBits, eBits> in,
                     char (&out)[(pBits + eBits) / 4 + 1 + 6]);

template<unsigned pBits, unsigned eBits>
void gfFPToBinString(fp<pBits, eBits> in,
                     char (&out)[(pBits + eBits) / 4 + 1 + 6]) {
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
bool gfExpAllSet(T f)
{
  unsigned long long cmp = (1 << f.eBits) - 1;
  return cmp == f.exponent;
}

template <typename T>
bool gfManAllSet(T f)
{
  unsigned long long cmp = (1 << f.pBits) - 1;
  return cmp == f.mantissa;
}

template <typename T>
bool gfIsNaN(T f)
{
  return (gfExpAllSet<T>(f) == true) &&
    (f.mantissa != 0);
}

template <typename T>
bool gfIsInf(T f)
{
  return (gfExpAllSet<T>(f) == true) &&
    (f.mantissa == 0);
}

template <typename T>
bool gfGetMantissaBit(T f, unsigned bitPos)
{
  assert(bitPos < f.pBits);
  unsigned long selector = 1 << bitPos;
  unsigned long bit = f.mantissa & selector;
  bit >>= bitPos;
  return bit;
}

template <typename T>
bool gfGetExponentBit(T f, unsigned bitPos)
{
  assert(bitPos < f.eBits);
  unsigned long selector = 1 << bitPos;
  unsigned long bit = f.exponent & selector;
  bit >>= bitPos;
  return bit;
}

template <typename TDest, typename TSrc>
TDest gfRoundNearest(TSrc src)
{
  TDest dest;
  dest.sign = src.sign;
  dest.mantissa = 0;
  dest.exponent = 0;
  /* Compute the exponents corresponding to 1.0 */
  unsigned long srcCenter = (1 << src.eBits - 1) - 1;
  unsigned long destCenter = (1 << dest.eBits - 1) - 1;
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
      dest.mantissa |= 1 << dest.pBits - 1;
      dest.mantissa >>= numShifts;
    }
    /* Otherwise it's just 0 */
  }
  else {
    /* Plausibly not going to infinity :) */
    dest.exponent = src.exponent - centerDiff;
    if(dest.pBits >= src.pBits) {
      /* And we are done */
      dest.mantissa = src.mantissa;
    }
    else {
      unsigned roundingBit = src.pBits - dest.pBits;
      dest.mantissa = src.mantissa >> roundingBit;
      unsigned long truncated = src.mantissa &
        ((1 << roundingBit) - 1);
      /* Check the first truncated bit to see if we
       * need to consider rounding up */
      if((truncated & (1 << roundingBit - 1)) > 0) {
        unsigned long trailing;
        trailing = truncated &
          ((1 << (roundingBit - 1)) - 1);
        /* Round up if trailing is nonzero or if
         * it is zero whatever direction makes the
         * 0'th bit of the mantissa 0 */
        if(trailing > 0 ||
           ((dest.mantissa & 1) == 1)) {
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

#endif
