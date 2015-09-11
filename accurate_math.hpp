
#ifndef _ACCURATE_MATH_HPP_
#define _ACCURATE_MATH_HPP_

#include <array>
#include <cmath>
#include <typeinfo>

#include <string.h>
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
  QUADT_INTERSECTPLANES,
  QUADT_INTERSECTPLANES_IM,
  QUADT_INTERSECTPLANES_RE,
  QUADT_PARALLELPLANES,
  QUADT_PARALLELPLANES_IM,
  QUADT_PARALLELPLANES_RE,
  QUADT_ELLIPSOID,
  QUADT_ELLIPSOID_IM,
  QUADT_ELLIPSOID_RE,
  QUADT_CONE,
  QUADT_CONE_IM,
  QUADT_CONE_RE,
  QUADT_CYLINDER_ELL,
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
int classifyCalcDetSign(
    const Geometry::Quadric<3, fptype> &quad) {
  int err = 0;
  constexpr const int numDetTerms = 17;
  /* Only 4 error causing multiplications occur per term */
  constexpr const int numDetProds = 4;
  /* The determinant is as follows:
   * d = c0 c1 c2 c3 - c2 c3 c4^2 - c1 c3  c5^2 -
   *     c1 c2 c6^2 + 2 c3 c4 c5 c7 - c0 c3 c7^2 +
   *     c6^2 c7^2 + 2 c2 c4 c6 c8 - 2 c5 c6 c7 c8 -
   *     c0 c2 c8^2 + c5^2 c8^2 + 2 c1 c5 c6 c9 -
   *     2 c4 c6 c7 c9 + 2 c4 c5 c8 c9 + 2 c0 c7 c8 c9 -
   *     c0 c1 c9^2 + c4^2 c9^2
   */
  constexpr const int detCoeffs[] = {1,  -1, -1, -1, 2, -1,
                                     1,  2,  -2, -1, 1, 2,
                                     -2, 2,  2,  -1, 1};
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
  constexpr const int guessedExtraPrec = 2;
  mpfr_t detTerm, detSum, tmpSum, extra, modAdd;
  mpfr_init2(detTerm, detTermPrec);
  mpfr_init2(detSum, guessedExtraPrec * detTermPrec);
  mpfr_init2(tmpSum, guessedExtraPrec * detTermPrec);
  mpfr_init2(extra, guessedExtraPrec * detTermPrec);
  mpfr_init2(modAdd, guessedExtraPrec * detTermPrec);
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
  mpfr_clear(tmpSum);
  mpfr_clear(extra);
  mpfr_clear(modAdd);
  int cmpZero = mpfr_cmp_si(detSum, 0);
  mpfr_clear(detSum);
  return cmpZero;
}

template <unsigned mtxDim>
void eliminateColumn(mpfr_t elems[mtxDim][mtxDim],
                     int column, bool rowsDone[mtxDim],
                     mpfr_t &coeff, mpfr_t &delta) {
  for(unsigned i = 0; i < mtxDim; i++) {
    int isZero = mpfr_cmp_d(elems[i][column], 0.0);
    if(isZero != 0 && rowsDone[i] == false) {
      for(unsigned j = 0; j < mtxDim; j++) {
        if(j == i) continue;
        isZero = mpfr_cmp_d(elems[j][column], 0.0);
        if(isZero) continue;
        mpfr_div(coeff, elems[j][column], elems[i][column],
                 MPFR_RNDN);
        for(unsigned k = column + 1; k < mtxDim; k++) {
          mpfr_mul(delta, elems[i][k], coeff, MPFR_RNDN);
          mpfr_sub(elems[j][k], elems[j][k], delta,
                   MPFR_RNDN);
        }
        mpfr_set_d(elems[j][column], 0.0, MPFR_RNDN);
      }
      rowsDone[i] = true;
      break;
    }
  }
}

template <typename fptype>
std::array<int, 2> classifyCalcRank(
    const Geometry::Quadric<3, fptype> &quad) {
  constexpr const int mtxDim = 4;
  constexpr const int precision = 512;
  constexpr const unsigned mtxVals[mtxDim][mtxDim] = {
      {0, 4, 5, 6},
      {4, 1, 7, 8},
      {5, 7, 2, 9},
      {6, 8, 9, 3}};
  mpfr_t elems[mtxDim][mtxDim];
  /* Just use Gaussian Elimination with a ridiculous
   * amount
   * of precision */
  for(unsigned i = 0; i < mtxDim; i++) {
    for(unsigned j = 0; j < mtxDim; j++) {
      mpfr_init2(elems[i][j], precision);
      unsigned coeffNum = mtxVals[i][j];
      mpfr_set_d(elems[i][j], quad.currentCoeffs[coeffNum],
                 MPFR_RNDN);
    }
  }
  mpfr_t coeff, delta;
  mpfr_inits2(precision, coeff, delta, (mpfr_ptr)NULL);
  bool rowsDone[mtxDim];
  memset(rowsDone, false, sizeof(rowsDone));
  for(int numRowsDone = 0; numRowsDone < (mtxDim - 1);
      numRowsDone++)
    eliminateColumn<mtxDim>(elems, numRowsDone, rowsDone,
                            coeff, delta);
  mpfr_clears(coeff, delta, (mpfr_ptr)NULL);
  std::array<int, 2> ranks;
  ranks[0] = 0;
  ranks[1] = 0;
  for(int i = 0; i < mtxDim; i++) {
    for(int j = 0; j < mtxDim; j++) {
      int isZero = mpfr_cmp_d(elems[i][j], 0.0);
      if(isZero != 0) {
        ranks[1]++;
        if(i < (mtxDim - 1) && j < (mtxDim - 1)) ranks[0]++;
      }
      mpfr_clear(elems[i][j]);
    }
  }
  return ranks;
}

template <typename fptype>
int classifyCalcEigenSign(
    const Geometry::Quadric<3, fptype> &quad) {
  /* Basically I'm computing the roots of the derivative
   * of the characteristic cubic of the matrix.
   * This is easier, and with the constant term,
   * sufficient to determine the signs of the eigenvalues.
   *
   * The cubic is of the following form:
   * -x^3 + (c0 + c1 + c2)x^2 +
   * (c4^2/4 + c5^2/4 + c7^2/4 - c0 c1 - c0 c2 - c1 c2)x +
   * c0 c1 c2 - c0 c7^2/4 - c1 c5^2/4 -
   * c2 c4^2/4 + c4 c5 c7/4
   *
   * The derivative is of the following form:
   * -3 x^2 + 2(c0 + c1 + c2)x +
   * ((c4^2 + c5^2 + c7^2)/4 - c0 c1 - c0 c2 - c1 c2)x
   */
  constexpr const int precision =
      fpconvert<fptype>::precision;
  constexpr const int numProds = 2;
  constexpr const int guessedExtraPrec = 2;
  constexpr const int sqrTermPrec =
      guessedExtraPrec * precision * numProds;
  constexpr const int numSqrCoeffs = 4;
  mpfr_t sqrCoeffs[numSqrCoeffs];
  for(int i = 0; i < numSqrCoeffs; i++)
    mpfr_init2(sqrCoeffs[i], sqrTermPrec);
  /* Coefficient 2: -3 */
  mpfr_set_si(sqrCoeffs[3], -3, MPFR_RNDN);
  /* Coefficient 1: 2(c0 + c1 + c2) */
  mpfr_set_d(sqrCoeffs[2], 2 * quad.currentCoeffs[0],
             MPFR_RNDN);
  int err;
  err = mpfr_add_d(sqrCoeffs[2], sqrCoeffs[2],
                   2 * quad.currentCoeffs[1], MPFR_RNDN);
  if(err) return -1;
  err = mpfr_add_d(sqrCoeffs[2], sqrCoeffs[2],
                   2 * quad.currentCoeffs[2], MPFR_RNDN);
  if(err) return -1;
  /* Coefficient 0: c4^2 / 4 + c5^2 / 4 + c7^2 / 4 -
   *                c0 c1 - c0 c2 - c1 c2
   */
  mpfr_set_d(sqrCoeffs[1], quad.currentCoeffs[4] / 2,
             MPFR_RNDN);
  err = mpfr_sqr(sqrCoeffs[1], sqrCoeffs[1], MPFR_RNDN);
  if(err) return -1;
  mpfr_t buf;
  mpfr_init2(buf, sqrTermPrec);

  mpfr_set_d(buf, quad.currentCoeffs[5] / 2, MPFR_RNDN);
  err = mpfr_sqr(buf, buf, MPFR_RNDN);
  if(err) return -1;
  err =
      mpfr_add(sqrCoeffs[1], sqrCoeffs[1], buf, MPFR_RNDN);
  if(err) return -1;

  mpfr_set_d(buf, quad.currentCoeffs[7] / 2, MPFR_RNDN);
  err = mpfr_sqr(buf, buf, MPFR_RNDN);
  if(err) return -1;
  err =
      mpfr_add(sqrCoeffs[1], sqrCoeffs[1], buf, MPFR_RNDN);
  if(err) return -1;

  mpfr_set_d(buf, -quad.currentCoeffs[0], MPFR_RNDN);
  err = mpfr_mul_d(buf, buf, quad.currentCoeffs[1],
                   MPFR_RNDN);
  if(err) return -1;
  err =
      mpfr_add(sqrCoeffs[1], sqrCoeffs[1], buf, MPFR_RNDN);
  if(err) return -1;

  mpfr_set_d(buf, -quad.currentCoeffs[0], MPFR_RNDN);
  err = mpfr_mul_d(buf, buf, quad.currentCoeffs[2],
                   MPFR_RNDN);
  if(err) return -1;
  err =
      mpfr_add(sqrCoeffs[1], sqrCoeffs[1], buf, MPFR_RNDN);
  if(err) return -1;

  mpfr_set_d(buf, -quad.currentCoeffs[1], MPFR_RNDN);
  err = mpfr_mul_d(buf, buf, quad.currentCoeffs[2],
                   MPFR_RNDN);
  if(err) return -1;
  err =
      mpfr_add(sqrCoeffs[1], sqrCoeffs[1], buf, MPFR_RNDN);
  if(err) return -1;
  mpfr_clear(buf);

  /* We have the coefficients of the derivative,
   * now find the roots.
   * This will let us determine the sign of the
   * eigenvalues.
   */
  for(int i = 0; i < numSqrCoeffs; i++)
    mpfr_clear(sqrCoeffs[i]);
  return 0;
}

template <typename fptype>
QuadType classifyRank_0(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_1_1(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_COINCIDENTPLANES;
}

template <typename fptype>
QuadType classifyRank_1_2(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_PARALLELPLANES;
}

template <typename fptype>
QuadType classifyRank_1_3(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_CYLINDER_PAR;
}

template <typename fptype>
QuadType classifyRank_1_4(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_2_1(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_2_2(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  if(eigenSign == 0)
    return QUADT_INTERSECTPLANES_RE;
  else
    return QUADT_INTERSECTPLANES_IM;
}

template <typename fptype>
QuadType classifyRank_2_3(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  if(eigenSign == 0)
    return QUADT_CYLINDER_HYP;
  else
    return QUADT_CYLINDER_ELL;
}

template <typename fptype>
QuadType classifyRank_2_4(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  if(detSign == 1 && eigenSign == 0)
    return QUADT_PARABOLOID_HYP;
  else if(detSign == -1 && eigenSign == 0)
    return QUADT_PARABOLOID_ELL;
  else
    return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_3_1(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_3_2(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_3_3(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  if(eigenSign == 1)
    return QUADT_CONE_IM;
  else
    return QUADT_CONE_RE;
}

template <typename fptype>
QuadType classifyRank_3_4(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  if(eigenSign == 0) {
    if(detSign == 1)
      return QUADT_HYPERBOLOID_ONE;
    else
      return QUADT_HYPERBOLOID_TWO;
  } else {
    if(detSign == -1)
      return QUADT_ELLIPSOID_RE;
    else
      return QUADT_ELLIPSOID_IM;
  }
}

/* Warning: Requires fptype to be recognized by genericfp
 */
template <typename fptype>
QuadType classifyQuadric(
    const Geometry::Quadric<3, fptype> &quad) {
  auto mtxRanks = classifyCalcRank(quad);
  int detSign = classifyCalcDetSign(quad);
  int eigenSign = classifyCalcEigenSign(quad);
  constexpr const int max4Ranks = 4;
  constexpr const int max3Ranks = 3;
  /* The array of function pointers maps the two rank values
   * to functions specific to that rank */
  using classFunc = QuadType (
      *)(int detSign, int eigenSign,
         const Geometry::Quadric<3, fptype> &quad);
  constexpr const classFunc classifiers[max3Ranks +
                                        1][max4Ranks +
                                           1] = {
      {&classifyRank_0<fptype>, &classifyRank_0<fptype>,
       &classifyRank_0<fptype>, &classifyRank_0<fptype>,
       &classifyRank_0<fptype>},
      {&classifyRank_0<fptype>, &classifyRank_1_1<fptype>,
       &classifyRank_1_2<fptype>, &classifyRank_1_3<fptype>,
       &classifyRank_1_4<fptype>},
      {&classifyRank_0<fptype>, &classifyRank_2_1<fptype>,
       &classifyRank_2_2<fptype>, &classifyRank_2_3<fptype>,
       &classifyRank_2_4<fptype>},
      {&classifyRank_0<fptype>, &classifyRank_3_1<fptype>,
       &classifyRank_3_2<fptype>, &classifyRank_3_3<fptype>,
       &classifyRank_3_4<fptype>}};
  QuadType quadclass =
      classifiers[mtxRanks[0]][mtxRanks[1]](
          detSign, eigenSign, quad);
  return quadclass;
}
};

#endif
