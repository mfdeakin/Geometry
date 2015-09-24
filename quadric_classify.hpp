
#ifndef _QUADRIC_CLASSIFY_HPP_
#define _QUADRIC_CLASSIFY_HPP_

#include <mpfr.h>

#include "geometry.hpp"
#include "quadrics.hpp"

namespace QuadricClassify {

enum QuadType {
  QUADT_ERROR = 0,

  QUADT_COINCIDENTPLANES = 1,
  QUADT_INTERSECTPLANES = 2,
  QUADT_INTERSECTPLANES_IM = 3,
  QUADT_INTERSECTPLANES_RE = 4,
  QUADT_PARALLELPLANES = 5,
  QUADT_PARALLELPLANES_IM = 6,
  QUADT_PARALLELPLANES_RE = 7,

  QUADT_ELLIPSOID = 8,
  QUADT_ELLIPSOID_IM = 9,
  QUADT_ELLIPSOID_RE = 10,

  QUADT_CONE = 11,
  QUADT_CONE_IM = 12,
  QUADT_CONE_RE = 13,

  QUADT_CYLINDER_ELL = 14,
  QUADT_CYLINDER_ELL_IM = 15,
  QUADT_CYLINDER_ELL_RE = 16,
  QUADT_CYLINDER_HYP = 17,
  QUADT_CYLINDER_PAR = 18,

  QUADT_PARABOLOID_ELL = 19,
  QUADT_PARABOLOID_HYP = 20,

  QUADT_HYPERBOLOID_ONE = 21,
  QUADT_HYPERBOLOID_TWO = 22,

  QUADT_ERRORINVALID
};

constexpr const char *const QuadTypeNames[] = {
    "Quadratic Classification Error",

    "Coincident Planes", "Intersecting Planes",
    "Imaginary Intersecting Planes",
    "Real Intersecting Planes", "Parallel Planes",
    "Imaginary Parallel Planes", "Real Parallel Planes",

    "Ellipsoid", "Imaginary Ellipsoid", "Real Ellipsoid",

    "Cone", "Imaginary Cone", "Real Cone",

    "Elliptic Cylinder", "Imaginary Elliptic Cylinder",
    "Real Elliptic Cylinder", "Hyperbolic Cylinder",
    "Parabolic Cylinder",

    "Elliptic Paraboloid", "Hyperbolic Paraboloid",

    "Hyperboloid of one sheet",
    "Hyperboloid of two sheets"};

template <typename fptype>
static int classifyCalcDetSign(
    const Geometry::Quadric<3, fptype> &quad) {
  constexpr const int numDetTerms = 17;
  /* Only 4 error causing multiplications occur per term */
  constexpr const int numDetProds = 4;
  /* The determinant is as follows:
   * d = c0 c1 c2 c3 - c2 c3 c4^2 - c1 c3 c5^2 -
   *     c1 c2 c6^2 + 2 c3 c4 c5 c7 - c0 c3 c7^2 +
   *     c6^2 c7^2 + 2 c2 c4 c6 c8 - 2 c5 c6 c7 c8 -
   *     c0 c2 c8^2 + c5^2 c8^2 + 2 c1 c5 c6 c9 -
   *     2 c4 c6 c7 c9 + 2 c4 c5 c8 c9 + 2 c0 c7 c8 c9 -
   *     c0 c1 c9^2 + c4^2 c9^2
   *
   * Note that the quadric coefficients i>3 have an extra
   * factor of 2, so compensate for that in the coefficient
   */
  constexpr const float detCoeffs[] = {
      1,        -1 / 4.0, -1 / 4.0, -1 / 4.0, 1 / 4.0,
      -1 / 4.0, 1 / 16.0, 1 / 4.0,  -1 / 8.0, -1 / 4.0,
      1 / 16.0, 1 / 4.0,  -1 / 8.0, -1 / 8.0, 1 / 4.0,
      -1 / 4.0, 1 / 16.0};
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
      GenericFP::fpconvert<fptype>::precision;
  constexpr const int detTermPrec = precision * numDetProds;
  constexpr const int guessedExtraPrec = 2;
  mpfr_t detTerm, detSum, tmpSum, extra, modAdd;
  mpfr_init2(detTerm, detTermPrec);
  mpfr_init2(detSum, guessedExtraPrec * detTermPrec);
  mpfr_init2(tmpSum, guessedExtraPrec * detTermPrec);
  mpfr_init2(extra, guessedExtraPrec * detTermPrec);
  mpfr_init2(modAdd, guessedExtraPrec * detTermPrec);
  mpfr_set_si(detSum, 0, MPFR_RNDN);
  mpfr_set_si(extra, 0, MPFR_RNDN);
  for(int i = 0; i < numDetTerms; i++) {
    mpfr_set_d(detTerm, detCoeffs[i], MPFR_RNDN);
    for(int j = 0; j < numDetProds; j++) {
      int coeffIdx = detProds[i][j];
      mpfr_mul_d(detTerm, detTerm, quad.coeff(coeffIdx),
                 MPFR_RNDN);
    }
    mpfr_sub(modAdd, detTerm, extra, MPFR_RNDN);
    mpfr_add(tmpSum, modAdd, detSum, MPFR_RNDN);
    mpfr_sub(extra, tmpSum, detSum, MPFR_RNDN);
    mpfr_sub(extra, extra, modAdd, MPFR_RNDN);
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
static void eliminateColumn(mpfr_t elems[mtxDim][mtxDim],
                            int column,
                            bool rowsDone[mtxDim],
                            mpfr_t &coeff, mpfr_t &delta) {
  for(unsigned i = 0; i < mtxDim; i++) {
    int isZero = mpfr_cmp_d(elems[i][column], 0.0);
    if(isZero != 0 && rowsDone[i] == false) {
      for(unsigned j = 0; j < mtxDim; j++) {
        if(j == i) continue;
        isZero = mpfr_cmp_d(elems[j][column], 0.0);
        if(isZero == 0) continue;
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
static std::array<int, 2> classifyCalcRank(
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
   * amount of precision */
  for(unsigned i = 0; i < mtxDim; i++) {
    for(unsigned j = 0; j < mtxDim; j++) {
      mpfr_init2(elems[i][j], precision);
      unsigned coeffNum = mtxVals[i][j];
      fptype factor = (coeffNum < 4) ? 1.0 : 0.5;
      mpfr_set_d(elems[i][j], quad.coeff(coeffNum) * factor,
                 MPFR_RNDN);
    }
  }
  mpfr_t coeff, delta;
  mpfr_inits2(precision, coeff, delta, (mpfr_ptr)NULL);
  bool rowsDone[mtxDim];
  memset(rowsDone, false, sizeof(rowsDone));
  for(int numRowsDone = 0; numRowsDone < mtxDim;
      numRowsDone++)
    eliminateColumn<mtxDim>(elems, numRowsDone, rowsDone,
                            coeff, delta);
  mpfr_clears(coeff, delta, (mpfr_ptr)NULL);
  std::array<int, 2> ranks;
  ranks[0] = 0;
  ranks[1] = 0;
  bool rowChecked;
  for(int i = 0; i < mtxDim; i++) {
    rowChecked = false;
    for(int j = 0; j < mtxDim; j++) {
      int isZero = mpfr_cmp_d(elems[i][j], 0.0);
      if(isZero != 0 && rowChecked == false) {
        ranks[1]++;
        if(i < (mtxDim - 1) && j < (mtxDim - 1)) ranks[0]++;
        rowChecked = true;
      }
      mpfr_clear(elems[i][j]);
    }
  }
  return ranks;
}

template <typename fptype>
static mpfr_ptr constructCubicCoeffs(
    const Geometry::Quadric<3, fptype> &quad,
    unsigned precision) {
  /* The cubic is of the following form:
   * -x^3 + (c0 + c1 + c2)x^2 +
   * (c4^2/4 + c5^2/4 + c7^2/4 - c0 c1 - c0 c2 - c1 c2)x +
   * c0 c1 c2 - c0 c7^2/4 - c1 c5^2/4 -
   * c2 c4^2/4 + c4 c5 c7/4
   */
  constexpr const unsigned numCoeffs = 4;
  mpfr_ptr coeffs = static_cast<mpfr_ptr>(
      malloc(sizeof(mpfr_t[numCoeffs])));
  for(unsigned i = 0; i < numCoeffs; i++)
    mpfr_init2(&coeffs[i], precision);
  /* Coefficient 3: -1 */
  mpfr_set_d(&coeffs[3], -1.0, MPFR_RNDN);
  /* Coefficient 2: c0 + c1 + c2 */
  mpfr_set_d(&coeffs[2], quad.coeff(0), MPFR_RNDN);
  mpfr_add_d(&coeffs[2], &coeffs[2], quad.coeff(1),
             MPFR_RNDN);
  mpfr_add_d(&coeffs[2], &coeffs[2], quad.coeff(2),
             MPFR_RNDN);
  /* Coefficient 1: c4^2/4 + c5^2/4 + c7^2/4 -
   *                c0 c1 - c0 c2 - c1 c2
   */
  mpfr_set_d(&coeffs[1], quad.coeff(4) / 2.0, MPFR_RNDN);
  mpfr_sqr(&coeffs[1], &coeffs[1], MPFR_RNDN);
  mpfr_t buf;
  mpfr_init2(buf, precision);

  mpfr_set_d(buf, quad.coeff(5) / 2.0, MPFR_RNDN);
  mpfr_sqr(buf, buf, MPFR_RNDN);
  mpfr_add(&coeffs[1], &coeffs[1], buf, MPFR_RNDN);

  mpfr_set_d(buf, quad.coeff(7) / 2.0, MPFR_RNDN);
  mpfr_sqr(buf, buf, MPFR_RNDN);
  mpfr_add(&coeffs[1], &coeffs[1], buf, MPFR_RNDN);

  mpfr_set_d(buf, -quad.coeff(0), MPFR_RNDN);
  mpfr_mul_d(buf, buf, quad.coeff(1), MPFR_RNDN);
  mpfr_add(&coeffs[1], &coeffs[1], buf, MPFR_RNDN);

  mpfr_set_d(buf, -quad.coeff(0), MPFR_RNDN);
  mpfr_mul_d(buf, buf, quad.coeff(2), MPFR_RNDN);
  mpfr_add(&coeffs[1], &coeffs[1], buf, MPFR_RNDN);

  mpfr_set_d(buf, -quad.coeff(1), MPFR_RNDN);
  mpfr_mul_d(buf, buf, quad.coeff(2), MPFR_RNDN);
  mpfr_add(&coeffs[1], &coeffs[1], buf, MPFR_RNDN);
  /* Coefficient 0: c0 c1 c2 + c4 c5 c7/4 - c0 c7^2/4 -
   *                c1 c5^2/4 - c2 c4^2/4
   */
  /* c0 c1 c2 */
  mpfr_set_d(&coeffs[0], quad.coeff(0), MPFR_RNDN);
  mpfr_mul_d(&coeffs[0], &coeffs[0], quad.coeff(1),
             MPFR_RNDN);
  mpfr_mul_d(&coeffs[0], &coeffs[0], quad.coeff(2),
             MPFR_RNDN);
  /* c4 c5 c7/4 */
  mpfr_set_d(buf, quad.coeff(4) / 4.0, MPFR_RNDN);
  mpfr_mul_d(buf, buf, quad.coeff(5), MPFR_RNDN);
  mpfr_mul_d(buf, buf, quad.coeff(7), MPFR_RNDN);
  mpfr_add(&coeffs[0], &coeffs[0], buf, MPFR_RNDN);
  /* -c0 c7^2/4 */
  mpfr_set_d(buf, quad.coeff(7), MPFR_RNDN);
  mpfr_sqr(buf, buf, MPFR_RNDN);
  mpfr_mul_d(buf, buf, -quad.coeff(0) / 4.0, MPFR_RNDN);
  mpfr_add(&coeffs[0], &coeffs[0], buf, MPFR_RNDN);
  /* -c1 c5^2/4 */
  mpfr_set_d(buf, quad.coeff(5), MPFR_RNDN);
  mpfr_sqr(buf, buf, MPFR_RNDN);
  mpfr_mul_d(buf, buf, -quad.coeff(1) / 4.0, MPFR_RNDN);
  mpfr_add(&coeffs[0], &coeffs[0], buf, MPFR_RNDN);
  /* -c2 c4^2/4 */
  mpfr_set_d(buf, quad.coeff(4), MPFR_RNDN);
  mpfr_sqr(buf, buf, MPFR_RNDN);
  mpfr_mul_d(buf, buf, -quad.coeff(2) / 4.0, MPFR_RNDN);
  mpfr_add(&coeffs[0], &coeffs[0], buf, MPFR_RNDN);
  mpfr_clear(buf);
  return coeffs;
}

static inline mpfr_ptr calcInflections(mpfr_ptr cubic,
                                       unsigned precision) {
  /* Compute the derivative and it's roots.
   * The input is as follows:
   * c3 x^3 + c2 x^2 + c1 x + c0
   * So the derivative is as follows:
   * 3 c3 x^2 + 2 c2 x + c1
   */
  /* Start with the discriminant
   * The discriminant is as follows:
   * (2 c2)^2 - 4(3 c3 c1)
   */
  mpfr_t disc;
  mpfr_init2(disc, precision);
  mpfr_mul(disc, &cubic[1], &cubic[3], MPFR_RNDN);
  mpfr_mul_d(disc, disc, -12.0, MPFR_RNDN);
  mpfr_t buf;
  mpfr_init2(buf, precision);
  mpfr_sqr(buf, &cubic[2], MPFR_RNDN);
  mpfr_mul_d(buf, buf, 4.0, MPFR_RNDN);
  mpfr_add(disc, disc, buf, MPFR_RNDN);
  /* Now compute the roots.
   * The roots are as follows:
   * -sign(c2) (|2 c2|+sqrt(disc)) / (3 c3)
   * -sign(c2) c1 / (|2 c2|+sqrt(disc))
   */
  mpfr_sqrt(disc, disc, MPFR_RNDN);
  constexpr const int numRoots = 2;
  mpfr_ptr roots = static_cast<mpfr_ptr>(
      malloc(sizeof(mpfr_t[numRoots])));
  for(int i = 0; i < numRoots; i++)
    mpfr_init2(&roots[i], precision);

  int c2Sign = mpfr_signbit(&cubic[2]);
  mpfr_abs(buf, &cubic[2], MPFR_RNDN);
  mpfr_mul_d(buf, buf, 2.0, MPFR_RNDN);
  mpfr_add(buf, buf, disc, MPFR_RNDN);
  if(c2Sign == 0) {
    mpfr_mul_d(buf, buf, -0.5, MPFR_RNDN);
  } else {
    mpfr_mul_d(buf, buf, 0.5, MPFR_RNDN);
  }
  mpfr_div(&roots[0], buf, &cubic[3], MPFR_RNDN);
  mpfr_div_d(&roots[0], &roots[0], 3.0, MPFR_RNDN);
  mpfr_div(&roots[1], &cubic[1], buf, MPFR_RNDN);
  mpfr_clear(disc);
  mpfr_clear(buf);
  return roots;
}

template <typename fptype>
static int classifyCalcEigenSign(
    const Geometry::Quadric<3, fptype> &quad) {
  /* Basically I'm computing the roots of the derivative
   * of the characteristic cubic of the matrix.
   * This is easier, and with the constant term,
   * sufficient to determine the signs of the eigenvalues.
   * Note that since this is a real symmetric matrix,
   * three real eigenvalues are guaranteed.
   *
   * The cubic is of the following form:
   * -x^3 + (c0 + c1 + c2)x^2 +
   * (c4^2/4 + c5^2/4 + c7^2/4 - c0 c1 - c0 c2 - c1 c2)x +
   * c0 c1 c2 - c0 c7^2/4 - c1 c5^2/4 -
   * c2 c4^2/4 + c4 c5 c7/4
   *
   * The derivative is of the following form:
   * -3 x^2 + 2(c0 + c1 + c2)x +
   * ((c4^2 + c5^2 + c7^2)/4 - c0 c1 - c0 c2 - c1 c2)
   */
  constexpr const int precision =
      GenericFP::fpconvert<fptype>::precision;
  constexpr const int numCubeProds = 3;
  constexpr const int guessedExtraPrec = 2;
  constexpr const int cubeTermPrec =
      guessedExtraPrec * precision * numCubeProds;
  constexpr const int numCubeCoeffs = 4;
  mpfr_ptr cubicCoeffs =
      constructCubicCoeffs(quad, cubeTermPrec);
  /* We have the coefficients of the derivative,
   * now find the roots.
   * This will let us determine the sign of the
   * eigenvalues.
   */
  constexpr const int numRoots = 2;
  mpfr_ptr roots =
      calcInflections(cubicCoeffs, cubeTermPrec);
  const int zeroVal = mpfr_cmp_d(&cubicCoeffs[0], 0.0);
  const int drootSigns[] = {mpfr_cmp_d(&roots[0], 0.0),
                            mpfr_cmp_d(&roots[1], 0.0)};
  for(int i = 0; i < numCubeCoeffs; i++)
    mpfr_clear(&cubicCoeffs[i]);
  free(cubicCoeffs);
  for(int i = 0; i < numRoots; i++) mpfr_clear(&roots[i]);
  free(roots);

  bool drootsPlus =
      (drootSigns[0] >= 0) && (drootSigns[1] >= 0);
  bool drootsMinus =
      (drootSigns[0] <= 0) && (drootSigns[1] <= 0);
  if((drootsPlus && zeroVal >= 0) ||
     (drootsMinus && zeroVal <= 0))
    return 1;
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
  else if(detSign == -1 && eigenSign == 1)
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
    if(detSign < 0)
      return QUADT_ELLIPSOID_RE;
    else
      return QUADT_ELLIPSOID_IM;
  }
}

/* Warning: Requires fptype to be recognized by genericfp
 */
template <typename fptype>
static QuadType classifyQuadric(
    const Geometry::Quadric<3, fptype> &quad) {
  auto mtxRanks = classifyCalcRank(quad);
  int detSign = classifyCalcDetSign(quad);
  int eigenSign = classifyCalcEigenSign(quad);
  constexpr const int max4Ranks = 4;
  constexpr const int max3Ranks = 3;
  printf(
      "\nSmall Matrix Rank: %d\n"
      "Large Matrix Rank: %d\n"
      "Determinant Sign: %d\n"
      "Eigenvalue Sign: %d\n",
      mtxRanks[0], mtxRanks[1], detSign, eigenSign);
  assert(mtxRanks[0] <= max3Ranks);
  assert(mtxRanks[1] <= max4Ranks);
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
}

#endif
