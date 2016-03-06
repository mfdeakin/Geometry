
#ifndef _QUADRIC_CLASSIFY_HPP_
#define _QUADRIC_CLASSIFY_HPP_

#include "mpreal.hpp"

#include "geometry.hpp"
#include "quadric.hpp"

#include "array.hpp"

namespace QuadricClassify {

enum QuadType {
  QUADT_ERROR = 0,
  QUADT_DEGENERATE = 1,

  QUADT_COINCIDENTPLANES = 2,
  QUADT_INTERSECTPLANES = 3,
  QUADT_INTERSECTPLANES_IM = 4,
  QUADT_INTERSECTPLANES_RE = 5,
  QUADT_PARALLELPLANES = 6,
  QUADT_PARALLELPLANES_IM = 7,
  QUADT_PARALLELPLANES_RE = 8,

  QUADT_ELLIPSOID = 9,
  QUADT_ELLIPSOID_IM = 10,
  QUADT_ELLIPSOID_RE = 11,

  QUADT_CONE = 12,
  QUADT_CONE_IM = 13,
  QUADT_CONE_RE = 14,

  QUADT_CYLINDER_ELL = 15,
  QUADT_CYLINDER_ELL_IM = 16,
  QUADT_CYLINDER_ELL_RE = 17,
  QUADT_CYLINDER_HYP = 18,
  QUADT_CYLINDER_PAR = 19,

  QUADT_PARABOLOID_ELL = 20,
  QUADT_PARABOLOID_HYP = 21,

  QUADT_HYPERBOLOID_ONE = 22,
  QUADT_HYPERBOLOID_TWO = 23,

  QUADT_PLANE = 24,

  QUADT_ERRORINVALID
};

constexpr const char *const QuadTypeNames[] = {
    "Quadratic Classification Error",
    "Degenerate",

    "Coincident Planes",
    "Intersecting Planes",
    "Imaginary Intersecting Planes",
    "Real Intersecting Planes",
    "Parallel Planes",
    "Imaginary Parallel Planes",
    "Real Parallel Planes",

    "Ellipsoid",
    "Imaginary Ellipsoid",
    "Real Ellipsoid",

    "Cone",
    "Imaginary Cone",
    "Real Cone",

    "Elliptic Cylinder",
    "Imaginary Elliptic Cylinder",
    "Real Elliptic Cylinder",
    "Hyperbolic Cylinder",
    "Parabolic Cylinder",

    "Elliptic Paraboloid",
    "Hyperbolic Paraboloid",

    "Hyperboloid of one sheet",
    "Hyperboloid of two sheets",

    "Plane"};

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
      {0, 1, 2, 3}, {2, 3, 4, 4}, {1, 3, 5, 5},
      {1, 2, 6, 6}, {3, 4, 5, 7}, {0, 3, 7, 7},
      {6, 6, 7, 7}, {2, 4, 6, 8}, {5, 6, 7, 8},
      {0, 2, 8, 8}, {5, 5, 8, 8}, {1, 5, 6, 9},
      {4, 6, 7, 9}, {4, 5, 8, 9}, {0, 7, 8, 9},
      {0, 1, 9, 9}, {4, 4, 9, 9}};
  const int precision =
      MathFuncs::MathFuncs<fptype>::getPrec(fptype(0.0));
  const int detTermPrec = precision * numDetProds;
  constexpr const int guessedExtraPrec = 2;
  const int prevPrec = mpfr::mpreal::get_default_prec();
  mpfr::mpreal::set_default_prec(guessedExtraPrec *
                                 detTermPrec);
  /* Use Kahan Summation to sum up the terms */
  mpfr::mpreal detTerm, detSum(0.0), tmpSum, extra(0.0),
      modAdd;
  for(int i = 0; i < numDetTerms; i++) {
    detTerm = detCoeffs[i];
    for(int j = 0; j < numDetProds; j++) {
      int coeffIdx = detProds[i][j];
      detTerm *= quad.coeff(coeffIdx);
    }
    modAdd = detTerm - extra;
    tmpSum = modAdd + detSum;
    extra = tmpSum - detSum;
    extra -= modAdd;
    detSum = tmpSum;
  }
  int cmpZero = 0;
  if(detSum < 0)
    cmpZero = -1;
  else if(detSum > 0)
    cmpZero = 1;
  return cmpZero;
}

template <unsigned mtxDim>
static void eliminateColumn(
    mpfr::mpreal elems[mtxDim][mtxDim], int column,
    bool rowsDone[mtxDim], mpfr::mpreal &coeff,
    mpfr::mpreal &delta) {
  for(unsigned i = 0; i < mtxDim; i++) {
    if(elems[i][column] != 0.0 && rowsDone[i] == false) {
      for(unsigned j = 0; j < mtxDim; j++) {
        if(j == i) continue;
        if(elems[j][column] == 0.0) continue;
        coeff = elems[j][column] / elems[i][column];
        for(unsigned k = column + 1; k < mtxDim; k++) {
          delta = elems[i][k] * coeff;
          elems[j][k] -= delta;
        }
        elems[j][column] = 0.0;
      }
      rowsDone[i] = true;
      break;
    }
  }
}

template <typename fptype>
static Array<int, 2> classifyCalcRank(
    const Geometry::Quadric<3, fptype> &quad) {
  constexpr const int mtxDim = 4;
  constexpr const int precision = 512;
  constexpr const unsigned mtxVals[mtxDim][mtxDim] = {
      {0, 4, 5, 6},
      {4, 1, 7, 8},
      {5, 7, 2, 9},
      {6, 8, 9, 3}};
  const int prevPrec = mpfr::mpreal::get_default_prec();
  mpfr::mpreal::set_default_prec(precision);
  mpfr::mpreal elems[mtxDim][mtxDim];
  /* Just use Gaussian Elimination with a ridiculous
   * amount of precision */
  for(unsigned i = 0; i < mtxDim; i++) {
    for(unsigned j = 0; j < mtxDim; j++) {
      unsigned coeffNum = mtxVals[i][j];
      fptype factor = (coeffNum < 4) ? 1.0 : 0.5;
      elems[i][j] = quad.coeff(coeffNum);
      elems[i][j] *= factor;
    }
  }
  mpfr::mpreal coeff, delta;
  bool rowsDone[mtxDim];
  memset(rowsDone, false, sizeof(rowsDone));
  for(int numRowsDone = 0; numRowsDone < mtxDim;
      numRowsDone++)
    eliminateColumn<mtxDim>(elems, numRowsDone, rowsDone,
                            coeff, delta);
  Array<int, 2> ranks;
  ranks[0] = 0;
  ranks[1] = 0;
  bool rowChecked;
  for(int i = 0; i < mtxDim; i++) {
    rowChecked = false;
    for(int j = 0; j < mtxDim; j++) {
      if(elems[i][j] != 0 && rowChecked == false) {
        ranks[1]++;
        if(i < (mtxDim - 1) && j < (mtxDim - 1)) ranks[0]++;
        rowChecked = true;
      }
      elems[i][j] = NAN;
    }
  }
  mpfr::mpreal::set_default_prec(prevPrec);
  return ranks;
}

template <typename fptype>
static mpfr::mpreal *constructCubicCoeffs(
    const Geometry::Quadric<3, fptype> &quad,
    unsigned precision) {
  const int prevPrec = mpfr::mpreal::get_default_prec();
  mpfr::mpreal::set_default_prec(precision);
  /* The cubic is of the following form:
   * -x^3 + (c0 + c1 + c2)x^2 +
   * (c4^2/4 + c5^2/4 + c7^2/4 - c0 c1 - c0 c2 - c1 c2)x +
   * c0 c1 c2 - c0 c7^2/4 - c1 c5^2/4 -
   * c2 c4^2/4 + c4 c5 c7/4
   */
  constexpr const unsigned numCoeffs = 4;
  mpfr::mpreal *coeffs = new mpfr::mpreal[numCoeffs];
  /* Coefficient 3: -1 */
  coeffs[3] = -1.0;
  /* Coefficient 2: c0 + c1 + c2 */
  coeffs[2] = quad.coeff(0) + quad.coeff(1) + quad.coeff(2);
  /* Coefficient 1: c4^2/4 + c5^2/4 + c7^2/4 -
   *                c0 c1 - c0 c2 - c1 c2
   */
  coeffs[1] = quad.coeff(4) * quad.coeff(4);
  coeffs[1] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      quad.coeff(5), quad.coeff(5), coeffs[1]);
  coeffs[1] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      quad.coeff(7), quad.coeff(7), coeffs[1]);
  coeffs[1] /= 4;
  coeffs[1] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      -quad.coeff(0), quad.coeff(1), coeffs[1]);
  coeffs[1] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      -quad.coeff(0), quad.coeff(2), coeffs[1]);
  coeffs[1] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      -quad.coeff(1), quad.coeff(2), coeffs[1]);
  /* Coefficient 0: c0 c1 c2 + c4 c5 c7/4 - c0 c7^2/4 -
   *                c1 c5^2/4 - c2 c4^2/4
   */
  /* c0 c1 c2 */
  coeffs[0] = mpfr::mpreal(quad.coeff(0)) * quad.coeff(1) *
              quad.coeff(2);
  /* c4 c5 c7/4 */
  coeffs[0] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      quad.coeff(4) / 4.0,
      mpfr::mpreal(quad.coeff(5)) * quad.coeff(7),
      coeffs[0]);
  /* -c0 c7^2/4 */
  coeffs[0] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      -quad.coeff(0) / 4.0,
      mpfr::mpreal(quad.coeff(7)) * quad.coeff(7),
      coeffs[0]);
  /* -c1 c5^2/4 */
  coeffs[0] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      -quad.coeff(1) / 4.0,
      mpfr::mpreal(quad.coeff(5)) * quad.coeff(5),
      coeffs[0]);
  /* -c2 c4^2/4 */
  coeffs[0] = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      -quad.coeff(2) / 4.0,
      mpfr::mpreal(quad.coeff(4)) * quad.coeff(4),
      coeffs[0]);
  mpfr::mpreal::set_default_prec(prevPrec);
  return coeffs;
}

static inline mpfr::mpreal *calcInflections(
    mpfr::mpreal *cubic, unsigned precision) {
  const int prevPrec = mpfr::mpreal::get_default_prec();
  mpfr::mpreal::set_default_prec(precision);
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
  mpfr::mpreal disc = -12.0 * cubic[1] * cubic[3];
  disc = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
      4.0 * cubic[2], cubic[2], disc);
  mpfr::mpreal width =
      MathFuncs::MathFuncs<mpfr::mpreal>::sqrt(disc);
  /* Now compute the roots.
   * The roots are as follows:
   * -sign(c2) (|2 c2|+sqrt(disc)) / (3 c3)
   * -sign(c2) c1 / (|2 c2|+sqrt(disc))
   */
  constexpr const int numRoots = 2;
  mpfr::mpreal *roots = new mpfr::mpreal[numRoots];

  bool c2Sign =
      MathFuncs::MathFuncs<mpfr::mpreal>::signbit(cubic[2]);
  roots[0] = MathFuncs::MathFuncs<mpfr::mpreal>::abs(
      2.0 * cubic[2]);
  roots[0] += width;
  if(c2Sign == 0)
    roots[0] /= -2.0;
  else
    roots[0] /= 2.0;
  roots[1] = cubic[1] / roots[0];
  roots[0] /= (3.0 * cubic[3]);

  mpfr::mpreal::set_default_prec(prevPrec);
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
  const int precision =
      MathFuncs::MathFuncs<fptype>::getPrec(fptype(0.0));
  constexpr const int numCubeProds = 3;
  constexpr const int guessedExtraPrec = 2;
  const int cubeTermPrec =
      guessedExtraPrec * precision * numCubeProds;
  constexpr const int numCubeCoeffs = 4;
  mpfr::mpreal *cubicCoeffs =
      constructCubicCoeffs(quad, cubeTermPrec);
  /* We have the coefficients of the derivative,
   * now find the roots.
   * This will let us determine the sign of the
   * eigenvalues.
   */
  constexpr const int numRoots = 2;
  mpfr::mpreal *roots =
      calcInflections(cubicCoeffs, cubeTermPrec);
  const int zeroVal = (cubicCoeffs[0] == 0.0);
  const float drootSigns[] = {
      MathFuncs::MathFuncs<float>::copysign(
          1.0, float(roots[0])),
      MathFuncs::MathFuncs<float>::copysign(
          1.0, float(roots[1]))};
  delete[] cubicCoeffs;
  delete[] roots;

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
QuadType classifyEllCylinder(
    const Geometry::Quadric<3, fptype> &quad) {
  /* This is an elliptic cylinder
   * We have it in the form c_{xx} x^2 + c_{xy} x y + ... +
   *                          c_z z + c_0
   * We want it in the form c_{xx} x^2 + c_{yy} y^2 + ... +
   *                          c_z z + c_0
   * Consider c_{xx} x^2 + c_{xy} x y + c_{yy} y^2
   * We can eliminate c_{xy} with a rotation in the xy plane
   * c_xx (x cos(theta) - y sin(theta))^2 +
   *   c_xy (x cos(theta) - y sin(theta)) *
   *     (x sin(theta) + y cos(theta)) +
   *   c_yy (x sin(theta) + y cos(theta))^2
   * This has two solutions:
   * x cos(theta) = y sin(theta)
   * x sin(theta) = -y cos(theta)
   */
  return QUADT_CYLINDER_ELL;
}

template <typename fptype>
QuadType classifyRank_Error(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_ERROR;
}

template <typename fptype>
QuadType classifyRank_0_0(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_DEGENERATE;
}

template <typename fptype>
QuadType classifyRank_0_2(
    int detSign, int eigenSign,
    const Geometry::Quadric<3, fptype> &quad) {
  return QUADT_PLANE;
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
    return classifyEllCylinder(quad);
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

static bool isImaginary(QuadType qt) {
  switch(qt) {
    case QUADT_INTERSECTPLANES_IM:
    case QUADT_PARALLELPLANES_IM:
    case QUADT_ELLIPSOID_IM:
    case QUADT_CONE_IM:
    case QUADT_CYLINDER_ELL_IM:
      return true;
    default:
      return false;
  }
}

static bool isReal(QuadType qt) {
  /* Assume unknown quadric types are real until
   * better tests are implemented */
  switch(qt) {
    case QUADT_COINCIDENTPLANES:
    case QUADT_INTERSECTPLANES:
    case QUADT_INTERSECTPLANES_RE:
    case QUADT_PARALLELPLANES:
    case QUADT_PARALLELPLANES_RE:
    case QUADT_ELLIPSOID:
    case QUADT_ELLIPSOID_RE:
    case QUADT_CONE:
    case QUADT_CONE_RE:
    case QUADT_CYLINDER_ELL:
    case QUADT_CYLINDER_ELL_RE:
    case QUADT_CYLINDER_HYP:
    case QUADT_CYLINDER_PAR:
    case QUADT_PARABOLOID_ELL:
    case QUADT_PARABOLOID_HYP:
    case QUADT_HYPERBOLOID_ONE:
    case QUADT_HYPERBOLOID_TWO:
    case QUADT_PLANE:
      return true;
    default:
      return false;
  }
}

template <typename fptype>
static QuadType classifyQuadric(
    const Geometry::Quadric<3, fptype> &quad) {
  auto mtxRanks = classifyCalcRank(quad);
  int detSign = classifyCalcDetSign(quad);
  int eigenSign = classifyCalcEigenSign(quad);
  constexpr const int max4Ranks = 4;
  constexpr const int max3Ranks = 3;
  assert(mtxRanks[0] <= max3Ranks);
  assert(mtxRanks[1] <= max4Ranks);
  /* The array of function pointers maps the two rank values
   * to functions specific to that rank */
  using classFunc = QuadType (
          *)(int detSign, int eigenSign,
             const Geometry::Quadric<3, fptype> &quad);
  constexpr const classFunc
      classifiers[max3Ranks + 1][max4Ranks + 1] = {
          {&classifyRank_0_0<fptype>,
           &classifyRank_Error<fptype>,
           &classifyRank_0_2<fptype>,
           &classifyRank_Error<fptype>,
           &classifyRank_Error<fptype>},
          {&classifyRank_Error<fptype>,
           &classifyRank_1_1<fptype>,
           &classifyRank_1_2<fptype>,
           &classifyRank_1_3<fptype>,
           &classifyRank_Error<fptype>},
          {&classifyRank_Error<fptype>,
           &classifyRank_Error<fptype>,
           &classifyRank_2_2<fptype>,
           &classifyRank_2_3<fptype>,
           &classifyRank_2_4<fptype>},
          {&classifyRank_Error<fptype>,
           &classifyRank_Error<fptype>,
           &classifyRank_Error<fptype>,
           &classifyRank_3_3<fptype>,
           &classifyRank_3_4<fptype>}};
  QuadType quadclass =
      classifiers[mtxRanks[0]][mtxRanks[1]](
          detSign, eigenSign, quad);
  return quadclass;
}
}

#endif
