
#ifndef _QUADRIC_HPP_
#define _QUADRIC_HPP_

#include <stdlib.h>
#include <assert.h>

#include <ostream>
#include <istream>
#include <iostream>
#include <memory>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "line.hpp"

#include "quadric_classify.hpp"

#include "array.hpp"

#include "polynomial.hpp"

#include "accurate_math.hpp"

#include "mathfuncs.hpp"

#include "mpreal.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Quadric : public Solid<dim, fptype> {
 public:
  struct QuadricData {
    typename Origin<dim, fptype>::OriginData origin;
    Array<fptype, Quadric<dim, fptype>::numCoeffs> coeffs;
  };

  CUDA_CALLABLE Quadric() : Solid<dim, fptype>() {
    for(int i = 0; i < numCoeffs; i++) setCoeff(i, NAN);
  }

  CUDA_CALLABLE Quadric(const Origin<dim, fptype> &origin)
      : Solid<dim, fptype>(origin) {
    for(int i = 0; i < numCoeffs; i++) setCoeff(i, NAN);
  }

  template <typename srctype>
  CUDA_CALLABLE Quadric(const Quadric<dim, srctype> &src)
      : Solid<dim, fptype>(src.origin) {
    for(int i = 0; i < numCoeffs; i++)
      setCoeff(i, src.coeff(i));
  }

  CUDA_CALLABLE Quadric(const QuadricData &src)
      : Solid<dim, fptype>(Origin<dim, fptype>(src.origin)),
        coeffs(src.coeffs) {}

  CUDA_CALLABLE virtual ~Quadric() {}

  /* d=0 corresponds to the first dimension,
   * d=1 the second, ...
   * d=dim the constant coefficient
   */
  CUDA_CALLABLE fptype coeff(int d1, int d2) const {
    int coeffNum = getCoeffPos(d1, d2);
    return coeffs[coeffNum];
  }

  CUDA_CALLABLE fptype setCoeff(int d1, int d2,
                                fptype val) {
    int coeffNum = getCoeffPos(d1, d2);
    return setCoeff(coeffNum, val);
  }

  /* Access the coefficients as a 1D array.
   * pos=0 to dim-1 are for the square terms
   * pos=dim is the constant term
   * pos=dim+1 to 2dim-1 is the first dim prod terms
   * pos=2dim is the first dim linear term
   * ...
   * pos=(dim+2)*(dim+1)/2-1 is the last dim linear term
   */
  CUDA_CALLABLE fptype coeff(int pos) const {
    assert(pos >= 0);
    assert(pos < numCoeffs);
    return coeffs[pos];
  }

  CUDA_CALLABLE fptype setCoeff(int pos, fptype val) {
    assert(pos >= 0);
    assert(pos < numCoeffs);
    coeffs[pos] = val;
    return coeffs[pos];
  }

  /* For finite quadric surfaces,
   * compute the point which has the minimum average
   * distance to the surface.
   * This will likely reduce the error of computations in
   * this space.
   * Currently only implemented for the real ellipsoidal
   * surface.
   */
  CUDA_CALLABLE void optimizeOrigin() {
    auto qt = QuadricClassify::classifyQuadric(*this);
    switch(qt) {
      case QuadricClassify::QUADT_ELLIPSOID:
      case QuadricClassify::QUADT_ELLIPSOID_IM:
      case QuadricClassify::QUADT_ELLIPSOID_RE:
        /* Calculate the center of the ellipsoid and use it
         * as the origin.
         * This optimizes the average distance and the
         * minimax distance, so should be optimal for
         * reducing error. */
        mpfr::mpreal rad(coeff(dim, dim));
        Array<fptype, dim> center;
        for(int i = 0; i < dim; i++) {
          fptype mult = coeff(i, i);
          center[i] = -coeff(i, dim) / 2.0 / mult;
          coeff(i, dim) = 0.0;
          mpfr::mpreal tmp(center[i]);
          tmp =
              MathFuncs::MathFuncs<mpfr::mpreal>::sqr(tmp);
          rad = MathFuncs::MathFuncs<mpfr::mpreal>::fma(
              -tmp, mult, rad);
        }
        coeff(dim, dim) = static_cast<fptype>(rad);
        this->origin(center);
        break;
    }
  }

  /* Returns a signed value which determines the side of the
   * surface that the point is on.
   * This is done by simply evaluating
   * c00 p0 p0 + c01 p0 p1 + ... + c11 p1 p1 + c12 p1 p2 +
   * ...
   */
  template <typename outtype>
  CUDA_CALLABLE outtype evaluatePoint(
      const Point<dim, fptype> &pt,
      fptype absPrecision = defAbsPrecision) const {
    /* Use Kahan summation to evaluate this more correctly
     */
    outtype ret = 0.0;
    outtype extra = 0.0;
    Point<dim, fptype> copy(pt);
    copy.shiftOrigin(this->origin);
    Vector<dim, fptype> quadOffset = copy.getOffset();
    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        int coeffNum = getCoeffPos(i, j);
        outtype product = coeff(coeffNum) * quadOffset(i) *
                              quadOffset(j) -
                          extra;
        outtype tmp = ret + product;
        extra = (tmp - ret) - product;
        ret = tmp;
      }
    }
    return ret;
  }

  /* Calculates the quadratic polynomial with roots
   * at the points at which the line's parameter
   * will intersect with the quadric surface */
  CUDA_CALLABLE Polynomial<2, fptype> calcLineDistPoly(
      const Line<dim, fptype> &line,
      fptype absPrecision = defAbsPrecision) const {
    /* Notation: The line is represented as $d t + p$
     * The component of the direction or point in the
     * i'th direction is given as d[i] and p[i],
     * respectively */
    auto lDir = line.getDirection();
    auto lInt = line.getIntercept().getOffset();
    /* The square coefficient is only composed of terms
     * from the direction terms.
     * This term is easily seen to be the sum of
     * c[i, j] d[i] d[j].
     */
    fptype sqCoeff = 0.0;
    /* The linear coefficient is composed of a combination
     * of direction and intercept terms.
     *
     */
    fptype linCoeff = 0.0;
    /* The constant term is only composed of terms
     * from the intercept coefficients.
     * This term is easily seen to be the sum of
     * 2 c[i, dim] p[i] p[j]
     * and c[dim, dim]
     */
    fptype constant = 0.0;
    constant += coeff(dim, dim);
    for(int i = 0; i < dim; i++) {
      sqCoeff =
          fptype(fptype(coeff(i, i) * lDir(i)) * lDir(i)) +
          sqCoeff;
      linCoeff =
          fptype(MathFuncs::MathFuncs<fptype>::fastMult2Pow(
                     coeff(i, i), 1) *
                 fptype(lDir(i) * lInt(i))) +
          fptype(fptype(coeff(i, dim) * lDir(i)) +
                 linCoeff);
      constant =
          fptype(fptype(coeff(i, dim) * lInt(i)) +
                 fptype(fptype(coeff(i, i) * lInt(i)) *
                        lInt(i))) +
          constant;
      for(int j = i + 1; j < dim; j++) {
        sqCoeff = fptype(fptype(coeff(i, j) * lDir(i)) *
                         lDir(j)) +
                  sqCoeff;
        linCoeff =
            fptype(fptype(fptype(coeff(i, j) * lDir(i)) *
                          lInt(j)) +
                   fptype(fptype(coeff(i, j) * lDir(j)) *
                          lInt(i))) +
            linCoeff;
        constant = fptype(fptype(coeff(i, j) * lInt(i)) *
                          lInt(j)) +
                   constant;
      }
    }
    Polynomial<2, fptype> poly;
    poly.set(0, constant);
    poly.set(1, linCoeff);
    poly.set(2, sqCoeff);
    return poly;
  }

  CUDA_CALLABLE Array<fptype, 2> calcLineDistToIntersect(
      Line<dim, fptype> line,
      fptype absPrecision = defAbsPrecision) const {
    fptype shiftDist =
        line.argPointMinDist(Point<dim, fptype>(
            this->origin, Vector<dim, fptype>()));
    line.shiftOrigin(this->origin);
    Polynomial<2, fptype> poly =
        calcLineDistPoly(line, absPrecision);
    auto roots = poly.calcRoots();
    for(int i = 0; i < 2; i++) roots[i] += shiftDist;
    return roots;
  }

  CUDA_CALLABLE Array<Point<dim, fptype>, 2>
  calcLineIntersect(
      const Line<dim, fptype> &line,
      fptype absPrecision = defAbsPrecision) const {
    Array<fptype, 2> roots =
        calcLineDistToIntersect(line, absPrecision);
    return Array<Point<dim, fptype>, 2>(
        {line.getPosAtDist(roots[0]),
         line.getPosAtDist(roots[1])});
  }

  CUDA_CALLABLE PointLocation
  ptLocation(const Point<dim, fptype> &pt,
             fptype absPrecision = defAbsPrecision) const {
    assert(absPrecision >= 0.0);
    fptype ptPos = evaluatePoint<fptype>(pt);
    if(ptPos < -absPrecision)
      return PT_INSIDE;
    else if(MathFuncs::MathFuncs<fptype>::fabs(ptPos) <
            absPrecision)
      return PT_ON;
    else
      return PT_OUTSIDE;
  }

  CUDA_CALLABLE Vector<dim, fptype> normal(
      const Point<dim, fptype> &pt) {
    Vector<dim, fptype> n;
    for(int i = 0; i < dim; i++) {
      fptype tmp = coeff(i, i) * pt.getOffset().get(i);
      for(int j = 0; j < dim; j++) {
        tmp += coeff(i, j) * pt.getOffset().get(j);
      }
      tmp += coeff(i, dim);
      n.set(i, tmp);
    }
    return n;
  }

  CUDA_CALLABLE Quadric<dim, fptype> operator=(
      const Quadric<dim, fptype> &q) {
    Solid<dim, fptype>::operator=(q);
    for(int i = 0; i < numCoeffs; i++) {
      setCoeff(i, q.coeff(i));
    }
    return *this;
  }

  CUDA_CALLABLE Quadric<dim, fptype> operator=(
      const Quadric<dim, fptype>::QuadricData &q) {
    this->origin = q.origin;
    coeffs = q.coeffs;
    return *this;
  }

  CUDA_CALLABLE bool operator==(
      const Quadric<dim, fptype> &q) const {
    for(int i = 0; i < numCoeffs; i++) {
      if(coeff(i) != q.coeff(i)) return false;
    }
    return true;
  }

  template <typename srctype>
  CUDA_CALLABLE bool operator==(
      const Quadric<dim, srctype> &q) const {
    using hPrec =
        typename Typecast::typecast<fptype,
                                    srctype>::higherPrec;
    for(int i = 0; i < numCoeffs; i++) {
      if(static_cast<hPrec>(coeff(i)) !=
         static_cast<hPrec>(q.coeff(i)))
        return false;
    }
    return true;
  }

  template <typename srctype>
  CUDA_CALLABLE bool operator!=(
      const Quadric<dim, srctype> &q) const {
    return !((*this) == q);
  }

#ifdef __CUDACC__
  std::shared_ptr<QuadricData> cudaCopy() const {
    QuadricData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    auto safeMem =
        std::shared_ptr<QuadricData>(cudaMem, cudaFree);
    err = cudaCopy(safeMem);
    return safeMem;
  }

  cudaError_t cudaCopy(
      std::shared_ptr<QuadricData> cudaMem) const {
    cudaError_t err =
        cudaMemcpy(&cudaMem->coeffs, &coeffs,
                   sizeof(coeffs), cudaMemcpyHostToDevice);
    err = this->origin.cudaCopy(&cudaMem->origin);
    return err;
  }

  cudaError_t cudaCopy(QuadricData *cudaMem) const {
    cudaError_t err =
        cudaMemcpy(&cudaMem->coeffs, &coeffs,
                   sizeof(coeffs), cudaMemcpyHostToDevice);
    err = this->origin.cudaCopy(&cudaMem->origin);
    return err;
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<QuadricData> cudaMem) {
    QuadricData buffer;
    cudaError_t err = cudaMemcpy(
        &buffer.coeffs, &cudaMem->coeffs,
        sizeof(buffer.coeffs), cudaMemcpyDeviceToHost);
    for(int i = 0; i < numCoeffs; i++)
      setCoeff(i, buffer.coeffs[i]);
    return err;
  }

  cudaError_t cudaRetrieve(QuadricData *cudaMem) {
    QuadricData buffer;
    cudaError_t err = cudaMemcpy(
        &buffer.coeffs, &cudaMem->coeffs,
        sizeof(buffer.coeffs), cudaMemcpyDeviceToHost);
    for(int i = 0; i < numCoeffs; i++)
      setCoeff(i, buffer.coeffs[i]);
    return err;
  }
#endif

  bool readFileMatrix(std::istream &is) {
    for(int i = 0; i < dim + 1; i++) {
      for(int j = 0; j < dim + 1; j++) {
        if(is.eof()) return false;
        fptype val;
        is >> val;
        if(!MathFuncs::MathFuncs<fptype>::isnan(
               coeff(i, j)))
          val += coeff(i, j);
        setCoeff(i, j, val);
      }
    }
    return true;
  }

  bool readFilePolynomial(std::istream &is) {
    for(int i = 0; i < numCoeffs; i++) {
      if(is.eof()) return false;
      fptype val;
      is >> val;
      setCoeff(i, val);
    }
    return true;
  }

  bool readFile(std::istream &is) {
    enum QuadricFormats {
      FMT_MATRIX = 'm',
      FMT_POLY = 'p',
    };
    char qFormat;
    is >> qFormat;
    switch(qFormat) {
      case FMT_MATRIX:
        return readFileMatrix(is);
      case FMT_POLY:
        return readFilePolynomial(is);
      default:
        return false;
    }
  }

  friend std::istream &operator>>(std::istream &is,
                                  Quadric<dim, fptype> &q) {
    q.readFile(is);
    return is;
  }

  friend std::ostream &operator<<(
      std::ostream &os, const Quadric<dim, fptype> &q) {
    bool once = false;
    for(unsigned i = 0; i < dim + 1; i++) {
      for(unsigned j = i; j < dim + 1; j++) {
        if(q.coeff(i, j) != 0) {
          if(once) {
            os << " + ";
          } else {
            once = !once;
          }
          os << q.coeff(i, j);
          if(i < dim) {
            os << " * x" << i;
            if(j < dim) {
              os << " * x" << j;
            }
          }
        }
      }
    }
    os << " = 0";
    return os;
  }

  template <int, typename>
  friend class Quadric;

  static constexpr const int numCoeffs =
      (dim + 2) * (dim + 1) / 2;

 private:
  CUDA_CALLABLE static int getCoeffPos(int d1, int d2) {
    /* Quadric Coefficient Ordering:
     * c1 x1^2 + c2 x2^2 + ... + cd +
     * c(d+1) x1 x2 + c(d+2) x1 x3 + ... + c(2d) x1 +
     * c(2d+1) x2 x3 + c(2d+2) x2 x4 + ... c(3d-1) x2 + ...
     */
    assert(0 <= d1);
    assert(d1 < dim + 1);
    assert(0 <= d2);
    assert(d2 < dim + 1);
    if(d1 == d2) {
      return d1;
    } else {
      int first = d1 < d2 ? d1 : d2;
      int second = d1 < d2 ? d2 : d1;
      int offset =
          first * dim - (first * (first - 1)) / 2 + dim + 1;
      return offset + second - first - 1;
    }
  }

  Array<fptype, numCoeffs> coeffs;
};
};

#endif
