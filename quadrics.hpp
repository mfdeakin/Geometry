
#ifndef _QUADRICS_HPP_
#define _QUADRICS_HPP_

#include <stdlib.h>
#include <assert.h>

#include <ostream>
#include <istream>
#include <iostream>
#include <array>
#include <memory>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "line.hpp"

#include "polynomial.hpp"

#include "accurate_math.hpp"

#include "mathfuncs.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Quadric : public Solid<dim, fptype> {
 public:
  struct QuadricData {
    typename Origin<dim, fptype>::OriginData origin;
    std::array<fptype, Quadric<dim, fptype>::numCoeffs>
        coeffs;
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
      : Solid<dim, fptype>(src.origin),
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

  CUDA_CALLABLE fptype
      setCoeff(int d1, int d2, fptype val) {
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

  CUDA_CALLABLE Polynomial<2, fptype> calcLineDistPoly(
      const Line<dim, fptype> &line,
      fptype absPrecision = defAbsPrecision) const {
    auto lDir = line.getDirection();
    auto lInt = line.getIntercept().getOffset();
    fptype sqCoeff = 0.0;
    fptype linCoeff = 0.0;
    fptype constant = coeff(dim, dim);
    for(int i = 0; i < dim; i++) {
      sqCoeff += coeff(i, i) * lDir(i) * lDir(i);
      linCoeff += 2 * coeff(i, i) * lDir(i) * lInt(i);
      linCoeff += coeff(i, dim) * lDir(i);
      constant += coeff(i, dim) * lInt(i);
      constant += coeff(i, i) * lInt(i) * lInt(i);
      for(int j = i + 1; j < dim; j++) {
        sqCoeff += coeff(i, j) * lDir(i) * lDir(j);
        linCoeff += coeff(i, j) * lDir(i) * lInt(j);
        linCoeff += coeff(i, j) * lDir(j) * lInt(i);
        constant += coeff(i, j) * lInt(i) * lInt(j);
      }
    }
    Polynomial<2, fptype> poly;
    poly.set(0, constant);
    poly.set(1, linCoeff);
    poly.set(2, sqCoeff);
    return poly;
  }

  CUDA_CALLABLE std::array<fptype, 2>
  calcLineDistToIntersect(
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

  CUDA_CALLABLE std::array<Point<dim, fptype>, 2>
  calcLineIntersect(
      const Line<dim, fptype> &line,
      fptype absPrecision = defAbsPrecision) const {
    std::array<fptype, 2> roots =
        calcLineDistToIntersect(line, absPrecision);
    return std::array<Point<dim, fptype>, 2>(
        {line.getPosAtDist(roots[0]),
         line.getPosAtDist(roots[1])});
  }

  CUDA_CALLABLE PointLocation ptLocation(
      const Point<dim, fptype> &pt,
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

  CUDA_CALLABLE Quadric<dim, fptype> operator=(
      const Quadric<dim, fptype> &q) {
    Solid<dim, fptype>::operator=(q);
    for(int i = 0; i < numCoeffs; i++)
      setCoeff(i, q.coeff(i));
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
    using hPrec = typename Typecast::typecast<
        fptype, srctype>::higherPrec;
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
    err = cudaCopy(cudaMem);
    return std::shared_ptr<QuadricData>(cudaMem, cudaFree);
  }

  cudaError_t cudaCopy(QuadricData *cudaMem) const {
    cudaError_t err =
        cudaMemcpy(&cudaMem->coeffs, &coeffs,
                   sizeof(coeffs), cudaMemcpyHostToDevice);
    err = this->origin.cudaCopy(&cudaMem->origin);
    return err;
  }
#endif

  bool readFileMatrix(std::istream &is) {
    for(int i = 0; i < dim + 1; i++) {
      for(int j = 0; j < dim + 1; j++) {
        if(is.eof()) return false;
        fptype val;
        is >> val;
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
        if(once)
          os << " + ";
        else
          once = !once;
        os << q.coeff(i, j);
        if(i < dim) {
          os << " * x" << i;
          if(j < dim) os << " * x" << j;
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
    assert(0 <= d1);
    assert(d1 < dim + 1);
    assert(0 <= d2);
    assert(d2 < dim + 1);
    if(d1 == d2) {
      return d1;
    } else {
      int first = MathFuncs::MathFuncs<fptype>::min(d1, d2);
      int second =
          MathFuncs::MathFuncs<fptype>::max(d1, d2);
      int offset =
          first * dim - (first * (first - 1)) / 2 + dim + 1;
      return offset + second - first - 1;
    }
  }

  std::array<fptype, numCoeffs> coeffs;
};
};

#endif
