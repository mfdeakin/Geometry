
#ifndef _QUADRICS_HPP_
#define _QUADRICS_HPP_

#include <stdlib.h>
#include <assert.h>

#include <ostream>
#include <iostream>
#include <memory>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "line.hpp"

#include "accurate_math.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Quadric : public Solid<dim, fptype> {
 public:
  Quadric(const Origin<dim, fptype> &origin)
      : Solid<dim, fptype>(origin),
        coeffs(new fptype[numCoeffs],
               std::default_delete<fptype[]>()) {
    for(int i = 0; i < numCoeffs; i++) coeff(i) = 0.0;
  }

  Quadric(const Quadric &src) : coeffs(src.currentCoeffs) {}

  virtual ~Quadric() {}

  /* d=0 corresponds to the first dimension,
   * d=1 the second, ...
   * d=dim the constant coefficient
   */
  fptype &coeff(int d1, int d2) const {
    int coeffNum = getCoeffPos(d1, d2);
    coeffs.get()[coeffNum] = coeffs.get()[coeffNum];
    return coeffs.get()[coeffNum];
  }

  /* Access the coefficients as a 1D array.
   * pos=0 to dim-1 are for the square terms
   * pos=dim is the constant term
   * pos=dim+1 to 2dim-1 is the first dim prod terms
   * pos=2dim is the first dim linear term
   * ...
   * pos=(dim+2)*(dim+1)/2-1 is the last dim linear term
   */
  fptype &coeff(int pos) const {
    assert(pos >= 0);
    assert(pos < numCoeffs);
    return coeffs.get()[pos];
  }

  fptype evaluatePoint(
      const Point<dim, fptype> &pt,
      fptype absPrecision = defAbsPrecision) const {
    /* Use Kahan summation to evaluate this more correctly
     */
    fptype ret = 0.0;
    fptype extra = 0.0;
    Point<dim, fptype> copy(pt);
    copy.shiftOrigin(this->origin);
    Vector<dim, fptype> quadOffset = copy.getOffset();
    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        int coeffNum = getCoeffPos(i, j);
        fptype product = coeffs.get()[coeffNum] *
                             quadOffset(i) * quadOffset(j) -
                         extra;
        fptype tmp = ret + product;
        extra = (tmp - ret) - product;
        ret = tmp;
      }
    }
    return ret;
  }

  std::array<Point<dim, fptype>, 2> calcLineIntersect(
      Line<dim, fptype> line,
      fptype absPrecision = defAbsPrecision) const {
    line.shiftOrigin(this->origin);
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
    std::cout << sqCoeff << " x^2 + " << linCoeff << " x + "
              << constant << "\n";
    std::array<fptype, 2> roots =
        AccurateMath::kahanQuadratic(sqCoeff, linCoeff,
                                     constant);
    return std::array<Point<dim, fptype>, 2>(
        {{line.getPosAtDist(roots[0]),
          line.getPosAtDist(roots[1])}});
  }

  PointLocation ptLocation(
      const Point<dim, fptype> &pt,
      fptype absPrecision = defAbsPrecision) const {
    assert(absPrecision >= 0.0);
    double ptPos = evaluatePoint(pt);
    if(ptPos < -absPrecision)
      return PT_INSIDE;
    else if(std::fabs(ptPos) < absPrecision)
      return PT_ON;
    else
      return PT_OUTSIDE;
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

 private:
  static constexpr int getCoeffPos(int d1, int d2) {
    assert(0 <= d1);
    assert(d1 < dim + 1);
    assert(0 <= d2);
    assert(d2 < dim + 1);
    if(d1 == d2) {
      return d1;
    } else {
      int first = std::min(d1, d2);
      int second = std::max(d1, d2);
      int offset =
          first * dim - (first * (first - 1)) / 2 + dim + 1;
      return offset + second - first - 1;
    }
  }

  static constexpr const int numCoeffs =
      (dim + 2) * (dim + 1) / 2;

  std::shared_ptr<fptype> coeffs;
};
};

#endif
