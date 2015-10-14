
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

#include "polynomial.hpp"

#include "accurate_math.hpp"

#include "mathfuncs.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Quadric : public Solid<dim, fptype> {
 public:
  Quadric(const Origin<dim, fptype> &origin)
      : Solid<dim, fptype>(origin),
        coeffs(new fptype[numCoeffs],
               std::default_delete<fptype[]>()) {
    for(int i = 0; i < numCoeffs; i++) setCoeff(i, 0.0);
  }

  Quadric(const Quadric<dim, fptype> &src)
      : Solid<dim, fptype>(src.origin),
        coeffs(src.coeffs) {}

  template <typename srctype>
  Quadric(const Quadric<dim, srctype> &src)
      : Solid<dim, fptype>(src.origin),
        coeffs(new fptype[numCoeffs],
               std::default_delete<fptype[]>()) {
    for(int i = 0; i < numCoeffs; i++)
      setCoeff(i, src.coeff(i));
  }

  virtual ~Quadric() {}

  /* d=0 corresponds to the first dimension,
   * d=1 the second, ...
   * d=dim the constant coefficient
   */
  fptype coeff(int d1, int d2) const {
    int coeffNum = getCoeffPos(d1, d2);
    return coeffs.get()[coeffNum];
  }

  fptype setCoeff(int d1, int d2, fptype val) {
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
  fptype coeff(int pos) const {
    assert(pos >= 0);
    assert(pos < numCoeffs);
    return coeffs.get()[pos];
  }

  fptype setCoeff(int pos, fptype val) {
    assert(pos >= 0);
    assert(pos < numCoeffs);
    if(coeffs.unique() == false) {
      fptype *newCoeffs = new fptype[numCoeffs];
      for(int i = 0; i < numCoeffs; i++)
        newCoeffs[i] = coeff(i);
      coeffs.reset(newCoeffs,
                   std::default_delete<fptype[]>());
    }
    coeffs.get()[pos] = val;
    return coeffs.get()[pos];
  }

  template <typename outtype>
  outtype evaluatePoint(
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
        outtype product = coeffs.get()[coeffNum] *
                              quadOffset(i) *
                              quadOffset(j) -
                          extra;
        outtype tmp = ret + product;
        extra = (tmp - ret) - product;
        ret = tmp;
      }
    }
    return ret;
  }

  Polynomial<2, fptype> calcLineDistPoly(
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

  std::array<fptype, 2> calcLineDistToIntersect(
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

  std::array<Point<dim, fptype>, 2> calcLineIntersect(
      const Line<dim, fptype> &line,
      fptype absPrecision = defAbsPrecision) const {
    std::array<fptype, 2> roots =
        calcLineDistToIntersect(line, absPrecision);
    return std::array<Point<dim, fptype>, 2>(
        {line.getPosAtDist(roots[0]),
         line.getPosAtDist(roots[1])});
  }

  PointLocation ptLocation(
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

  Quadric<dim, fptype> operator=(
      const Quadric<dim, fptype> &q) {
    Solid<dim, fptype>::operator=(q);
    coeffs = q.coeffs;
    return *this;
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

 private:
  static constexpr int getCoeffPos(int d1, int d2) {
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

  static constexpr const int numCoeffs =
      (dim + 2) * (dim + 1) / 2;

  std::shared_ptr<fptype> coeffs;
};
};

#endif
