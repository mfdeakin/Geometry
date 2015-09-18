
#ifndef _QUADRICS_HPP_
#define _QUADRICS_HPP_

#include <stdlib.h>
#include <assert.h>

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
      : Solid<dim, fptype>(origin) {
    initialCoeffs = new fptype[numCoeffs];
    currentCoeffs = new fptype[numCoeffs];
    refs = new int;
    *refs = 1;
  }

  Quadric(const Quadric &src) {
    refs = src.refs;
    (*refs)++;
    initialCoeffs = src.initialCoeffs;
    currentCoeffs = src.currentCoeffs;
  }

  virtual ~Quadric() {
    (*refs)--;
    if(*refs == 0) {
      delete initialCoeffs;
      delete currentCoeffs;
      delete refs;
      initialCoeffs = NULL;
      currentCoeffs = NULL;
      refs = NULL;
    }
  }

  /* d=0 corresponds to the first dimension,
   * d=1 the second, ...
   * d=dim the constant coefficient
   */
  fptype &coeff(int d1, int d2) const {
    int coeffNum = getCoeffPos(d1, d2);
    return currentCoeffs[coeffNum];
  }

  /* d=0 corresponds to the first dimension,
   * d=1 the second, ...
   * d=dim the constant coefficient
   */
  fptype &coeff(int pos) const {
    assert(pos >= 0);
    assert(pos < numCoeffs);
    return currentCoeffs[pos];
  }

  fptype evaluatePoint(
      const Point<dim, fptype> &pt,
      fptype absPrecision = defAbsPrecision) const {
    /* Use Kahan summation to evaluate this more correctly
     */
    fptype ret = 0.0;
    fptype extra = 0.0;
    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < dim; j++) {
        int coeffNum = getCoeffPos(i, j);
        fptype product =
            currentCoeffs[coeffNum] * pt(i) * pt(j) - extra;
        fptype tmp = ret + product;
        extra = (tmp - ret) - product;
        ret = tmp;
      }
    }
    return ret;
  }

  fptype calcLineIntersect(
      const Line<dim, fptype> &line,
      fptype absPrecision = defAbsPrecision) const {
    auto ld = line.getDir();
    fptype sqCoeff = 0;
    for(int i = 0; i < dim; i++) {
      sqCoeff += currentCoeffs[i] * ld(i) * ld(i);
      for(int j = i + 1; j < dim; j++) {
        int coeffPos = getCoeffPos(i, j);
        sqCoeff += currentCoeffs[coeffPos] * ld(i) * ld(j);
      }
    }
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

  fptype *initialCoeffs;
  fptype *currentCoeffs;
  int *refs;
};
};

#endif
