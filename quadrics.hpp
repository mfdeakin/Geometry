
#ifndef _QUADRICS_HPP_
#define _QUADRICS_HPP_

#include <stdlib.h>
#include <assert.h>

#include "origin.hpp"
#include "point.hpp"
#include "line.hpp"

namespace Geometry {

template <unsigned dim, typename fptype>
class Quadric : public Solid<dim, fptype> {
 public:
  Quadric(const Origin<dim, fptype> &origin)
      : Solid<dim, fptype>(origin) {
    initialCoeffs = new fptype[numCoeffs];
    currentCoeffs = new fptype[numCoeffs];
    refs = new unsigned;
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

  fptype getCoeff(int d1, int d2) {}

  template <typename ptfptype>
  fptype evaluatePoint(
      const Point<dim, ptfptype> &pt,
      fptype absPrecision = defAbsPrecision) const {
    /* Use Kahan summation to evaluate this more correctly
     */
    fptype ret = 0.0;
    fptype extra = 0.0;
    for(unsigned i = 0; i < dim; i++) {
      for(unsigned j = 0; j < dim; j++) {
        unsigned coeffNum = std::fma(i, dim, j);
        fptype product =
            currentCoeffs[coeffNum] * pt(i) * pt(j) - extra;
        fptype tmp = ret + product;
        extra = (tmp - ret) - product;
        ret = tmp;
      }
    }
    return ret;
  }

  template <typename linefptype>
  fptype calcLineIntersect(
      const Line<dim, linefptype> &line,
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
      const Point<dim, float> &pt,
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

  PointLocation ptLocation(const Point<dim, double> &pt,
                           fptype absPrecision) const {
    double ptPos = evaluatePoint(pt);
    if(ptPos < 0)
      return PT_INSIDE;
    else if(ptPos == 0)
      return PT_ON;
    else
      return PT_OUTSIDE;
  }

 private:
  int getCoeffPos(int d1, int d2) {
    if(d1 == d2)
      return d1;
    else if(d2 >= 0 && d2 < dim)
      return d1 * dim - d1 - ((d1 - 1) * d1) / 2 - 1 +
             (d2 - d1) + dim;
    else if(d1 >= 0 && d1 < dim)
      return dim * (dim + 1) / 2 + d1;
    else
      return numCoeffs - 1;
  }

  static constexpr unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;

  fptype *initialCoeffs;
  fptype *currentCoeffs;
  unsigned *refs;
};
};

#endif
