
#ifndef _QUADRICS_HPP_
#define _QUADRICS_HPP_

#include <stdlib.h>
#include <assert.h>

#include "point.hpp"

namespace Geometry {

template <unsigned _dim, typename fptype>
class Quadric : public Geometry<_dim, fptype> {
 public:
  const static unsigned numCoeffs =
      (_dim + 2) * (_dim + 1) / 2;

  Quadric() {
    initialCoeffs = new fptype[numCoeffs];
    currentCoeffs = new fptype[numCoeffs];
    refs = new unsigned;
    *refs = 1;
  }

  Quadric(const Quadric &src) {
    refs = src.refs;
    (*refs)++;
    currentOrigin = src.currentOrigin;
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

  template <typename ptfptype>
  fptype evaluatePoint(const Point<_dim, ptfptype> &pt,
                       fptype absPrecision = defAbsPrecision) {
    /* Use Kahan summation to evaluate this more correctly
     */
    fptype ret = 0.0;
    fptype extra = 0.0;
    for(unsigned i = 0; i < _dim; i++) {
      for(unsigned j = 0; j < _dim; j++) {
        unsigned coeffNum = std::fma(i, _dim, j);
        fptype product =
            currentCoeffs[coeffNum] * pt(i) * pt(j) - extra;
        fptype tmp = ret + product;
        extra = (tmp - ret) - product;
        ret = tmp;
      }
    }
    return ret;
  }
  
  PointLocation ptLocation(const Point<_dim, float> &pt,
                           fptype absPrecision = defAbsPrecision) {
    assert(absPrecision >= 0.0);
    double ptPos = evaluatePoint(pt);
    if(ptPos < -absPrecision)
      return PT_INSIDE;
    else if(std::fabs(ptPos) < absPrecision)
      return PT_ON;
    else
      return PT_OUTSIDE;
  }

  PointLocation ptLocation(const Point<_dim, double> &pt,
                           fptype absPrecision) {
    double ptPos = evaluatePoint(pt);
    if(ptPos < 0)
      return PT_INSIDE;
    else if(ptPos == 0)
      return PT_ON;
    else
      return PT_OUTSIDE;
  }

 private:
  Point<_dim, fptype> currentOrigin;
  fptype *initialCoeffs;
  fptype *currentCoeffs;
  unsigned *refs;
};
};

#endif
