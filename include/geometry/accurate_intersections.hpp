
#ifndef _ACCURATE_INTERSECTIONS_HPP_
#define _ACCURATE_INTERSECTIONS_HPP_

#include "genericfp.hpp"
#include "geometry.hpp"
#include "line.hpp"
#include "polynomial.hpp"
#include "quadric.hpp"

#include <cmath>
#include <iostream>
#include <list>
#include <memory>
#include <utility>

namespace Geometry {

template <int dim, typename fptype, typename cmpAlg>
class IntersectionBase {
 public:
  Quadric<dim, fptype> q;
  Line<dim, fptype> l;
  fptype intPos, otherIntPos;
  fptype absErrMargin;
  mutable int numIP;

  IntersectionBase()
      : q(),
        intPos(NAN),
        otherIntPos(NAN),
        absErrMargin(defAbsPrecision),
        numIP(0) {}

  IntersectionBase(const Quadric<dim, fptype> &quad,
                   const Line<dim, fptype> &line,
                   fptype intPos = NAN,
                   fptype otherIntPos = NAN,
                   fptype absErrMargin = defAbsPrecision)
      : q(quad),
        l(line),
        intPos(intPos),
        otherIntPos(otherIntPos),
        absErrMargin(absErrMargin),
        numIP(0) {}

  template <typename otherCmpAlg>
  IntersectionBase(
      const IntersectionBase<dim, fptype, otherCmpAlg> &i)
      : q(i.q),
        l(i.l),
        intPos(i.intPos),
        otherIntPos(i.otherIntPos),
        absErrMargin(i.absErrMargin),
        numIP(i.numIP) {}

  template <typename otherCmpAlg>
  IntersectionBase<dim, fptype, cmpAlg> operator=(
      const IntersectionBase<dim, fptype, otherCmpAlg> &i) {
    q = i.q;
    l = i.l;
    intPos = i.intPos;
    otherIntPos = i.otherIntPos;
    absErrMargin = i.absErrMargin;
    numIP = i.numIP;
    return *this;
  }

  bool isUndetermined(const fptype &i1,
                      const fptype &i2) const {
    fptype norm1 = l.getDirection().norm();
    fptype norm2 = l.getDirection().norm();
    return MathFuncs::MathFuncs<fptype>::fabs(i1 * norm2 -
                                              i2 * norm1) <
           absErrMargin * norm1 * norm2;
  }

  int incPrecCount() const { return numIP; }

  static bool cmp(
      const IntersectionBase<dim, fptype, cmpAlg> &lhs,
      const IntersectionBase<dim, fptype, cmpAlg> &rhs) {
    fptype diff = lhs.compare(rhs);
    if(diff < 0) return true;
    return false;
  }

  template <typename otherCmpAlg>
  fptype compare(
      const IntersectionBase<dim, fptype, otherCmpAlg> &i)
      const {
    fptype delta = intPos - i.intPos;
    if(isUndetermined(intPos, i.intPos) && i.q != q &&
       !MathFuncs::MathFuncs<fptype>::isnan(otherIntPos) &&
       !MathFuncs::MathFuncs<fptype>::isnan(
           i.otherIntPos)) {
      numIP++;
      return static_cast<cmpAlg>(*this).accurateCompare(
          static_cast<cmpAlg>(i));
    } else {
      return delta;
    }
  }
};

template <int dim, typename fptype, typename cmpAlg>
std::shared_ptr<std::list<cmpAlg>> sortIntersections(
    const Line<dim, fptype> &line,
    std::list<Quadric<dim, fptype>> &quads,
    const fptype &absErrMargin = defAbsPrecision) {
  std::shared_ptr<std::list<cmpAlg>> inter(
      new std::list<cmpAlg>);
  for(auto q : quads) {
    auto intPos = q.calcLineDistToIntersect(line);
    for(int k = 0; k < 2; k++) {
      if(!MathFuncs::MathFuncs<fptype>::isnan(intPos[k])) {
        cmpAlg i(q, line, intPos[k], intPos[1 - k],
                 absErrMargin);
        inter->push_back(i);
      }
    }
  }
  inter->sort([](cmpAlg &lhs, cmpAlg &rhs) {
    return lhs.compare(rhs) < 0.0;
  });
  /* TODO: Fix this */
  // quicksortInt(inter->begin(), inter->end());
  return inter;
}
}

#endif
