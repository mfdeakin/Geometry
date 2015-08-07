
#ifndef _POINT_HPP_
#define _POINT_HPP_

#include "geometry.hpp"
#include "origin.hpp"
#include "vector.hpp"
#include "line.hpp"

#include <assert.h>
#include <typeinfo>
#include <cmath>

#include "accurate_math.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Point : public Solid<dim, fptype> {
 public:
  Point(const Origin<dim, fptype> &origin,
        const Vector<dim, fptype> &offset)
      : Solid<dim, fptype>(origin), offset(offset) {}

  Point(const Point<dim, fptype> &src)
      : Solid<dim, fptype>(src.origin),
        offset(src.offset) {}

  virtual ~Point(){};

  fptype getPos(int dimension) const {
    assert(dimension < this->dim);
    return NAN;
  }

  fptype setPos(int dimension, fptype value) {
    assert(dimension < this->dim);
    return NAN;
  }

  fptype operator()(int dimension) const {
    return getPos(dimension);
  }

  fptype operator()(int dimension, fptype pos) {
    return setPos(dimension, pos);
  }

  Vector<dim, fptype> pointOffset(
      const Point<dim, fptype> &other) const {
    assert(other.origin == this->origin);
    return offset - other.offset;
  }

  fptype distToOrigin() const { return offset.norm(); }

  PointLocation ptLocation(
      const Point<dim, float> &pt,
      fptype absPrecision = defAbsPrecision) const {
    // We have no knowledge of the precision, so downcast
    Vector<dim, float> delta = pointOffset(pt);
    float dist = delta.norm();
    if(dist > absPrecision)
      return PT_OUTSIDE;
    else
      return PT_INSIDE;
  }

  PointLocation ptLocation(
      const Point<dim, double> &pt,
      fptype absPrecision = defAbsPrecision) const {
    bool isEq = true;
    if(isEq) {
      return PT_INSIDE;
    } else {
      return PT_OUTSIDE;
    }
  }

  template <int, typename>
  friend class Line;

  template <int, typename>
  friend class Point;

 protected:
  Vector<dim, fptype> offset;
};
};

#endif
