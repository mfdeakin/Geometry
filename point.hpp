
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
class Point : public Solid<dim, fptype>,
              public Origin<dim, fptype> {
 public:
  template <typename srctype>
  Point(const Origin<dim, srctype> &origin,
        const Vector<dim, srctype> &offset)
      : Solid<dim, fptype>(origin), offset(offset) {}

  template <typename srctype>
  Point(const Point<dim, srctype> &src)
      : Solid<dim, fptype>(src.origin),
        offset(src.offset) {}

  virtual ~Point(){};

  virtual void shiftOrigin(
      const Origin<dim, fptype> &newOrigin) {
    Vector<dim, fptype> v =
        this->origin.calcOffset(newOrigin);
    Solid<dim, fptype>::shiftOrigin(newOrigin);
  }

  virtual Origin<dim, fptype> getOrigin() const {
    return this->origin;
  }

  virtual Vector<dim, fptype> getOffset() const {
    return offset;
  }

  Point<dim, fptype> addVector(
      Vector<dim, fptype> v) const {
    Vector<dim, fptype> newOffset;
    for(int i = 0; i < dim; i++) {
      fptype newCoord = v(i) + offset(i);
      newOffset.set(i, newCoord);
    }
    return Point<dim, fptype>(this->origin, newOffset);
  }

  fptype distToOrigin() const { return offset.norm(); }

  PointLocation ptLocation(
      const Point<dim, fptype> &pt,
      fptype absPrecision = defAbsPrecision) const {
    // We have no knowledge of the precision, so downcast
    Vector<dim, float> delta = this->calcOffset(pt);
    float dist = delta.norm();
    if(dist > absPrecision)
      return PT_OUTSIDE;
    else
      return PT_INSIDE;
  }

  virtual Vector<dim, fptype> globalOffset() const {
    Vector<dim, fptype> o =
        offset + this->origin.globalOffset();
    return o;
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
