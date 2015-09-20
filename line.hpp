
#ifndef _LINE_HPP_
#define _LINE_HPP_

#include <stdlib.h>
#include <assert.h>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "vector.hpp"
#include "typecast.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Line : public Solid<dim, fptype> {
 public:
  Line(const Point<dim, fptype> &intercept,
       const Vector<dim, fptype> &direction)
      : Solid<dim, fptype>(
            static_cast<const Origin<dim, fptype> &>(
                intercept)),
        intercept(intercept, Vector<dim, fptype>(
                                 std::array<fptype, dim>(
                                     {{0.0, 0.0, 0.0}}))),
        dir(direction) {}

  Line(const Line &src)
      : Solid<dim, fptype>(src.origin),
        intercept(src.intercept),
        dir(src.dir) {}

  virtual ~Line() {}

  virtual void shiftOrigin(
      const Origin<dim, fptype> &newOrigin) {
    /* Shift the origin by changing the intercept to be the
     * minimal perpendicular to the direction.
     * Do this by computing the perpendicular direction
     * of the offset and it's magnitude.
     * This leaves the only possible catastrophic
     * cancellation in the computation of the offset */
    Vector<dim, fptype> offset =
        newOrigin.offset(this->origin);
    auto perpDirs = dir.calcOrthogonals();
    Solid<dim, fptype>::shiftOrigin(newOrigin);
  }

  virtual PointLocation ptLocation(
      const Point<dim, fptype> &test,
      fptype absPrecision = defAbsPrecision) const {
    auto ptDir = test.pointOffset(intercept);
    fptype offsetLen = ptDir.dot(ptDir);
    fptype dist = ptDir.dot(dir);
    fptype perpDist = std::abs(offsetLen - dist * dist);
    if(perpDist < absPrecision)
      return PT_INSIDE;
    else if(perpDist == absPrecision)
      return PT_ON;
    else
      return PT_OUTSIDE;
  }

  Vector<dim, fptype> getDirection() const { return dir; }

  template <int, typename>
  friend class Line;

 protected:
  Point<dim, fptype> intercept;
  Vector<dim, fptype> dir;
};
};

#endif
