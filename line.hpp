
#ifndef _LINE_HPP_
#define _LINE_HPP_

#include <stdlib.h>
#include <assert.h>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "vector.hpp"

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
