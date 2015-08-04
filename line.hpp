
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
      : Solid<dim, fptype>(intercept.origin),
        intercept(intercept),
        dir(direction) {}

  Line(const Line &src)
      : Solid<dim, fptype>(src.origin),
        intercept(src.intercept),
        dir(src.dir) {}

  virtual ~Line() {}

  virtual PointLocation ptLocation(
      const Point<dim, float> &test,
      fptype absPrecision = defAbsPrecision) {
    return PT_OUTSIDE;
  }

  virtual PointLocation ptLocation(
      const Point<dim, double> &test,
      fptype absPrecision = defAbsPrecision) {
    return PT_OUTSIDE;
  }

  template <int, typename>
  friend class Line;

 protected:
  Point<dim, fptype> intercept;
  Vector<dim, fptype> dir;
};
};

#endif
