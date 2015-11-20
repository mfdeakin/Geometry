
#ifndef _GEOMETRY_HPP_
#define _GEOMETRY_HPP_

#include "cudadef.h"

namespace Geometry {

enum PointLocation {
  PT_INSIDE = -1,
  PT_ON = 0,
  PT_OUTSIDE = 1
};

template <int, typename>
class Origin;

template <int, typename>
class Vector;

template <int, typename>
class Solid;

template <int, typename>
class Point;

template <int, typename>
class Line;

template <int dim, typename fptype>
class Quadric;

constexpr float defAbsPrecision = 9.5367431640625e-7;

template <int _dim, typename fptype>
class GeometryBase {
 public:
  CUDA_CALLABLE virtual ~GeometryBase(){};

  static_assert(_dim >= 0,
                "The dimension of a geometric object "
                "cannot be negative!");
  static constexpr const int dim = _dim;
};

/* A solid is a well defined geometric object which
 * is positioned relative to an origin,
 * and which has a defined interior and exterior,
 * and possibly surface
 */
template <int dim, typename fptype>
class Solid : public GeometryBase<dim, fptype> {
 public:
  CUDA_CALLABLE Solid()
      : origin(Origin<dim, fptype>::uOrigin()) {}

  template <typename srctype>
  CUDA_CALLABLE Solid(const Solid<dim, srctype> &s)
      : origin(s.origin) {}

  template <typename srctype>
  CUDA_CALLABLE Solid(const Origin<dim, srctype> &o)
      : origin(o) {}

  CUDA_CALLABLE virtual ~Solid(){};

  CUDA_CALLABLE Solid<dim, fptype> &operator=(
      const Solid<dim, fptype> &s) {
    origin = s.origin;
    return *this;
  }

  CUDA_CALLABLE virtual PointLocation ptLocation(
      const Point<dim, fptype> &test,
      fptype absPrecision = defAbsPrecision) const = 0;

  CUDA_CALLABLE virtual void shiftOrigin(
      const Origin<dim, fptype> &newOrigin) {
    origin = newOrigin;
  }

  CUDA_CALLABLE Origin<dim, fptype> getOrigin() const {
    return origin;
  }

  template <int, typename>
  friend class Solid;

 protected:
  Origin<dim, fptype> origin;
};
};

#endif
