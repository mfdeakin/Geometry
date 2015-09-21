
#ifndef _ORIGIN_HPP_
#define _ORIGIN_HPP_

#include "geometry.hpp"
#include "vector.hpp"

#include <array>

namespace Geometry {

template <int dim, typename fptype>
class Origin : public GeometryBase<dim, fptype> {
 public:
  Origin() {}

  Origin(const Origin<dim, fptype> &src)
      : globalCoords(src.globalCoords) {}

  Origin(const std::array<fptype, dim> &globalPos)
      : globalCoords(globalPos) {}

  template <typename srctype>
  bool operator==(const Origin<dim, srctype> &other) const {
    return globalOffset() == other.globalOffset();
  }

  Vector<dim, fptype> calcOffset(
      const Origin<dim, fptype> &other) const {
    return globalOffset() - other.globalOffset();
  }

  virtual Vector<dim, fptype> globalOffset() const {
    return globalCoords;
  }

  template <int, typename>
  friend class Origin;

 private:
  Vector<dim, fptype> globalCoords;
};
}

#endif
