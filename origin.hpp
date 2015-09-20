
#ifndef _ORIGIN_HPP_
#define _ORIGIN_HPP_

#include "geometry.hpp"
#include "vector.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Origin : public GeometryBase<dim, fptype> {
 public:
  Origin() {
    for(int i = 0; i < dim; i++) globalCoords[i] = 0.0;
  }

  template <typename srctype>
  Origin(const Origin<dim, srctype> &src) {
    for(int i = 0; i < dim; i++)
      globalCoords[i] = src.globalCoords[i];
  }

  template <typename srctype>
  bool operator==(const Origin<dim, srctype> &other) const {
    bool isEq = true;
    for(int i = 0; i < dim; i++)
      isEq &= (globalCoords[i] == other.globalCoords[i]);
    return isEq;
  }

  Vector<dim, fptype> offset(const Origin<dim, fptype> &other) const {
    Vector<dim, fptype> off;
    for(int i = 0; i < dim; i++) {
      fptype delta = globalCoords[i] - other.globalCoords[i];
      off(i, delta);
    }
    return off;
  }
  
  template <int, typename>
  friend class Origin;

 private:
  fptype globalCoords[dim];
};
}

#endif
