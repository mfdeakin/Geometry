
#ifndef _ORIGIN_HPP_
#define _ORIGIN_HPP_

#include "geometry.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Origin : public Geometry<dim, fptype> {
 public:
  Origin() {
    for(int i = 0; i < dim; i++) globalCoords[i] = 0.0;
  }

  template <typename srctype>
  Origin(const Origin<dim, srctype> &src) {
    for(int i = 0; i < dim; i++)
      globalCoords[i] = src.globalCoords[i];
  }

  bool operator==(const Origin &other) const {
    bool isEq = true;
    for(int i = 0; i < dim; i++)
      isEq &= (globalCoords[i] == other.globalCoords[i]);
    return isEq;
  }

  template <int, typename>
  friend class Origin;
 private:
  fptype globalCoords[dim];
};
}

#endif