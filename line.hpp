
#ifndef _LINE_HPP_
#define _LINE_HPP_

#include <stdlib.h>
#include <assert.h>
#include "geometry.hpp"
#include "accurate_math.hpp"

namespace Geometry {

template <unsigned _dim, typename fptype>
class Line : public Geometry<_dim, fptype> {
public:
  Line() {
  }
  
  virtual ~Line() {
  }
  
  Line(const Line &src) {
    currentOrigin = src.currentOrigin;
  }

  virtual PointLocation
  ptLocation(const Point<_dim, float> &test) {
    return PT_OUTSIDE;
  }
  
  virtual PointLocation
  ptLocation(const Point<_dim, double> &test) {
    return PT_OUTSIDE;
  }
  
  template<unsigned, typename>
  friend class Line;
protected:
  Point<_dim, fptype> points[2];
  Point<_dim, fptype> currentOrigin;
};

};

#endif
