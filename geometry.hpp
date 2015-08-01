
#ifndef _GEOMETRY_HPP_
#define _GEOMETRY_HPP_

namespace Geometry {

template <unsigned, typename>
class Point;

enum PointLocation {
  PT_INSIDE = -1,
  PT_ON = 0,
  PT_OUTSIDE = 1
};

template <unsigned _dim, typename fptype>
class Geometry {
public:

  virtual ~Geometry() {};

  virtual PointLocation
  ptLocation(const Point<_dim, float> &test,
             fptype absPrecision = 9.5367431640625e-7) = 0;
  
  virtual PointLocation
  ptLocation(const Point<_dim, double> &test,
             fptype absPrecision = 9.5367431640625e-7) = 0;
  
  const unsigned dim = _dim;
};

};

#endif
