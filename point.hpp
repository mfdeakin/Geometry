
#ifndef _POINT_HPP_
#define _POINT_HPP_

#include "geometry.hpp"

#include <assert.h>
#include <typeinfo>
#include <cmath>

#include "accurate_math.hpp"

namespace Geometry {

template <unsigned _dim, typename fptype>
class Point : public Geometry<_dim, fptype> {
 public:
  Point() {
    for(unsigned i = 0; i < _dim; i++) {
      setPos(i, 0.0);
    }
  };

  virtual ~Point(){};

  fptype getPos(unsigned dimension) const {
    assert(dimension < this->dim);
    return pos[dimension];
  }

  fptype setPos(unsigned dimension, fptype value) {
    assert(dimension < this->dim);
    pos[dimension] = value;
    return pos[dimension];
  }

  fptype operator()(unsigned dimension) const {
    return getPos(dimension);
  }

  fptype operator()(unsigned dimension, fptype pos) {
    return setPos(dimension, pos);
  }

  Point<_dim, fptype> subtract(
      const Point<_dim, fptype> &rhs) {
    Point<_dim, fptype> diff;
    for(unsigned i = 0; i < _dim; i++) {
      fptype delta = getPos(i) - rhs.getPos(i);
      diff(i, delta);
    }
    return diff;
  }

  Point<_dim, fptype> add(const Point<_dim, fptype> &rhs) {
    Point<_dim, fptype> sum;
    for(unsigned i = 0; i < _dim; i++) {
      fptype total = getPos(i) + rhs.getPos(i);
      sum(i, total);
    }
    return sum;
  }

  fptype distToOrigin() {
    fptype sum = kahanSum<fptype, _dim, twoNormSum>(pos);
    return sqrt(sum);
  }

  PointLocation ptLocation(
      const Point<_dim, float> &pt,
      fptype absPrecision = 9.5367431640625e-7) {
    // We have no knowledge of the precision, so downcast
    Point<_dim, float> diff = subtract(pt);
    float dist = diff.distToOrigin();
    if(dist > absPrecision)
      return PT_OUTSIDE;
    else
      return PT_INSIDE;
  }

  PointLocation ptLocation(
      const Point<_dim, double> &pt,
      fptype absPrecision = 9.5367431640625e-7) {
    bool isEq = true;
    if(typeid(fptype) == typeid(float)) {
      for(unsigned i = 0; i < this->dim && isEq; i++) {
        isEq = isEq && (pos[i] == ((float)pt.pos[i]));
      }
    } else {
      for(unsigned i = 0; i < this->dim && isEq; i++) {
        isEq = isEq && (pos[i] == pt.pos[i]);
      }
    }
    if(isEq) {
      return PT_INSIDE;
    } else {
      return PT_OUTSIDE;
    }
  }

  template <unsigned, typename>
  friend class Geometry;
  template <unsigned, typename>
  friend class Point;

 protected:
  fptype pos[_dim];
};
};

#endif
