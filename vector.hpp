
#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "geometry.hpp"

#include <cmath>

namespace Geometry {

template <int dim, typename fptype>
class Vector : public Geometry<dim, fptype> {
public:
  Vector() {
    for(int i = 0; i < dim; i++)
      offset[i] = 0.0;
  }
  
  Vector(const Vector<dim, fptype> &src) {
    for(int i = 0; i < dim; i++)
      offset[i] = src.offset[i];
  }
  
  Vector<dim, fptype> add(const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> sum;
    for(int i = 0; i < dim; i++)
      sum.offset[i] = offset[i] + rhs.offset[i];
    return sum;
  }
  
  Vector<dim, fptype> subtract(const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> diff;
    for(int i = 0; i < dim; i++)
      diff.offset[i] = offset[i] - rhs.offset[i];
    return diff;
  }
  
  Vector<dim, fptype> operator -(const Vector<dim, fptype> &rhs) const {
    return subtract(rhs);
  }
  
  Vector<dim, fptype> operator +(const Vector<dim, fptype> &rhs) const {
    return add(rhs);
  }
  
  fptype dot(const Vector<dim, fptype> &rhs) const {
    fptype sum = 0.0;
    for(int i = 0; i < dim; i++)
      sum = std::fma(offset[i], rhs.offset[i], sum);
    return sum;
  }
  
  Vector<dim, fptype> normalize() const {
    fptype mag = norm();
    Vector<dim, fptype> normalized;
    for(int i = 0; i < dim; i++)
      normalized.offset[i] = offset[i] / mag;
  }
  
  fptype norm() const {
    fptype normSq = dot(*this);
    return std::sqrt(normSq);
  }

  Vector<dim, fptype> scale(fptype scalar) const {
    Vector<dim, fptype> scaled;
    for(int i = 0; i < dim; i++)
      scaled.offset[i] = offset[i] * scalar;
    return scaled;
  }
  
private:
  fptype offset[dim];
};

};

#endif
