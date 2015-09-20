
#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "geometry.hpp"

#include <assert.h>
#include <array>
#include <cmath>

namespace Geometry {

template <int dim, typename fptype>
class Vector : public GeometryBase<dim, fptype> {
 public:
  Vector() {
    for(int i = 0; i < dim; i++) offset[i] = 0.0;
  }

  Vector(std::array<fptype, dim> v) {
    for(int i = 0; i < dim; i++) offset[i] = v[i];
  }

  template <typename srctype>
  Vector(const Vector<dim, srctype> &src) {
    for(int i = 0; i < dim; i++) offset[i] = src.offset[i];
  }

  fptype get(int dimension) const {
    assert(dimension >= 0);
    assert(dimension < dim);
    return offset[dimension];
  }

  fptype set(int dimension, fptype val) {
    assert(dimension >= 0);
    assert(dimension < dim);
    offset[dimension] = val;
    return offset[dimension];
  }

  Vector<dim, fptype> add(
      const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> sum;
    for(int i = 0; i < dim; i++)
      sum.offset[i] = offset[i] + rhs.offset[i];
    return sum;
  }

  Vector<dim, fptype> subtract(
      const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> diff;
    for(int i = 0; i < dim; i++)
      diff.offset[i] = offset[i] - rhs.offset[i];
    return diff;
  }

  fptype operator()(int dimension) const {
    return get(dimension);
  }

  fptype operator()(int dimension, fptype val) {
    return set(dimension, val);
  }

  Vector<dim, fptype> operator*(fptype scale) {
    return this->scale(scale);
  }

  Vector<dim, fptype> operator+(
      const Vector<dim, fptype> &rhs) const {
    return add(rhs);
  }

  Vector<dim, fptype> operator-(
      const Vector<dim, fptype> &rhs) const {
    return subtract(rhs);
  }

  fptype dot(const Vector<dim, fptype> &rhs) const {
    fptype sum = 0.0;
    for(int i = 0; i < dim; i++)
      sum = std::fma(get(i), rhs.get(i), sum);
    return sum;
  }

  template <int otherDim>
  Vector<dim + otherDim - 1, fptype> convZPad(
      const Vector<otherDim, fptype> &rhs) {
    Vector<dim + otherDim - 1, fptype> convd;
    for(int i = 0; i < dim; i++) {
      for(int j = 0; j < otherDim; j++) {
        fptype csum =
            std::fma(get(i), rhs.get(j), convd.get(i + j));
        convd.set(i + j, csum);
      }
    }
    return convd;
  }

  Vector<dim, fptype> normalize() const {
    fptype mag = norm();
    Vector<dim, fptype> normalized;
    for(int i = 0; i < dim; i++)
      normalized.offset[i] = offset[i] / mag;
    return normalized;
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

  std::array<Vector<dim, fptype>, dim - 1> calcOrthogonals()
      const {
    std::array<Vector<dim, fptype>, dim - 1> basis;
    int firstNonZero = -1;
    for(int d = 0, vecNum = 0; d < dim; d++) {
      if(offset[d] == 0.0) {
        basis[vecNum].offset[d] = 1.0;
        vecNum++;
      }
      else if(firstNonZero == -1) {
        firstNonZero = d;
      }
      else {
        basis[vecNum].offset[d] = -offset[firstNonZero];
        basis[vecNum].offset[firstNonZero] = offset[d];
        vecNum++;
      }
    }
    return basis;
  }

  template <int, typename>
  friend class Vector;

 private:
  fptype offset[dim];
};
};

#endif
