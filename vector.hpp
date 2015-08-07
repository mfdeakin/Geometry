
#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "geometry.hpp"

#include <assert.h>
#include <array>
#include <cmath>

namespace Geometry {

template <int dim, typename fptype>
class Vector : public Geometry<dim, fptype> {
 public:
  Vector() {
    for(int i = 0; i < dim; i++) offset[i] = 0.0;
  }

  Vector(const Vector<dim, fptype> &src) {
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
    fptype totalProd = 1.0;
    fptype newMag = 0.0;
    int subtractIdx = -1;
    int vecNum = 0;
    std::array<Vector<dim, fptype>, dim - 1> basis;
    /* O(dim) computation of a basis for orthogonal vectors
     */
    for(int i = 0; i < dim; i++) {
      if(offset[i] != 0.0) {
        totalProd *= offset[i];
        subtractIdx = i;
        newMag += 1.0 / offset[i] / offset[i];
      } else {
        basis[vecNum].offset[i] = 1.0;
        vecNum++;
      }
    }
    /* Ignore degenerate cases when all components are 0 */
    assert(subtractIdx != -1);
    for(int vComp = 0; vecNum < dim - 1; vComp++) {
      if(offset[vComp] == 0.0) continue;
      fptype vecMag =
          std::sqrt(newMag +
                    (dim - 2) * (dim - 2) / offset[vComp] /
                        offset[vComp]) * totalProd;
      for(int i = 0; i < dim; i++) {
        if(i != vComp && offset[i] != 0.0)
          basis[vecNum].offset[i] =
              totalProd / offset[i] / vecMag;
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
