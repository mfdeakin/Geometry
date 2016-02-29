
#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "geometry.hpp"

#include <assert.h>
#include <memory>
#include <ostream>

#include "array.hpp"
#include "mathfuncs.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Vector : public GeometryBase<dim, fptype> {
 public:
  struct VectorData {
    Array<fptype, dim> vector;
    // For efficiency,
    // this is only computed when explicitly asked for.
    // Otherwise, the first value will contain an NAN
    Array<Array<fptype, dim>, dim - 1> orthonorms;
  };

  CUDA_CALLABLE Vector() {
    for(int i = 0; i < dim; i++) {
      offset.vector[i] = 0.0;
    }
    offset.orthonorms[0][0] = NAN;
  }

  CUDA_CALLABLE Vector(const Array<fptype, dim> &v) {
    offset.vector = v;
    offset.orthonorms[0][0] = NAN;
  }

  CUDA_CALLABLE Vector(const VectorData &v) {
    offset.vector = v.vector;
    if(!MathFuncs::MathFuncs<fptype>::isnan(
           v.orthonorms[0][0]))
      offset.orthonorms = v.orthonorms;
    else
      offset.orthonorms[0][0] = fptype(NAN);
  }

  template <typename srctype>
  CUDA_CALLABLE Vector(const Vector<dim, srctype> &src) {
    for(int i = 0; i < dim; i++) {
      offset.vector[i] = fptype(src.offset.vector[i]);
      /* We can't assume that the orthogonal vectors can
       * just be cast to the current type, unfortunately.
       * Thus, initialize them to NAN
       */
      for(int j = 0; j < dim - 1; j++) {
        offset.orthonorms[j][i] = fptype(NAN);
      }
    }
  }

  CUDA_CALLABLE fptype get(int dimension) const {
    return offset.vector[dimension];
  }

  CUDA_CALLABLE fptype set(int dimension, fptype val) {
    offset.vector[dimension] = val;
    return offset.vector[dimension];
  }

  CUDA_CALLABLE Vector<dim, fptype> add(
      const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> sum;
    for(int i = 0; i < dim; i++)
      sum.offset.vector[i] =
          offset.vector[i] + rhs.offset.vector[i];
    return sum;
  }

  CUDA_CALLABLE Vector<dim, fptype> subtract(
      const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> diff;
    for(int i = 0; i < dim; i++)
      diff.offset.vector[i] =
          offset.vector[i] - rhs.offset.vector[i];
    return diff;
  }

  CUDA_CALLABLE fptype operator()(int dimension) const {
    return get(dimension);
  }

  CUDA_CALLABLE fptype operator()(int dimension,
                                  fptype val) {
    return set(dimension, val);
  }

  CUDA_CALLABLE Vector<dim, fptype> operator*(
      fptype scale) const {
    return this->scale(scale);
  }

  CUDA_CALLABLE Vector<dim, fptype> operator+(
      const Vector<dim, fptype> &rhs) const {
    return add(rhs);
  }

  CUDA_CALLABLE Vector<dim, fptype> operator-(
      const Vector<dim, fptype> &rhs) const {
    return subtract(rhs);
  }

  CUDA_CALLABLE Vector<dim, fptype> operator*=(
      fptype scale) {
    for(unsigned i = 0; i < dim; i++)
      offset.vector[i] *= scale;
    return *this;
  }

  CUDA_CALLABLE Vector<dim, fptype> operator+=(
      const Vector<dim, fptype> &rhs) {
    for(unsigned i = 0; i < dim; i++)
      set(i, get(i) + rhs.get(i));
    offset.orthonorms[0][0] = fptype(NAN);
    return *this;
  }

  CUDA_CALLABLE Vector<dim, fptype> operator-=(
      const Vector<dim, fptype> &rhs) {
    for(unsigned i = 0; i < dim; i++)
      set(i, get(i) - rhs.get(i));
    offset.orthonorms[0][0] = fptype(NAN);
    return *this;
  }

  CUDA_CALLABLE bool operator==(
      const Vector<dim, fptype> &rhs) const {
    for(unsigned i = 0; i < dim; i++)
      if(get(i) != rhs.get(i)) return false;
    return true;
  }

  CUDA_CALLABLE bool operator!=(
      const Vector<dim, fptype> &rhs) const {
    return !((*this) == rhs);
  }

  CUDA_CALLABLE Vector<dim, fptype> operator=(
      const Vector<dim, fptype> &v) {
    offset.vector = v.offset.vector;
    if(v.orthosComputed())
      offset.orthonorms = v.offset.orthonorms;
    else
      offset.orthonorms[0][0] = NAN;
    return *this;
  }

  CUDA_CALLABLE Vector<dim, fptype> operator=(
      const VectorData &v) {
    offset.vector = v.vector;
    if(v.orthosComputed())
      offset.orthonorms = v.orthonorms;
    else
      offset.orthonorms[0][0] = NAN;
    return *this;
  }

  CUDA_CALLABLE VectorData copyData() const {
    VectorData v;
    v.vector = offset.vector;
    if(orthosComputed())
      v.orthonorms = offset.orthonorms;
    else
      v.orthonorms[0][0] = NAN;
    return v;
  }

  CUDA_CALLABLE VectorData &copyData(VectorData &v) const {
    v.vector = offset.vector;
    if(orthosComputed())
      v.orthonorms = offset.orthonorms;
    else
      v.orthonorms[0][0] = NAN;
    return v;
  }

  CUDA_CALLABLE fptype
  dot(const Vector<dim, fptype> &rhs) const {
    fptype sum = 0.0;
    for(int i = 0; i < dim; i++)
      sum = MathFuncs::MathFuncs<fptype>::fma(
          get(i), rhs.get(i), sum);
    return sum;
  }

  template <int otherDim>
  CUDA_CALLABLE Vector<dim + otherDim - 1, fptype> convZPad(
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

  CUDA_CALLABLE Vector<dim, fptype> normalize() const {
    fptype mag = norm();
    if(mag == 0.0)
      // Can't normalize a 0 vector
      return *this;
    return scale(fptype(1.0) / mag);
  }

  CUDA_CALLABLE fptype normSquare() const {
    fptype normSq = 0.0;
    for(int i = 0; i < dim; i++) {
      normSq += offset.vector[i] * offset.vector[i];
    }
    return normSq;
  }

  CUDA_CALLABLE fptype norm() const {
    return MathFuncs::MathFuncs<fptype>::sqrt(normSquare());
  }

  CUDA_CALLABLE Vector<dim, fptype> scale(
      fptype scalar) const {
    if(scalar == 0.0) {
      return Vector<dim, fptype>();
    } else {
      Vector<dim, fptype> scaled(*this);
      for(int i = 0; i < dim; i++)
        scaled.set(i, scaled.get(i) * scalar);
      return scaled;
    }
  }

  CUDA_CALLABLE bool orthosComputed() const {
    return !MathFuncs::MathFuncs<fptype>::isnan(
        offset.orthonorms[0][0]);
  }

  CUDA_CALLABLE Vector<dim, fptype> getOrthogonal(
      int orthoIdx) {
    assert(orthoIdx >= 0);
    assert(orthoIdx < dim - 1);
    if(!orthosComputed()) {
      calcOrthogonals();
    }
    return returnOrthogonal(orthoIdx);
  }

  friend std::ostream &operator<<(
      std::ostream &os, const Vector<dim, fptype> &v) {
    os << "(";
    if(v.dim > 0) os << v.offset.vector[0];
    for(unsigned i = 1; i < dim; i++)
      os << ", " << v.offset.vector[i];
    os << ")";
    return os;
  }

  CUDA_CALLABLE void copy(VectorData *mem) const {
    copyData(*mem);
  }

  CUDA_CALLABLE void copy(
      std::shared_ptr<VectorData> mem) const {
    copyData(*mem);
  }

#ifdef __CUDACC__
  std::shared_ptr<VectorData> cudaCopy() const {
    VectorData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    cudaMemcpy(cudaMem, &offset, sizeof(offset),
               cudaMemcpyHostToDevice);
    return std::shared_ptr<VectorData>(cudaMem, cudaFree);
  }

  cudaError_t cudaCopy(VectorData *cudaMem) const {
    return cudaMemcpy(cudaMem, &offset, sizeof(offset),
                      cudaMemcpyHostToDevice);
  }

  cudaError_t cudaCopy(
      std::shared_ptr<VectorData> cudaMem) const {
    return cudaMemcpy(cudaMem.get(), &offset,
                      sizeof(offset),
                      cudaMemcpyHostToDevice);
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<VectorData> cudaMem) {
    return cudaMemcpy(&offset, cudaMem.get(),
                      sizeof(offset),
                      cudaMemcpyDeviceToHost);
  }

  cudaError_t cudaRetrieve(VectorData *cudaMem) {
    return cudaMemcpy(&offset, cudaMem, sizeof(offset),
                      cudaMemcpyDeviceToHost);
  }
#endif

  template <int, typename>
  friend class Vector;

 private:
  CUDA_CALLABLE Vector<dim, fptype> returnOrthogonal(
      int orthoIdx) const {
    return Vector<dim, fptype>(offset.orthonorms[orthoIdx]);
  }

  CUDA_CALLABLE void calcOrthogonals() {
    /* Use a Householder reflection to compute the
     * orthonormal vectors.
     * This works because the Householder matrix is
     * unitary.
     * See the following answer for more details:
     * https://math.stackexchange.com/questions/710103/algorithm-to-find-an-orthogonal-basis-orthogonal-to-a-given-vector
     */
    fptype n = MathFuncs::MathFuncs<fptype>::copysign(
        norm(), get(0));
    assert(n != 0.0);
    Array<fptype, dim> w(offset.vector);
    w[0] += n;
    fptype wNormSq = 0.0;
    for(int i = 0; i < dim; i++) {
      wNormSq = MathFuncs::MathFuncs<fptype>::fma(
          w[i], w[i], wNormSq);
      for(int j = 0; j < dim - 1; j++) {
        offset.orthonorms[j][i] = 0.0;
      }
    }
    for(int i = 0; i < dim - 1; i++) {
      offset.orthonorms[i][i + 1] = 1.0;
      for(int j = 0; j < dim; j++) {
        fptype updated = MathFuncs::MathFuncs<fptype>::fma(
            2 * w[i + 1], -w[j] / wNormSq,
            offset.orthonorms[i][j]);
        offset.orthonorms[i][j] = updated;
      }
    }
  }

  VectorData offset;
};
};

#endif
