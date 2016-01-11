
#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "geometry.hpp"

#include <assert.h>
#include <array>
#include <memory>
#include <ostream>

#include "array.hpp"
#include "mathfuncs.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Vector : public GeometryBase<dim, fptype> {
 public:
  typedef Array<fptype, dim> VectorData;

  CUDA_CALLABLE Vector() {
    for(int i = 0; i < dim; i++) offset[i] = 0.0;
  }

  CUDA_CALLABLE Vector(const VectorData &v) {
    for(int i = 0; i < dim; i++) offset[i] = v[i];
  }

  template <typename srctype>
  CUDA_CALLABLE Vector(const Vector<dim, srctype> &src) {
    for(int i = 0; i < dim; i++)
      offset[i] = (fptype)src.offset[i];
  }

  CUDA_CALLABLE fptype get(int dimension) const {
    assert(dimension >= 0);
    assert(dimension < dim);
    return offset[dimension];
  }

  CUDA_CALLABLE fptype set(int dimension, fptype val) {
    assert(dimension >= 0);
    assert(dimension < dim);
    offset[dimension] = val;
    return offset[dimension];
  }

  CUDA_CALLABLE Vector<dim, fptype> add(
      const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> sum;
    for(int i = 0; i < dim; i++)
      sum.offset[i] = offset[i] + rhs.offset[i];
    return sum;
  }

  CUDA_CALLABLE Vector<dim, fptype> subtract(
      const Vector<dim, fptype> &rhs) const {
    Vector<dim, fptype> diff;
    for(int i = 0; i < dim; i++)
      diff.offset[i] = offset[i] - rhs.offset[i];
    return diff;
  }

  CUDA_CALLABLE fptype operator()(int dimension) const {
    return get(dimension);
  }

  CUDA_CALLABLE fptype
  operator()(int dimension, fptype val) {
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
      set(i, get(i) * scale);
    return *this;
  }

  CUDA_CALLABLE Vector<dim, fptype> operator+=(
      const Vector<dim, fptype> &rhs) {
    for(unsigned i = 0; i < dim; i++)
      set(i, get(i) + rhs.get(i));
    return *this;
  }

  CUDA_CALLABLE Vector<dim, fptype> operator-=(
      const Vector<dim, fptype> &rhs) {
    for(unsigned i = 0; i < dim; i++)
      set(i, get(i) - rhs.get(i));
    return *this;
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
    Vector<dim, fptype> normalized;
    if(mag == 0.0)
      // Can't normalize a 0 vector
      return normalized;
    for(int i = 0; i < dim; i++)
      normalized.offset[i] = offset[i] / mag;
    return normalized;
  }

  CUDA_CALLABLE fptype norm() const {
    fptype normSq = dot(*this);
    return MathFuncs::MathFuncs<fptype>::sqrt(normSq);
  }

  CUDA_CALLABLE Vector<dim, fptype> scale(
      fptype scalar) const {
    Vector<dim, fptype> scaled;
    for(int i = 0; i < dim; i++)
      scaled.offset[i] = offset[i] * scalar;
    return scaled;
  }

  CUDA_CALLABLE std::array<Vector<dim, fptype>, dim - 1>
  calcOrthogonals() const {
    /* Use a Householder reflection to compute the
     * orthonormal vectors.
     * This works because the Householder matrix is unitary.
     * See the following answer for more details:
     * https://math.stackexchange.com/questions/710103/algorithm-to-find-an-orthogonal-basis-orthogonal-to-a-given-vector
     */
    fptype n = MathFuncs::MathFuncs<fptype>::copysign(
        norm(), offset[0]);
    assert(n != 0.0);
    Vector<dim, fptype> w(*this);
    w.set(0, n + offset[0]);
    const fptype wNormSq = w.dot(w);
    std::array<Vector<dim, fptype>, dim - 1> basis;
    for(int i = 0; i < dim - 1; i++) {
      basis[i].offset[i + 1] = 1.0;
      for(int j = 0; j < dim; j++) {
        fptype updated = MathFuncs::MathFuncs<fptype>::fma(
            2 * w.get(i + 1), -w.get(j) / wNormSq,
            basis[i].offset[j]);
        basis[i].set(j, updated);
      }
    }
    return basis;
  }

  friend std::ostream &operator<<(
      std::ostream &os, const Vector<dim, fptype> &v) {
    os << "(";
    if(v.dim > 0) os << v.offset[0];
    for(unsigned i = 1; i < dim; i++)
      os << ", " << v.offset[i];
    os << ")";
    return os;
  }

#ifdef __CUDACC__
  std::shared_ptr<const VectorData> cudaCopy() const {
    VectorData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    cudaMemcpy(cudaMem, &offset, sizeof(offset),
               cudaMemcpyHostToDevice);
    return std::shared_ptr<const VectorData>(cudaMem,
                                             cudaFree);
  }

  cudaError_t cudaCopy(VectorData *cudaMem) const {
    return cudaMemcpy(cudaMem, &offset, sizeof(offset),
                      cudaMemcpyHostToDevice);
  }

  cudaError_t cudaCopy(
											 std::shared_ptr<VectorData> cudaMem) const {
    return cudaMemcpy(cudaMem.get(), &offset, sizeof(offset),
											cudaMemcpyHostToDevice);
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<VectorData> cudaMem) {
		return cudaMemcpy(&offset, cudaMem.get(), sizeof(offset),
											cudaMemcpyDeviceToHost);
  }

  cudaError_t cudaRetrieve(
      VectorData *cudaMem) {
		return cudaMemcpy(&offset, cudaMem, sizeof(offset),
											cudaMemcpyDeviceToHost);
  }
#endif

  template <int, typename>
  friend class Vector;

 private:
  VectorData offset;
};
};

#endif
