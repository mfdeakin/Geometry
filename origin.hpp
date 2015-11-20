
#ifndef _ORIGIN_HPP_
#define _ORIGIN_HPP_

#include "geometry.hpp"
#include "vector.hpp"

#include <array>

namespace Geometry {

template <int dim, typename fptype>
class Origin : public GeometryBase<dim, fptype> {
 public:
  struct OriginData {
    typename Vector<dim, fptype>::VectorData v;
  };

  CUDA_CALLABLE Origin() {}

  CUDA_CALLABLE Origin(const Origin<dim, fptype> &src)
      : globalCoords(src.globalCoords) {}

  template <typename srctype>
  CUDA_CALLABLE Origin(const Origin<dim, srctype> &src)
      : globalCoords(src.globalCoords) {}

  CUDA_CALLABLE Origin(
      const std::array<fptype, dim> &globalPos)
      : globalCoords(globalPos) {}

  template <typename srctype>
  CUDA_CALLABLE bool operator==(
      const Origin<dim, srctype> &other) const {
    return globalOffset() == other.globalOffset();
  }

  CUDA_CALLABLE Vector<dim, fptype> calcOffset(
      const Origin<dim, fptype> &other) const {
    return globalOffset() - other.globalOffset();
  }

  CUDA_CALLABLE virtual Vector<dim, fptype> globalOffset()
      const {
    return globalCoords;
  }

  template <int d, typename f>
  friend class Origin;

  static const Origin<dim, fptype> &uOrigin() {
    static const Origin<dim, fptype> universeOrigin;
    return universeOrigin;
  }

#ifdef __CUDACC__
  std::shared_ptr<const OriginData> cudaCopy() const {
    OriginData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    err = globalCoords.cudaCopy(&cudaMem->v);
    return std::shared_ptr<const OriginData>(cudaMem,
                                             cudaFree);
  }

  cudaError_t cudaCopy(OriginData *cudaMem) const {
    return globalCoords.cudaCopy(&cudaMem->v);
  }
#endif

 private:
  Vector<dim, fptype> globalCoords;
};
}

#endif
