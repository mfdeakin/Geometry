
#ifndef _POINT_HPP_
#define _POINT_HPP_

#include "geometry.hpp"
#include "origin.hpp"
#include "vector.hpp"
#include "line.hpp"

#include <assert.h>
#include <typeinfo>
#include <cmath>

#include "accurate_math.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Point : public Solid<dim, fptype>,
              public Origin<dim, fptype> {
 public:
  struct PointData {
    typename Origin<dim, fptype>::OriginData o;
    typename Vector<dim, fptype>::VectorData v;
  };

  CUDA_CALLABLE Point() : Solid<dim, fptype>() {
    for(int i = 0; i < dim; i++) offset.set(i, NAN);
  }

  template <typename srctype>
  CUDA_CALLABLE Point(const Vector<dim, srctype> &offset)
      : Solid<dim, fptype>(), offset(offset) {}

  template <typename srctype>
  CUDA_CALLABLE Point(const Origin<dim, srctype> &origin,
                      const Vector<dim, srctype> &offset)
      : Solid<dim, fptype>(origin), offset(offset) {}

  template <typename srctype>
  CUDA_CALLABLE Point(const Point<dim, srctype> &src)
      : Solid<dim, fptype>(src.origin),
        offset(src.offset) {}

  CUDA_CALLABLE virtual ~Point(){};

  CUDA_CALLABLE virtual void shiftOrigin(
      const Origin<dim, fptype> &newOrigin) {
    Vector<dim, fptype> v =
        this->origin.calcOffset(newOrigin);
    Solid<dim, fptype>::shiftOrigin(newOrigin);
  }

  CUDA_CALLABLE virtual const Origin<dim, fptype>
      &getOrigin() const {
    return this->origin;
  }

  CUDA_CALLABLE virtual const Vector<dim, fptype>
      &getOffset() const {
    return offset;
  }

  CUDA_CALLABLE Point<dim, fptype> addVector(
      Vector<dim, fptype> v) const {
    Vector<dim, fptype> newOffset;
    for(int i = 0; i < dim; i++) {
      fptype newCoord = v(i) + offset(i);
      newOffset.set(i, newCoord);
    }
    return Point<dim, fptype>(this->origin, newOffset);
  }

  CUDA_CALLABLE Vector<dim, fptype> ptDiff(
      Point<dim, fptype> rhs) {
    rhs.shiftOrigin(this->origin);
    return offset - rhs.offset;
  }

  CUDA_CALLABLE fptype distToOrigin() const {
    return offset.norm();
  }

  CUDA_CALLABLE PointLocation
  ptLocation(const Point<dim, fptype> &pt,
             fptype absPrecision = defAbsPrecision) const {
    // We have no knowledge of the precision, so downcast
    Vector<dim, float> delta = this->calcOffset(pt);
    float dist = delta.norm();
    if(dist > absPrecision)
      return PT_OUTSIDE;
    else
      return PT_INSIDE;
  }

  CUDA_CALLABLE virtual Vector<dim, fptype> globalOffset()
      const {
    Vector<dim, fptype> o =
        offset + this->origin.globalOffset();
    return o;
  }

  CUDA_CALLABLE Point<dim, fptype> operator=(
      Point<dim, fptype> p) {
    Solid<dim, fptype>::operator=(p);
    offset = p.offset;
    return *this;
  }

  CUDA_CALLABLE Point<dim, fptype> operator=(PointData p) {
    Solid<dim, fptype>::operator=(p.origin);
    offset = p.v;
    return *this;
  }

  CUDA_CALLABLE PointData copyData() const {
    PointData p;
    this->origin.copyData(p.o);
    offset.copyData(p.v);
    return p;
  }

  CUDA_CALLABLE PointData &copyData(PointData &p) const {
    this->origin.copyData(p.o);
    offset.copyData(p.v);
    return p;
  }

#ifdef __CUDACC__
  std::shared_ptr<PointData> cudaCopy() const {
    PointData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    auto safeMem =
        std::shared_ptr<PointData>(cudaMem, cudaFree);
    err = cudaCopy(safeMem);
    return safeMem;
  }

  cudaError_t cudaCopy(
      std::shared_ptr<PointData> cudaMem) const {
    cudaError_t err = this->origin.cudaCopy(&cudaMem->o);
    err = this->offset.cudaCopy(&cudaMem->v);
    return err;
  }

  cudaError_t cudaCopy(PointData *cudaMem) const {
    cudaError_t err = this->origin.cudaCopy(&cudaMem->o);
    err = this->offset.cudaCopy(&cudaMem->v);
    return err;
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<PointData> cudaMem) {
    cudaError_t err =
        this->origin.cudaRetrieve(&cudaMem->o);
    err = this->offset.cudaRetrieve(&cudaMem->v);
    return err;
  }

  cudaError_t cudaRetrieve(PointData *cudaMem) {
    cudaError_t err =
        this->origin.cudaRetrieve(&cudaMem->o);
    err = this->offset.cudaRetrieve(&cudaMem->v);
    return err;
  }
#endif

  friend std::ostream &operator<<(
      std::ostream &os, const Point<dim, fptype> &p) {
    auto pos = p.globalOffset();
    os << "<";
    if(p.dim > 0) os << pos(0);
    for(unsigned i = 1; i < dim; i++) os << ", " << pos(i);
    os << ">";
    return os;
  }

  template <int, typename>
  friend class Line;

  template <int, typename>
  friend class Point;

 protected:
  Vector<dim, fptype> offset;
};
};

#endif
