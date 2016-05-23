
#ifndef _LINE_HPP_
#define _LINE_HPP_

#include <stdlib.h>
#include <assert.h>
#include <ostream>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "vector.hpp"

#include "typecast.hpp"

#include "mathfuncs.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Line : public Solid<dim, fptype> {
 public:
  struct LineData {
    typename Vector<dim, fptype>::VectorData d;
    typename Point<dim, fptype>::PointData p;
  };

  CUDA_CALLABLE Line() : Solid<dim, fptype>(), intercept() {
    for(int i = 0; i < dim; i++) dir.set(i, NAN);
  }

  CUDA_CALLABLE Line(const Point<dim, fptype> &intercept,
                     const Vector<dim, fptype> &direction)
      : Solid<dim, fptype>(intercept.origin),
        intercept(intercept),
        dir(direction) {
    assert(dir.norm() != 0.0);
  }

  template <typename srctype>
  CUDA_CALLABLE Line(const Line<dim, srctype> &src)
      : Solid<dim, fptype>(src.origin),
        intercept(src.intercept),
        dir(src.dir) {}

  CUDA_CALLABLE virtual ~Line() {}

  CUDA_CALLABLE virtual void shiftOrigin(
      const Origin<dim, fptype> &newOrigin) {
    /* Shift the origin by changing the intercept to be the
     * minimal perpendicular to the direction.
     * Do this by computing the perpendicular direction
     * of the offset and it's magnitude.
     * This leaves the only possible catastrophic
     * cancellation in the computation of the offset */
    Vector<dim, fptype> delta(
        intercept.calcOffset(newOrigin));
    Vector<dim, fptype> interceptOff;
    for(int i = 0; i < dim - 1; i++) {
      Vector<dim, fptype> p =
          dir.getOrthogonal(i).normalize();
      fptype scale = p.dot(delta);
      interceptOff += p * scale;
    }
    intercept = Point<dim, fptype>(newOrigin, interceptOff);
    Solid<dim, fptype>::shiftOrigin(newOrigin);
  }

  CUDA_CALLABLE virtual fptype argPointMinDist(
      const Point<dim, fptype> &pt) {
    Vector<dim, fptype> delta(pt.calcOffset(intercept));
    return delta.dot(dir);
  }

  CUDA_CALLABLE virtual PointLocation ptLocation(
      const Point<dim, fptype> &test,
      fptype absPrecision = defAbsPrecision) const {
    Vector<dim, fptype> ptDir(test.calcOffset(intercept));
    fptype offsetLen = ptDir.dot(ptDir);
    fptype dist = ptDir.dot(dir);
    fptype perpDist = MathFuncs::MathFuncs<fptype>::fabs(
        offsetLen - dist * dist);
    if(perpDist < absPrecision)
      return PT_INSIDE;
    else if(perpDist == absPrecision)
      return PT_ON;
    else
      return PT_OUTSIDE;
  }

  CUDA_CALLABLE Point<dim, fptype> getPosAtDist(
      fptype dist) const {
    Vector<dim, fptype> pDist;
    const auto p0 = intercept.getOffset();
    for(int i = 0; i < dim; i++) {
      fptype coord = std::fma(dir(i), dist, p0(i));
      pDist.set(i, coord);
    }
    return Point<dim, fptype>(this->origin, pDist);
  }

  CUDA_CALLABLE Point<dim, fptype> getIntercept() const {
    return intercept;
  }

  CUDA_CALLABLE Vector<dim, fptype> getDirection() const {
    return dir;
  }

  CUDA_CALLABLE Line<dim, fptype> operator=(
      const Line<dim, fptype>::LineData &l) {
    this->origin = l.p.o;
    intercept = l.p;
    dir = l.d;
    return *this;
  }

  CUDA_CALLABLE bool operator==(
      const Line<dim, fptype> &l) const {
    return (this->origin == l.origin &&
            intercept == l.intercept && dir == l.dir);
  }

  CUDA_CALLABLE LineData copyData() const {
    LineData l;
    dir.copyData(l.d);
    intercept.copyData(l.p);
  }

  CUDA_CALLABLE LineData &copyData(LineData &l) const {
    dir.copyData(l.d);
    intercept.copyData(l.p);
    return l;
  }

#ifdef __CUDACC__
  std::shared_ptr<LineData> cudaCopy() const {
    LineData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    auto safeMem =
        std::shared_ptr<LineData>(cudaMem, cudaFree);
    err = cudaCopy(safeMem);
    return safeMem;
  }

  cudaError_t cudaCopy(
      std::shared_ptr<LineData> cudaMem) const {
    cudaError_t err = dir.cudaCopy(&cudaMem->d);
    err = this->intercept.cudaCopy(&cudaMem->p);
    return err;
  }

  cudaError_t cudaCopy(LineData *cudaMem) const {
    cudaError_t err = dir.cudaCopy(&cudaMem->d);
    err = this->intercept.cudaCopy(&cudaMem->p);
    return err;
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<LineData> cudaMem) {
    cudaError_t err =
        this->intercept.cudaRetrieve(&cudaMem->p);
    err = this->dir.cudaRetrieve(&cudaMem->d);
    return err;
  }

  cudaError_t cudaRetrieve(LineData *cudaMem) {
    cudaError_t err =
        this->intercept.cudaRetrieve(&cudaMem->p);
    err = this->dir.cudaRetrieve(&cudaMem->d);
    return err;
  }
#endif

  friend std::ostream &operator<<(
      std::ostream &os, const Line<dim, fptype> &l) {
    os << l.dir << " * t + " << l.intercept;
    return os;
  }

  template <int, typename>
  friend class Line;

 protected:
  Point<dim, fptype> intercept;
  Vector<dim, fptype> dir;
};
};

#endif
