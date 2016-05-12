
#ifndef _ACCURATE_INTERSECTIONS_APPROX_HPP_
#define _ACCURATE_INTERSECTIONS_APPROX_HPP_

#include "accurate_intersections.hpp"

namespace Geometry {

template <int dim, typename fptype>
class IntersectionApproximate
    : public IntersectionBase<
          dim, fptype,
          IntersectionApproximate<dim, fptype>> {
 public:
  IntersectionApproximate()
      : IntersectionBase<
            dim, fptype,
            IntersectionApproximate<dim, fptype>>() {}

  IntersectionApproximate(
      const Quadric<dim, fptype> &quad,
      const Line<dim, fptype> &line, fptype intPos = NAN,
      fptype otherIntPos = NAN,
      fptype absErrMargin = defAbsPrecision)
      : IntersectionBase<
            dim, fptype,
            IntersectionApproximate<dim, fptype>>(
            quad, line, intPos, otherIntPos, absErrMargin) {
  }

  template <typename cmpAlg>
  IntersectionApproximate(
      const IntersectionBase<dim, fptype, cmpAlg> &i)
      : IntersectionBase<
            dim, fptype,
            IntersectionApproximate<dim, fptype>>(i) {}

  IntersectionApproximate<dim, fptype> operator=(
      const IntersectionApproximate<dim, fptype> &i) {
    IntersectionBase<dim, fptype,
                     IntersectionApproximate<dim, fptype>>::
    operator=(i);
    return *this;
  }

  fptype accurateCompare(
      const IntersectionBase<
          dim, fptype, IntersectionApproximate<dim, fptype>>
          &i) const {
    return this->intPos - i.intPos;
  }
};
}

#endif
