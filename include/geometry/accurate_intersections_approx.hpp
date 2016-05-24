
#ifndef _ACCURATE_INTERSECTIONS_APPROX_HPP_
#define _ACCURATE_INTERSECTIONS_APPROX_HPP_

#include "accurate_intersections.hpp"

namespace Geometry {

template <int dim, typename inptype, typename calctype>
class IntersectionApproximate
    : public IntersectionBase<
          dim, inptype,
          IntersectionApproximate<dim, inptype, calctype>> {
 public:
  IntersectionApproximate()
      : IntersectionBase<dim, inptype,
                         IntersectionApproximate<
                             dim, inptype, calctype>>() {
    calcIntPos = calctype(NAN);
  }

  IntersectionApproximate(
      const Quadric<dim, inptype> &quad,
      const Line<dim, inptype> &line, inptype intPos = NAN,
      inptype otherIntPos = NAN,
      inptype absErrMargin = defAbsPrecision)
      : IntersectionBase<dim, inptype,
                         IntersectionApproximate<
                             dim, inptype, calctype>>(
            quad, line, intPos, otherIntPos, absErrMargin) {
    calcIntPos = calctype(NAN);
  }

  template <typename cmpAlg>
  IntersectionApproximate(
      const IntersectionBase<dim, inptype, cmpAlg> &i)
      : IntersectionBase<dim, inptype,
                         IntersectionApproximate<
                             dim, inptype, calctype>>(i) {
    calcIntPos = calctype(NAN);
  }

  IntersectionApproximate<dim, inptype, calctype> operator=(
      const IntersectionApproximate<dim, inptype, calctype>
          &i) {
    IntersectionBase<
        dim, inptype,
        IntersectionApproximate<dim, inptype, calctype>>::
    operator=(i);
    calcIntPos = calctype(NAN);
    return *this;
  }

  bool rootReady() const {
    return !MathFuncs::MathFuncs<calctype>::isnan(
        calcIntPos);
  }

  calctype getRoot() const {
    if(!rootReady()) {
      Line<dim, calctype> line(this->l);
      Quadric<dim, calctype> quad(this->q);

      auto poly = quad.calcLineDistPoly(line);
      auto roots = poly.calcRoots();

      if((this->intPos < this->otherIntPos &&
          roots[0] > roots[1]) ||
         (this->intPos > this->otherIntPos &&
          roots[0] < roots[1])) {
        calcIntPos = roots[1];
      } else {
        calcIntPos = roots[0];
      }
    }
    return calcIntPos;
  }

  inptype accurateCompare(
      const IntersectionApproximate<dim, inptype, calctype>
          &i) const {
    return getRoot() - i.getRoot();
  }

 protected:
  mutable calctype calcIntPos;
};
}

#endif
