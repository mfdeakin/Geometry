
#ifndef _ACCURATE_INTERSECTIONS_INCPREC_HPP_
#define _ACCURATE_INTERSECTIONS_INCPREC_HPP_

#include "accurate_intersections.hpp"
#include "mpreal.hpp"

namespace Geometry {

template <int dim, typename fptype>
class IntersectionIncreasedPrec
    : public IntersectionBase<
          dim, fptype,
          IntersectionIncreasedPrec<dim, fptype>> {
 public:
  mutable mpfr::mpreal mpIntPos;

  IntersectionIncreasedPrec()
      : IntersectionBase<
            dim, fptype,
            IntersectionIncreasedPrec<dim, fptype>>() {
    mpIntPos = mpfr::mpreal(NAN);
  }

  IntersectionIncreasedPrec(
      const Quadric<dim, fptype> &quad,
      const Line<dim, fptype> &line, fptype intPos = NAN,
      fptype otherIntPos = NAN,
      fptype absErrMargin = defAbsPrecision)
      : IntersectionBase<
            dim, fptype,
            IntersectionIncreasedPrec<dim, fptype>>(
            quad, line, intPos, otherIntPos, absErrMargin) {
    mpIntPos = mpfr::mpreal(NAN);
  }

  template <typename cmpAlg>
  IntersectionIncreasedPrec(
      const IntersectionBase<dim, fptype, cmpAlg> &i)
      : IntersectionBase<
            dim, fptype,
            IntersectionIncreasedPrec<dim, fptype>>(i) {
    mpIntPos = mpfr::mpreal(NAN);
  }

  IntersectionIncreasedPrec<dim, fptype> operator=(
      const IntersectionIncreasedPrec<dim, fptype> &i) {
    IntersectionBase<dim, fptype, IntersectionIncreasedPrec<
                                      dim, fptype>>::
    operator=(i);
    return *this;
  }

  bool rootReady() const {
    return !MathFuncs::MathFuncs<mpfr::mpreal>::isnan(
        mpIntPos);
  }

  mpfr::mpreal getRoot() const {
    if(!rootReady()) {
      const unsigned prevPrec =
          mpfr::mpreal::get_default_prec();
      constexpr const unsigned machPrec =
          GenericFP::fpconvert<fptype>::precision;
      mpfr::mpreal::set_default_prec(machPrec);
      Line<dim, mpfr::mpreal> line(this->l);
      Quadric<dim, mpfr::mpreal> quad(this->q);

      constexpr const unsigned coeffPrec = 3 * machPrec;
      mpfr::mpreal::set_default_prec(coeffPrec);
      auto poly = quad.calcLineDistPoly(line);

      constexpr const unsigned rootPrec = 8 * coeffPrec;
      mpfr::mpreal::set_default_prec(rootPrec);
      auto roots = poly.calcRoots();

      mpfr::mpreal::set_default_prec(prevPrec);
      if((this->intPos < this->otherIntPos &&
          roots[0] > roots[1]) ||
         (this->intPos > this->otherIntPos &&
          roots[0] < roots[1])) {
        mpIntPos = roots[1];
      } else {
        mpIntPos = roots[0];
      }
    }
    return mpIntPos;
  }

  fptype accurateCompare(
      const IntersectionIncreasedPrec<dim, fptype> &i)
      const {
    mpfr::mpreal delta = getRoot() - i.getRoot();
    if(mpfr::iszero(delta)) {
      return fptype(0.0);
    } else {
      int sign =
          MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
              delta);
      if(sign == 1) {
        return fptype(-1.0);
      } else {
        return fptype(1.0);
      }
    }
  }
};
}

#endif
