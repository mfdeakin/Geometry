
#include "accurate_intersections.hpp"
#include "mpreal.hpp"

namespace Geometry {

template <int dim, typename fptype>
class IntersectionIncreasedPrec
    : IntersectionBase<
          dim, fptype,
          IntersectionIncreasedPrec<dim, fptype>> {
 public:
  mutable Geometry::Array<mpfr::mpreal, 2> roots;

  IntersectionIncreasedPrec() : IntersectionBase() {}

  IntersectionIncreasedPrec(
      const Quadric<dim, fptype> &quad,
      const Line<dim, fptype> &line, fptype intPos = NAN,
      fptype otherIntPos = NAN,
      fptype absErrMargin = defAbsPrecision)
      : IntersectionBase(quad, line, intPos, otherIntPos,
                         absErrMargin, 0) {
    roots[1] = mpfr::mpreal(NAN);
  }

  template <typename cmpAlg>
  IntersectionIncreasedPrec(
      const IntersectionBase<dim, fptype, cmpAlg> &i)
      : IntersectionBase(i) {
    roots[1] = mpfr::mpreal(NAN);
  }

  IntersectionIncreasedPrec<dim, fptype> operator=(
      const IntersectionIncreasedPrec<dim, fptype> &i) {
    IntersectionBase<dim, fptype,
                     IntersectionApproximate<dim, fptype>>::
    operator=(i);
    return *this;
  }

  mpfr::mpreal getRoot(int rootNum) {}

  fptype accurateCompare(
      const IntersectionBase<
          dim, fptype, IntersectionApproximate<dim, fptype>>
          &i) const {
    const unsigned prevPrec =
        mpfr::mpreal::get_default_prec();
    constexpr const unsigned machPrec =
        GenericFP::fpconvert<fptype>::precision;
    mpfr::mpreal::set_default_prec(machPrec);

    constexpr const unsigned coeffPrec = 3 * machPrec;
    mpfr::mpreal::set_default_prec(coeffPrec);

    mpfr::mpreal::set_default_prec(prevPrec);
    return intPos - i.intPos;
  }
};
}
