
#ifndef _ACCURATE_INTERSECTIONS_RESULTANT_HPP_
#define _ACCURATE_INTERSECTIONS_RESULTANT_HPP_

#include "accurate_intersections.hpp"

namespace Geometry {

template <int dim, typename fptype>
class IntersectionResultant
    : public IntersectionBase<
          dim, fptype, IntersectionResultant<dim, fptype>> {
 public:
  static constexpr const int numPartialProds = 6;
  mutable Array<mpfr::mpreal, numPartialProds> partialProds;

  IntersectionResultant()
      : IntersectionBase<dim, fptype, IntersectionResultant<
                                          dim, fptype>>() {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  IntersectionResultant(
      const Quadric<dim, fptype> &quad,
      const Line<dim, fptype> &line, fptype intPos = NAN,
      fptype otherIntPos = NAN,
      fptype absErrMargin = defAbsPrecision)
      : IntersectionBase<dim, fptype, IntersectionResultant<
                                          dim, fptype>>(
            quad, line, intPos, otherIntPos, absErrMargin) {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  template <typename cmpAlg>
  IntersectionResultant(
      const IntersectionBase<dim, fptype, cmpAlg> &i)
      : IntersectionBase<dim, fptype, IntersectionResultant<
                                          dim, fptype>>(i) {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  IntersectionResultant<dim, fptype> operator=(
      const IntersectionResultant<dim, fptype> &i) {
    IntersectionBase<dim, fptype,
                     IntersectionResultant<dim, fptype>>::
    operator=(i);
    partialProds = i.partialProds;
    return *this;
  }

  bool ppReady() const {
    return !MathFuncs::MathFuncs<mpfr::mpreal>::isnan(
        partialProds[numPartialProds - 1]);
  }

  const mpfr::mpreal &getDiscriminant() const {
    return getPartialProd(1);
  }

  const mpfr::mpreal &getPartialProd(int idx) const {
    if(!ppReady()) {
      const unsigned prevPrec =
          mpfr::mpreal::get_default_prec();
      /* This is less efficient than we'd like,
       * but is the only way to implement this in a
       * manner which is agnostic floating point type
       */
      constexpr const unsigned machPrec =
          GenericFP::fpconvert<fptype>::precision;
      constexpr const unsigned coeffPrec = 3 * machPrec;
      constexpr const int partPrec = 2 * coeffPrec;
      constexpr const int discPrec = 2 * partPrec;
      mpfr::mpreal::set_default_prec(coeffPrec);
      Quadric<dim, mpfr::mpreal> quad(this->q);
      Line<dim, mpfr::mpreal> line(this->l);
      Polynomial<2, mpfr::mpreal> lqInt =
          quad.calcLineDistPoly(line);
      partialProds[2] =
          mpfr::mult(lqInt.get(2), lqInt.get(2), partPrec);
      partialProds[3] =
          mpfr::mult(lqInt.get(0), lqInt.get(1), partPrec);
      partialProds[4] =
          mpfr::mult(lqInt.get(0), lqInt.get(2), partPrec);
      /* This is also the discriminant */
      partialProds[1] =
          mpfr::sub(mpfr::sqr(lqInt.get(1), partPrec),
                    partialProds[4], partPrec);
      partialProds[0] =
          mpfr::mult(lqInt.get(0), lqInt.get(0), partPrec);
      partialProds[5] =
          mpfr::mult(lqInt.get(1), lqInt.get(2), partPrec);
      mpfr::mpreal::set_default_prec(prevPrec);
    }
    return partialProds[idx];
  }

  mpfr::mpreal resultantDet(
      const IntersectionResultant<dim, fptype> &i) const {
    assert(this->l == i.l);
    /* The fastest way to compute our determinant
     * 1.0(0 0 5 5)+1.0(2 2 3 3)+1.0(1 1 3 5)+1.0(4 4 0 2)+
     * -1.0(1 2 3 4)-1.0(0 1 4 5)-2.0(0 2 3 5)
     * 35 FLOPs
     * (0 5)((0 5)-2.0(2 3)-(1 4))+(2 3)((2 3)-(1 4))+
     * (1 1 3 5)+(4 4 0 2)
     * 18 FLOPs, which is nice,
     * but we can do better with caching
     *
     * ((4 4)-(3 5))(0 2)+((1 1)-(0 2))(3 5)
     * Amortize like so:
     * (0 0)(5 5)+(2 2)(3 3)+((1 1)-(0 2))(3 5)+
     * ((4 4)-(3 5))(0 2)-(1 2)(3 4)-(0 1)(4 5)
     * where all partial products (0 0), (2 2), etc. are
     * computed only once.
     * This brings us to 11 FLOPs per comparison,
     * plus some small constant for the initial computation
     */
    const unsigned prevPrec =
        mpfr::mpreal::get_default_prec();
    constexpr const unsigned machPrec =
        GenericFP::fpconvert<fptype>::precision;
    constexpr const unsigned coeffPrec = 3 * machPrec;
    constexpr const unsigned partPrec = 2 * coeffPrec;
    constexpr const unsigned detPrec = 2 * partPrec;
    /* Now compute the larger terms */
    mpfr::mpreal detTerms[] = {
        //(0 0)(5 5)
        mpfr::mult(getPartialProd(0), i.getPartialProd(2),
                   detPrec),
        //(2 2)(3 3)
        mpfr::mult(getPartialProd(2), i.getPartialProd(0),
                   detPrec),
        //((1 1)-(0 2))(3 5)
        mpfr::mult(getPartialProd(1), i.getPartialProd(4),
                   detPrec),
        //(0 2)((4 4)-(3 5))
        mpfr::mult(getPartialProd(4), i.getPartialProd(1),
                   detPrec),
        //-(1 2)(3 4)
        -mpfr::mult(getPartialProd(5), i.getPartialProd(3),
                    detPrec),
        //-(0 1)(4 5)
        -mpfr::mult(getPartialProd(3), i.getPartialProd(5),
                    detPrec)};
    constexpr const int numTerms =
        sizeof(detTerms) / sizeof(detTerms[0]);
    /* There are only 4 terms to sum,
     * which isn't enough for compensated summation to be
     * worthwhile in my experience
     */
    mpfr::mpreal det(detTerms[0]);
    for(int i = 1; i < numTerms; i++) {
      det += detTerms[i];
    }
    mpfr::mpreal::set_default_prec(prevPrec);
    return det;
  }

  fptype accurateCompare(
      const IntersectionBase<
          dim, fptype, IntersectionResultant<dim, fptype>>
          &i) const {
    /* Since this is the only undetermined root,
     * the resultant will determine the sign
     */
    const mpfr::mpreal resultant = resultantDet(i);
    if(mpfr::iszero(resultant)) {
      return fptype(0.0);
    }
    /* So just compute the signs of the resultant terms.
     * Since the sign flips with each negative, just
     * use XOR on one bit to keep track of the final sign
     */
    constexpr const int numPartialTerms = 3;
    fptype partResTerms[numPartialTerms] = {
        this->intPos - i.otherIntPos,
        this->otherIntPos - i.intPos,
        this->otherIntPos - i.otherIntPos};
    int numNeg = resultant < 0;
    for(int j = 0; j < numPartialTerms; j++) {
      if(partResTerms[j] < 0.0)
        numNeg ^= 1;
      else if(partResTerms[j] == 0.0)
        /* Zero, so the sign of this difference can't be
         * determined with this method
         */
        return 0.0;
    }
    /* The sign of the product of the resultant and the
     * known differences of roots has the sign of the
     * final difference of the root
     */
    if(numNeg == 0) {
      return 1.0;
    } else {
      return -1.0;
    }
  }
};
}

#endif
