
#ifndef _ACCURATE_INTERSECTIONS_MULTI_HPP_
#define _ACCURATE_INTERSECTIONS_MULTI_HPP_

#include "accurate_intersections.hpp"

namespace Geometry {

template <int dim, typename fptype>
class IntersectionResultantMulti
    : public IntersectionBase<
          dim, fptype,
          IntersectionResultantMulti<dim, fptype>> {
 public:
  static constexpr const int numPartialProds = 6;
  mutable Polynomial<2, mpfr::mpreal> lqInt;
  mutable Array<mpfr::mpreal, numPartialProds> partialProds;

  IntersectionResultantMulti()
      : IntersectionBase<
            dim, fptype,
            IntersectionResultantMulti<dim, fptype>>() {
    partialProds[0] = mpfr::mpreal(NAN);
  }

  IntersectionResultantMulti(
      const Quadric<dim, fptype> &quad,
      const Line<dim, fptype> &line, fptype intPos = NAN,
      fptype otherIntPos = NAN,
      fptype absErrMargin = defAbsPrecision)
      : IntersectionBase<
            dim, fptype,
            IntersectionResultantMulti<dim, fptype>>(
            quad, line, intPos, otherIntPos, absErrMargin) {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  template <typename cmpAlg>
  IntersectionResultantMulti(
      const IntersectionBase<dim, fptype, cmpAlg> &i)
      : IntersectionBase<
            dim, fptype,
            IntersectionResultantMulti<dim, fptype>>(i) {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  IntersectionResultantMulti<dim, fptype> operator=(
      const IntersectionResultantMulti<dim, fptype> &i) {
    IntersectionBase<dim, fptype,
                     IntersectionResultantMulti<
                         dim, fptype>>::operator=(i);
    partialProds = i.partialProds;
    return *this;
  }

  bool intermediatesReady() const {
    return !MathFuncs::MathFuncs<mpfr::mpreal>::isnan(
        partialProds[numPartialProds - 1]);
  }

  static constexpr unsigned coeffPrec() {
    return GenericFP::fpconvert<fptype>::precision * 3;
  }

  static constexpr unsigned partialProdPrec() {
    return coeffPrec() * 2;
  }

  void cacheIntermediates() const {
    const unsigned prevPrec =
        mpfr::mpreal::get_default_prec();
    /* This is less efficient than we'd like,
     * but is the only way to implement this in a
     * manner which is agnostic floating point type
     */
    mpfr::mpreal::set_default_prec(coeffPrec());
    Quadric<dim, mpfr::mpreal> quad(this->q);
    Line<dim, mpfr::mpreal> line(this->l);
    lqInt = quad.calcLineDistPoly(line);
    partialProds[2] = mpfr::mult(lqInt.get(2), lqInt.get(2),
                                 partialProdPrec());
    partialProds[3] = mpfr::mult(lqInt.get(0), lqInt.get(1),
                                 partialProdPrec());
    partialProds[4] = mpfr::mult(lqInt.get(0), lqInt.get(2),
                                 partialProdPrec());
    /* This is also the discriminant */
    partialProds[1] = mpfr::sub(
        mpfr::sqr(lqInt.get(1), partialProdPrec()),
        partialProds[4], partialProdPrec());
    partialProds[0] = mpfr::mult(lqInt.get(0), lqInt.get(0),
                                 partialProdPrec());
    partialProds[5] = mpfr::mult(lqInt.get(1), lqInt.get(2),
                                 partialProdPrec());
    mpfr::mpreal::set_default_prec(prevPrec);
  }

  fptype discDiffSign(
      const IntersectionResultantMulti &i) const {
    Polynomial<2, mpfr::mpreal> p1 = getLQIntPoly(),
                                p2 = i.getLQIntPoly();
    mpfr::mpreal scaledDisc1 = mpfr::mult(
                     getPartialProd(1),
                     mpfr::sqr(p2.get(2),
                               partialProdPrec()),
                     4 * coeffPrec()),
                 scaledDisc2 = mpfr::mult(
                     i.getPartialProd(1),
                     mpfr::sqr(p1.get(2),
                               partialProdPrec()),
                     4 * coeffPrec());
  }

  bool discZero() const {
    mpfr::mpreal disc = getPartialProd(1);
    return MathFuncs::MathFuncs<mpfr::mpreal>::iszero(disc);
  }

  const mpfr::mpreal &getPartialProd(int idx) const {
    if(!intermediatesReady()) {
      cacheIntermediates();
    }
    return partialProds[idx];
  }

  const Polynomial<2, mpfr::mpreal> &getLQIntPoly() const {
    if(!intermediatesReady()) {
      cacheIntermediates();
    }
    return lqInt;
  }

  const fptype centralDiffSign(
      const IntersectionResultantMulti &i) const {
    /* The sign of the central differences can be computed
     * exactly as
     * sign(b_2 a_1 - b_1 a_2) sign(a_1) sign(a_2)
     */
    constexpr const int termPrec = 2 * coeffPrec();
    const Polynomial<2, mpfr::mpreal> p1 = getLQIntPoly(),
                                      p2 = i.getLQIntPoly();
    mpfr::mpreal sign1 =
        mpfr::mult(p1.get(1), p2.get(2), termPrec) -
        mpfr::mult(p1.get(2), p2.get(1), termPrec);
    if(MathFuncs::MathFuncs<mpfr::mpreal>::iszero(sign1)) {
      return fptype(0.0);
    }
    int negCounter =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(sign1);
    negCounter ^=
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            p1.get(2));
    negCounter ^=
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            p2.get(2));
    if(negCounter) {
      return fptype(-1.0);
    }
    return fptype(1.0);
  }

  fptype evalOtherCenterSign(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    constexpr const unsigned termPrec = 3 * coeffPrec();
    const Polynomial<2, mpfr::mpreal> p1 = getLQIntPoly(),
                                      p2 = i.getLQIntPoly();
    mpfr::mpreal terms[] = {
        mpfr::mult(p1.get(2),
                   mpfr::sqr(p2.get(1), 2 * coeffPrec()),
                   termPrec),
        mpfr::mult(-p1.get(1),
                   mpfr::mult(p2.get(1), p2.get(2),
                              2 * coeffPrec()),
                   termPrec)
            << 1,
        mpfr::mult(p1.get(0),
                   mpfr::sqr(p2.get(2), 2 * coeffPrec()),
                   termPrec)};
    mpfr::mpreal result = terms[0] + terms[1] + terms[2];
    if(MathFuncs::MathFuncs<mpfr::mpreal>::iszero(result)) {
      return fptype(0.0);
    } else if(MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
                  result) == 0) {
      return fptype(1.0);
    } else {
      return fptype(-1.0);
    }
  }

  mpfr::mpreal resultantDet(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
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
    constexpr const unsigned detPrec =
        2 * partialProdPrec();
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

  fptype accurateCompare_EqualDiscs(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    /* For equal discriminants,
     * we can just compare the discriminants
     * to the difference in centers
     */
    bool rootSign1 = MathFuncs::MathFuncs<fptype>::signbit(
        this->intPos - this->otherIntPos);
    bool rootSign2 = MathFuncs::MathFuncs<fptype>::signbit(
        i.intPos - i.otherIntPos);
    fptype cdSign = centralDiffSign(i);
    /* For the case of ++,
     * the sign of the central difference gives the result
     * For the case of --,
     * the negative sign of the central difference gives the
     * result
     * For the cases of -+, +-,
     * comparing the squared central difference to the sum
     * of the discriminants gives the result
     */
    if(rootSign1 == rootSign2) {
      if(MathFuncs::MathFuncs<fptype>::iszero(cdSign)) {
        /* If the centers are on top of each other,
         * then the roots must be as well
         */
        return 0.0;
      }
      if(cdSign > 0.0) {
        return fptype(1.0);
      } else {
        return fptype(-1.0);
      }
    } else {
      /* The central difference squared divided by four
       * compared to the discriminant will give us the
       * result
       */
      Polynomial<2, mpfr::mpreal> p1 = getLQIntPoly(),
                                  p2 = i.getLQIntPoly();
      constexpr const int intermediatePrec =
          partialProdPrec() + 5;
      constexpr const int maxPrec =
          i.partialProdPrec() + partialProdPrec() + 5;
      mpfr::mpreal terms[] = {
          /* -2 a1 b1 a2 b2 */
          mpfr::mult(-getPartialProd(5),
                     i.getPartialProd(5) << 2,
                     2 * partialProdPrec()),
          /* a1^2 b2^2 */
          mpfr::mult(
              getPartialProd(2),
              mpfr::sqr(p2.get(1), 2 * i.coeffPrec()),
              2 * i.coeffPrec() + partialProdPrec()),
          /* a2^2 (-3 b1^2 + 16 a1 c1) */
          mpfr::mult(
              i.getPartialProd(2),
              mpfr::sum(mpfr::mult(-3.0, getPartialProd(1),
                                   partialProdPrec() + 2),
                        mpfr::mult(19.0, getPartialProd(4),
                                   intermediatePrec),
                        intermediatePrec),
              maxPrec)};
      constexpr const int lastTerm =
          sizeof(terms) / sizeof(terms[0]) - 1;
      mpfr::mpreal sum(terms[lastTerm]);
      for(int i = lastTerm - 1; i >= 0; i--) {
        sum += terms[i];
      }
      return fptype(sum);
    }
  }

  fptype accurateCompare_OneRepeated(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    /* We just need to compare the squared difference of
     * centers divided by four to the other discriminant
     */
    bool rootSign = MathFuncs::MathFuncs<fptype>::signbit(
        i.intPos - i.otherIntPos);
    fptype centralDiff = centralDiffSign(i);
    bool cdSign =
        MathFuncs::MathFuncs<fptype>::signbit(centralDiff);
    /* If they're not on the same side of the center, the
     * result must be the sign of the central difference
     */
    if(cdSign != rootSign) {
      if(cdSign == 1 && rootSign == 0) {
        return fptype(-1.0);
      } else {
        return fptype(1.0);
      }
    } else {
    }
    return fptype(NAN);
  }

  fptype accurateCompare_TwoRepeated(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    fptype centralDiff = centralDiffSign(i);
    if(mpfr::iszero(centralDiff)) {
      return fptype(0.0);
    } else {
      bool signbit =
          MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
              centralDiff);
      if(signbit) {
        return fptype(-1.0);
      } else {
        return fptype(1.0);
      }
    }
  }

  fptype accurateCompare_One(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
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
    for(int i = 0; i < numPartialTerms; i++) {
      if(partResTerms[i] < 0.0)
        numNeg ^= 1;
      else if(partResTerms[i] == 0.0)
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

  fptype accurateCompare_EqCenter(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    return fptype(NAN);
  }

  fptype accurateCompare_Multi(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    fptype centralDiff = centralDiffSign(i);
    if(mpfr::iszero(centralDiff)) {
      return accurateCompare_EqCenter(i);
    }
    /* Now use the central difference and the sides being
     * used for the roots to determine the case we're
     * solving
     */
    int centralDiffNeg =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            centralDiff);
    int intNeg = MathFuncs::MathFuncs<fptype>::signbit(
        this->intPos - this->otherIntPos);
    int otherIntNeg = MathFuncs::MathFuncs<fptype>::signbit(
        i.intPos - i.otherIntPos);
    if(centralDiffNeg != intNeg &&
       centralDiffNeg != otherIntNeg) {
      /* This case is obvious with a picture */
      if(centralDiffNeg) {
        return fptype(-1.0);
      } else {
        return fptype(1.0);
      }
    }

    /* First dimension is the central difference sign (+,-);
     * Second dimension is the resultant <0, =0, >0;
     * Third dimension is whether the central difference
     * squared is less than the larger discriminant (<, >=);
     * Fourth dimension is the choice of root pairs
     * from ++, +-, -+, --
     */
    constexpr const fptype caseTable[2][3][2][4] = {
        {/* central difference sign is positive */
         {
             /* resultant >0 */
             {/* central difference squared is less than
                 the larger discriminant */
              fptype(-1.0), fptype(1.0), fptype(-1.0),
              fptype(1.0)},
             {/* central difference squared is at least
               * the larger discriminant */
              fptype(1.0), fptype(1.0), fptype(1.0),
              fptype(1.0)},
         },
         {
             /* resultant =0 */
             {/* central difference squared is less than
                 the larger discriminant */
              fptype(0.0), fptype(1.0), fptype(-1.0),
              fptype(1.0)},
             {/* central difference squared is at least
               * the larger discriminant */
              fptype(1.0), fptype(1.0), fptype(0.0),
              fptype(1.0)},
         },
         {
             /* resultant <0 */
             {/* central difference squared is less than
                 the larger discriminant */
              fptype(1.0), fptype(1.0), fptype(-1.0),
              fptype(-1.0)},
             {/* central difference squared is at least
               * the larger discriminant */
              fptype(1.0), fptype(1.0), fptype(-1.0),
              fptype(-1.0)},
         }},
        {/* central difference sign is negative */
         {
             /* resultant >0 */
             {/* central difference squared is less than
               * the larger discriminant */
              fptype(-1.0), fptype(1.0), fptype(-1.0),
              fptype(1.0)},
             {/* central difference squared is at least
               * the larger discriminant */
              fptype(-1.0), fptype(-1.0), fptype(-1.0),
              fptype(-1.0)},
         },
         {
             /* resultant =0 */
             {/* central difference squared is less than
               * the larger discriminant */
              fptype(-1.0), fptype(1.0), fptype(-1.0),
              fptype(0.0)},
             {/* central difference squared is at least
               * the larger discriminant */
              fptype(-1.0), fptype(0.0), fptype(-1.0),
              fptype(-1.0)},
         },
         {
             /* resultant <0 */
             {/* central difference squared is less than
                 the larger discriminant */
              fptype(-1.0), fptype(1.0), fptype(-1.0),
              fptype(-1.0)},
             {/* central difference squared is at least
               * the larger discriminant */
              fptype(-1.0), fptype(1.0), fptype(-1.0),
              fptype(-1.0)},
         }}};
    const mpfr::mpreal resultant = resultantDet(i);
    int resultantDim =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            resultant) *
        2;
    if(mpfr::iszero(resultant)) {
      resultantDim = 1;
    }
    return fptype(NAN);
  }

  fptype accurateCompare(
      const IntersectionResultantMulti<dim, fptype> &i)
      const {
    fptype ddSign = discDiffSign(i);
    if(MathFuncs::MathFuncs<fptype>::iszero(ddSign)) {
      fptype ret = accurateCompare_EqualDiscs(i);
      return ret;
    } else if(ddSign > 0.0) {
      return -i.accurateCompare(*this);
    } else if(discZero()) {
      if(i.discZero()) {
        fptype ret = accurateCompare_TwoRepeated(i);
        return ret;
      } else {
        fptype ret = accurateCompare_OneRepeated(i);
        return ret;
      }
    }
    if(this->isUndetermined(this->intPos, i.otherIntPos) ||
       this->isUndetermined(this->otherIntPos, i.intPos) ||
       this->isUndetermined(this->otherIntPos,
                            i.otherIntPos)) {
      fptype ret = accurateCompare_Multi(i);
      return ret;
    } else {
      fptype ret = accurateCompare_One(i);
      return ret;
    }
  }
};
}

#endif
