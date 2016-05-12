
#include "accurate_intersections.hpp"

namespace Geometry {

template <int dim, typename fptype>
class IntersectionResultantMulti
    : public IntersectionBase<
          dim, fptype,
          IntersectionResultantMulti<dim, fptype>> {
 public:
  static constexpr const int numPartialProds = 6;
  mutable Array<mpfr::mpreal, numPartialProds> partialProds;
  mutable mpfr::mpreal scaledDisc;
  mutable mpfr::mpreal center;

  IntersectionResultantMulti() : IntersectionBase() {
    partialProds[0] = mpfr::mpreal(NAN);
  }

  IntersectionResultantMulti(
      const Quadric<dim, fptype> &quad,
      const Line<dim, fptype> &line, fptype intPos = NAN,
      fptype otherIntPos = NAN,
      fptype absErrMargin = defAbsPrecision)
      : IntersectionBase(quad, line, intPos, otherIntPos,
                         absErrMargin, 0) {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  template <typename cmpAlg>
  IntersectionResultantMulti(
      const IntersectionBase<dim, fptype, cmpAlg> &i)
      : IntersectionBase(i) {
    partialProds[numPartialProds - 1] = mpfr::mpreal(NAN);
  }

  IntersectionResultantMulti<dim, fptype> operator=(
      const IntersectionResultantMulti<dim, fptype> &i) {
    IntersectionBase<dim, fptype,
                     IntersectionResultantMulti<
                         dim, fptype>>::operator=(i);
    IntersectionBase<dim, fptype,
                     IntersectionResultant<dim, fptype>>::
    operator=(i);
    center = i.center;
    scaledDisc = i.scaledDisc;
    partialProds = i.partialProds;
    return *this;
  }

  bool ppReady() const {
    return !MathFuncs::MathFuncs<mpfr::mpreal>::isnan(
        partialProds[numPartialProds - 1]);
  }

  const mpfr::mpreal &getCenter() const {
    if(!ppReady()) {
      getPartialProd(0);
    }
    return center;
  }

  const mpfr::mpreal &getDiscriminant() const {
    if(!ppReady()) {
      getPartialProd(1);
    }
    return scaledDisc;
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
      Quadric<dim, mpfr::mpreal> quad(q);
      Line<dim, mpfr::mpreal> line(l);
      Polynomial<2, mpfr::mpreal> lqInt =
          quad.calcLineDistPoly(line);
      center =
          mpfr::div(lqInt.get(1), lqInt.get(2), partPrec);
      scaledDisc = mpfr::sub(
          mpfr::sqr(center, discPrec),
          mpfr::div(lqInt.get(0), lqInt.get(2), partPrec),
          discPrec);
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
      const IntersectionBase<dim, fptype, true> &i) const {
    assert(l == i.l);
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

  fptype accurateCompare_EqualDiscs(
      const IntersectionBase<dim, fptype, true> &i) const {
    /* For equal discriminants,
     * we can just compare the discriminants
     * to the difference in centers
     */
    const mpfr::mpreal centralDiff =
        getCenter() - i.getCenter();
    bool rootSign1 = MathFuncs::MathFuncs<fptype>::signbit(
        intPos - otherIntPos);
    bool rootSign2 = MathFuncs::MathFuncs<fptype>::signbit(
        i.intPos - i.otherIntPos);
    bool cdSign =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            centralDiff);
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
      if(mpfr::iszero(centralDiff)) {
        /* If the centers are on top of each other,
         * then the roots must be as well
         */
        return 0.0;
      }
      if(cdSign == 0) {
        return fptype(1.0);
      } else {
        return fptype(-1.0);
      }
    } else {
      /* The central difference squared divided by four
       * compared to the discriminant will give us the
       * result
       */
      mpfr::mpreal centralDiffSqr =
          mpfr::sqr(centralDiff, -1) >> 2;
      /* This must be computed to the full precision to
       * verify it's not zero
       */
      mpfr::mpreal cmpSign =
          mpfr::sub(getDiscriminant(), centralDiffSqr,
                    centralDiffSqr.getPrecision());
      if(mpfr::iszero(cmpSign)) {
        /* The central difference is twice the square root
         * of the discriminant, so the roots are on top of
         * each other
         */
        return fptype(0.0);
      }
      bool cmpDistSign =
          MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
              cmpSign);
      bool finalSign = rootSign1 ^ cmpDistSign;
      if(finalSign == 0) {
        return fptype(1.0);
      } else {
        return fptype(-1.0);
      }
    }
    return fptype(NAN);
  }

  fptype accurateCompare_OneRepeated(
      const IntersectionBase<dim, fptype, true> &i) const {
    /* We just need to compare the squared difference of
     * centers divided by four to the other discriminant
     */
    bool rootSign = MathFuncs::MathFuncs<fptype>::signbit(
        i.intPos - i.otherIntPos);
    mpfr::mpreal centralDiff = getCenter() - i.getCenter();
    bool cdSign =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            centralDiff);
    /* If they're not on the same side of the center, the
     * result must be the sign of the central difference
     */
    if(cdSign != rootSign) {
      if(cdSign == 1 && rootSign == 0) {
        return fptype(-1.0);
      } else {
        return fptype(1.0);
      }
    }
    mpfr::mpreal centraldiscdiff =
        mpfr::sub(mpfr::sqr(centralDiff, -1),
                  i.getDiscriminant(), -1);
    if(mpfr::iszero(centraldiscdiff)) {
      return fptype(0.0);
    }
    bool cmpSign =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            centraldiscdiff);
    if(cmpSign == rootSign) {
      return fptype(1.0);
    } else {
      return fptype(-1.0);
    }
  }

  fptype accurateCompare_TwoRepeated(
      const IntersectionBase<dim, fptype, true> &i) const {
    mpfr::mpreal centralDiff =
        mpfr::sub(getCenter(), i.getCenter(), -1);
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
      const IntersectionBase<dim, fptype, true> &i) const {
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
        intPos - i.otherIntPos, otherIntPos - i.intPos,
        otherIntPos - i.otherIntPos};
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
      const IntersectionBase<dim, fptype, true> &i) const {
    if(getDiscriminant() == i.getDiscriminant()) {
      return fptype(0.0);
    } else {
      /* The discriminant must be less than the other
       * intersections discriminant, so the root we are
       * comparing against determines the sign
       */
      if(i.intPos < i.otherIntPos) {
        return fptype(1.0);
      } else {
        return fptype(-1.0);
      }
    }
  }

  fptype accurateCompare_Multi(
      const IntersectionBase<dim, fptype, true> &i) const {
    mpfr::mpreal centralDiff = getCenter() - i.getCenter();
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
        intPos - otherIntPos);
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
    const unsigned cdSqrPrec =
        centralDiff.getPrecision() * 2;
    mpfr::mpreal cdSqr = mpfr::sqr(centralDiff, cdSqrPrec);
    int cdSqrCmpDim = cdSqr >= i.getDiscriminant();
    int rootSignDim = 2 * intNeg + otherIntNeg;
    return caseTable[centralDiffNeg][resultantDim]
                    [cdSqrCmpDim][rootSignDim];
  }

  fptype accurateCompare(
      const IntersectionBase<dim, fptype, true> &i) const {
    /* First determine which root is more likely */
    if(getDiscriminant() > i.getDiscriminant()) {
      return -i.accurateCompare(*this);
    }
    if(getDiscriminant() == i.getDiscriminant()) {
      fptype ret = accurateCompare_EqualDiscs(i);
      return ret;
    }
    if(mpfr::iszero(getDiscriminant())) {
      if(mpfr::iszero(i.getDiscriminant())) {
        fptype ret = accurateCompare_TwoRepeated(i);
        return ret;
      } else {
        fptype ret = accurateCompare_OneRepeated(i);
        return ret;
      }
    }
    if(isUndetermined(intPos, i.otherIntPos) ||
       isUndetermined(otherIntPos, i.intPos) ||
       isUndetermined(otherIntPos, i.otherIntPos)) {
      fptype ret = accurateCompare_Multi(i);
      return ret;
    } else {
      fptype ret = accurateCompare_One(i);
      return ret;
    }
  }
};
}
