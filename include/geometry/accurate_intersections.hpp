
#include "geometry.hpp"
#include "quadric.hpp"
#include "line.hpp"
#include "polynomial.hpp"
#include "genericfp.hpp"

#include <cmath>
#include <list>
#include <memory>
#include <utility>
#include <iostream>

#include "mpreal.hpp"

namespace Geometry {

template <int dim, typename fptype, bool isMachPrec>
class IntersectionBase;

template <int dim, typename fptype>
class IntersectionBase<dim, fptype, true> {
 public:
  Quadric<dim, fptype> q;
  Line<dim, fptype> l;
  fptype intPos, otherIntPos;
  fptype absErrMargin;
  mutable int numIP;

  static constexpr const int numPartialProds = 6;
  mutable Array<mpfr::mpreal, numPartialProds> partialProds;
  mpfr::mpreal center;

  IntersectionBase()
      : q(),
        intPos(NAN),
        otherIntPos(NAN),
        absErrMargin(defAbsPrecision),
        numIP(0) {
    partialProds[0] = mpfr::mpreal(NAN);
  }

  IntersectionBase(const Quadric<dim, fptype> &quad,
                   const Line<dim, fptype> &line,
                   fptype intPos = NAN,
                   fptype otherIntPos = NAN,
                   fptype absErrMargin = defAbsPrecision)
      : q(quad),
        l(line),
        intPos(intPos),
        otherIntPos(otherIntPos),
        absErrMargin(absErrMargin),
        numIP(0) {
    partialProds[0] = mpfr::mpreal(NAN);
  }

  IntersectionBase(
      const IntersectionBase<dim, fptype, true> &i)
      : q(i.q),
        l(i.l),
        intPos(i.intPos),
        otherIntPos(i.otherIntPos),
        absErrMargin(i.absErrMargin),
        numIP(i.numIP) {
    partialProds[0] = mpfr::mpreal(NAN);
  }

  IntersectionBase<dim, fptype, true> operator=(
      const IntersectionBase<dim, fptype, true> &i) {
    q = i.q;
    l = i.l;
    intPos = i.intPos;
    otherIntPos = i.otherIntPos;
    absErrMargin = i.absErrMargin;
    partialProds = i.partialProds;
    numIP = i.numIP;
    return *this;
  }

  bool ppReady() const {
    return !MathFuncs::MathFuncs<mpfr::mpreal>::isnan(
        partialProds[0]);
  }

  const mpfr::mpreal &getCenter() const {
    if(MathFuncs::MathFuncs<mpfr::mpreal>::isnan(center)) {
      center =
          mpfr::div(lqInt.get(1), lqInt.get(2), partPrec);
    }
    return center;
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
      mpfr::mpreal::set_default_prec(coeffPrec);
      Quadric<dim, mpfr::mpreal> quad(q);
      Line<dim, mpfr::mpreal> line(l);
      Polynomial<2, mpfr::mpreal> lqInt =
          quad.calcLineDistPoly(line);
      partialProds[2] =
          mpfr::mult(lqInt.get(2), lqInt.get(2), partPrec);
      partialProds[3] =
          mpfr::mult(lqInt.get(0), lqInt.get(1), partPrec);
      partialProds[4] =
          mpfr::mult(lqInt.get(0), lqInt.get(2), partPrec);
      partialProds[5] =
          mpfr::mult(lqInt.get(1), lqInt.get(2), partPrec);
      /* This is also the discriminant */
      partialProds[1] = mpfr::sub(
          mpfr::mult(lqInt.get(1), lqInt.get(1), partPrec),
          partialProds[4], partPrec);
      partialProds[0] =
          mpfr::mult(lqInt.get(0), lqInt.get(0), partPrec);
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
        /* If the ccenters are on top of each other,
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
          mpfr::sub(getPartialProd(1), centralDiffSqr,
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
    }
    mpfr::mpreal centraldiscdiff =
        mpfr::sub(mpfr::sqr(centralDiff, -1),
                  i.getPartialProd(1), -1);
    if(mpfr::iszero(centraldiscdiff)) {
      return fptype(0.0);
    }
    bool cmpSign =
        MathFuncs::MathFuncs<mpfr::mpreal>::signbit(
            centraldiscdiff);
    if(cmpSign == rootSign) {
      return fptype(1.0);
    }
    else {
      return fptype(-1.0);
    }
  }

  fptype accurateCompare(
      const IntersectionBase<dim, fptype, true> &i) const {
    /* Keep track of the number of precision increases
     * required for statistics
     */
    numIP++;
    /* First determine which root is more likely */
    if(getPartialProd(1) > i.getPartialProd(1)) {
      return -i.accurateCompare(*this);
    }
    if(getPartialProd(1) == i.getPartialProd(1)) {
      return accurateCompare_EqualDiscs(i);
    }
    if(mpfr::iszero(getPartialProd(1))) {
      if(mpfr::iszero(i.getPartialProd(1))) {
        return accurateCompare_TwoRepeated(i);
      }
      else {
        return accurateCompare_OneRepeated(i);
      }
    }
    mpfr::mpreal det = resultantDet(i);
    /* Now compute the signs of the other resultant terms.
     * Since the sign flips with each negative, just
     * use XOR on one bit to keep track of the final sign
     */
    constexpr const int numPartialTerms = 3;
    fptype partResTerms[numPartialTerms] = {
        intPos - i.otherIntPos, otherIntPos - i.intPos,
        otherIntPos - i.otherIntPos};
    int numNeg = det < 0;
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

  fptype compare(
      const IntersectionBase<dim, fptype, true> &i) const {
    fptype delta = intPos - i.intPos;
    if(MathFuncs::MathFuncs<fptype>::fabs(delta) <
           absErrMargin &&
       i.q != q &&
       !MathFuncs::MathFuncs<fptype>::isnan(otherIntPos) &&
       !MathFuncs::MathFuncs<fptype>::isnan(
           i.otherIntPos)) {
      fptype cmp = accurateCompare(i);
      if(cmp == 0.0)
        return delta;
      else
        return cmp;
    } else
      return delta;
  }

  int incPrecCount() const { return numIP; }

  static bool cmp(
      const IntersectionBase<dim, fptype, true> &lhs,
      const IntersectionBase<dim, fptype, true> &rhs) {
    fptype diff = lhs.compare(rhs);
    if(diff < 0) return true;
    return false;
  }
};

template <int dim, typename fptype>
class IntersectionBase<dim, fptype, false> {
 public:
  Quadric<dim, fptype> q;
  Line<dim, fptype> l;
  fptype intPos, otherIntPos;
  fptype absErrMargin;
  mutable int numIP;

  IntersectionBase() {}

  IntersectionBase(const Quadric<dim, fptype> &quad,
                   const Line<dim, fptype> &line,
                   fptype intPos = NAN,
                   fptype otherIntPos = NAN,
                   fptype absErrMargin = defAbsPrecision)
      : q(quad),
        l(line),
        intPos(intPos),
        otherIntPos(otherIntPos),
        absErrMargin(absErrMargin),
        numIP(0) {}

  IntersectionBase(
      const IntersectionBase<dim, fptype, false> &i)
      : q(i.q),
        l(i.l),
        intPos(i.intPos),
        otherIntPos(i.otherIntPos),
        absErrMargin(i.absErrMargin),
        numIP(i.numIP) {}

  IntersectionBase<dim, fptype, false> operator=(
      const IntersectionBase<dim, fptype, false> &i) {
    q = i.q;
    l = i.l;
    intPos = i.intPos;
    otherIntPos = i.otherIntPos;
    absErrMargin = i.absErrMargin;
    numIP = i.numIP;
    return *this;
  }

  fptype compare(
      const IntersectionBase<dim, fptype, false> &i) const {
    numIP++;
    fptype delta = intPos - i.intPos;
    return delta;
  }

  int incPrecCount() const { return numIP; }

  static bool cmp(
      const IntersectionBase<dim, fptype, false> &lhs,
      const IntersectionBase<dim, fptype, false> &rhs) {
    fptype diff = lhs.compare(rhs);
    if(diff < 0) return true;
    return false;
  }
};

template <int dim, typename fptype>
class Intersection
    : public IntersectionBase<
          dim, fptype,
          !std::is_same<fptype, mpfr::mpreal>::value> {
 public:
  struct IntersectionData {
    typename Quadric<dim, fptype>::QuadricData quad;
    typename Line<dim, fptype>::LineData line;
    fptype intPos;
    fptype otherIntPos;
    fptype absErrMargin;
  };
  Intersection()
      : IntersectionBase<
            dim, fptype,
            !std::is_same<fptype, mpfr::mpreal>::value>() {}

  Intersection(const Quadric<dim, fptype> &quad,
               const Line<dim, fptype> &line,
               fptype intPos = NAN,
               fptype otherIntPos = NAN,
               fptype absErrMargin = defAbsPrecision)
      : IntersectionBase<
            dim, fptype,
            !std::is_same<fptype, mpfr::mpreal>::value>(
            quad, line, intPos, otherIntPos, absErrMargin) {
  }

  Intersection(
      const IntersectionBase<
          dim, fptype,
          !std::is_same<fptype, mpfr::mpreal>::value> &i)
      : IntersectionBase<
            dim, fptype,
            !std::is_same<fptype, mpfr::mpreal>::value>(i) {
  }

  Intersection(const Intersection<dim, fptype> &i)
      : IntersectionBase<
            dim, fptype,
            !std::is_same<fptype, mpfr::mpreal>::value>(
            static_cast<IntersectionBase<
                dim, fptype,
                !std::is_same<fptype,
                              mpfr::mpreal>::value>>(i)) {}

#ifdef __CUDACC__
  std::shared_ptr<IntersectionData> cudaCopy() const {
    IntersectionData *cudaMem = NULL;
    cudaError_t err =
        cudaMalloc(&cudaMem, sizeof(*cudaMem));
    auto safeMem = std::shared_ptr<IntersectionData>(
        cudaMem, cudaFree);
    err = cudaCopy(cudaMem);
    return safeMem;
  }

  cudaError_t cudaCopy(
      std::shared_ptr<IntersectionData> cudaMem) const {
    cudaError_t err = cudaMemcpy(
        &cudaMem->intPos, &this->intPos,
        sizeof(this->intPos), cudaMemcpyHostToDevice);
    err = cudaMemcpy(
        &cudaMem->otherIntPos, &this->otherIntPos,
        sizeof(this->otherIntPos), cudaMemcpyHostToDevice);
    err = cudaMemcpy(
        &cudaMem->absErrMargin, &this->otherIntPos,
        sizeof(this->absErrMargin), cudaMemcpyHostToDevice);
    err = this->q.cudaCopy(&cudaMem->quad);
    err = this->l.cudaCopy(&cudaMem->line);
    return err;
  }

  cudaError_t cudaCopy(IntersectionData *cudaMem) const {
    cudaError_t err = cudaMemcpy(
        &cudaMem->intPos, &this->intPos,
        sizeof(this->intPos), cudaMemcpyHostToDevice);
    err = cudaMemcpy(
        &cudaMem->otherIntPos, &this->otherIntPos,
        sizeof(this->otherIntPos), cudaMemcpyHostToDevice);
    err = cudaMemcpy(
        &cudaMem->absErrMargin, &this->absErrMargin,
        sizeof(this->absErrMargin), cudaMemcpyHostToDevice);
    err = this->q.cudaCopy(&cudaMem->quad);
    err = this->l.cudaCopy(&cudaMem->line);
    return err;
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<IntersectionData> cudaMem) {
    cudaError_t err = cudaMemcpy(
        &this->intPos, &cudaMem->intPos,
        sizeof(this->intPos), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(
        &this->otherIntPos, &cudaMem->otherIntPos,
        sizeof(this->otherIntPos), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(
        &this->absErrMargin, &cudaMem->absErrMargin,
        sizeof(this->absErrMargin), cudaMemcpyDeviceToHost);
    err = this->q.cudaRetrieve(&cudaMem->quad);
    err = this->l.cudaRetrieve(&cudaMem->line);
    return err;
  }

  cudaError_t cudaRetrieve(IntersectionData *cudaMem) {
    cudaError_t err = cudaMemcpy(
        &this->intPos, &cudaMem->intPos,
        sizeof(this->intPos), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(
        &this->otherIntPos, &cudaMem->otherIntPos,
        sizeof(this->otherIntPos), cudaMemcpyDeviceToHost);
    err = cudaMemcpy(
        &this->absErrMargin, &cudaMem->absErrMargin,
        sizeof(this->absErrMargin), cudaMemcpyDeviceToHost);
    err = this->q.cudaRetrieve(&cudaMem->quad);
    err = this->l.cudaRetrieve(&cudaMem->line);
    return err;
  }
#endif

  friend std::ostream &operator<<(
      std::ostream &os,
      const Intersection<dim, fptype> &i) {
    os << i.q << "\n" << i.l << "\n(" << i.intPos << ", "
       << i.otherIntPos << ")";
    return os;
  }

 private:
  static constexpr const bool isMPFR =
      std::is_same<fptype, mpfr::mpreal>::value;
};

template <typename Iter>
void quicksortInt(Iter start, Iter end) {
  const auto pivot = *start;
  Iter bot = start, top = end;
  bool prevPivot = false;
  /* First rearrange the elements [start, end)
   * into lesser and greater parts, [start, p), [p+1, end)
   */
  for(; bot != top;) {
    if(!prevPivot)
      bot++;
    else
      top--;
    prevPivot = !prevPivot;
    for(; pivot.compare(*bot) < 0.0 && bot != top; bot++) {
    }
    for(; pivot.compare(*top) > 0.0 && bot != top; top--) {
    }
    if(bot == top) break;
    std::iter_swap(bot, top);
  }
  if(bot != end) {
    /* Now rearrange the part larger than the pivot */
    quicksortInt(bot, end);
  }
  /* Put the pivot into the correct place,
   * and then rearrange the part smaller than the pivot */
  bot--;
  std::iter_swap(start, bot);
  if(start != bot) quicksortInt(start, bot);
}

template <int dim, typename fptype>
std::shared_ptr<std::list<Intersection<dim, fptype>>>
sortIntersections(const Line<dim, fptype> &line,
                  std::list<Quadric<dim, fptype>> &quads,
                  fptype absErrMargin = defAbsPrecision) {
  using IP = Intersection<dim, fptype>;
  std::shared_ptr<std::list<IP>> inter(new std::list<IP>);
  for(auto q : quads) {
    auto intPos = q.calcLineDistToIntersect(line);
    for(int k = 0; k < 2; k++) {
      if(!MathFuncs::MathFuncs<fptype>::isnan(intPos[k])) {
        Intersection<dim, fptype> i(q, line, intPos[k],
                                    intPos[1 - k],
                                    absErrMargin);
        inter->push_back(i);
      }
    }
  }
  inter->sort([](IP &lhs, IP &rhs) {
    return lhs.compare(rhs) < 0.0;
  });
  /* TODO: Fix this */
  // quicksortInt(inter->begin(), inter->end());
  return inter;
}
}
