
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
  int numIP;

  IntersectionBase()
      : q(),
        intPos(NAN),
        otherIntPos(NAN),
        absErrMargin(defAbsPrecision),
        numIP(0) {}

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
      const IntersectionBase<dim, fptype, true> &i)
      : q(i.q),
        l(i.l),
        intPos(i.intPos),
        otherIntPos(i.otherIntPos),
        absErrMargin(i.absErrMargin),
        numIP(i.numIP) {}

  IntersectionBase<dim, fptype, true> operator=(
      const IntersectionBase<dim, fptype, true> &i) {
    q = i.q;
    l = i.l;
    intPos = i.intPos;
    otherIntPos = i.otherIntPos;
    absErrMargin = i.absErrMargin;
    return *this;
  }

  mpfr::mpreal resultantDet(
      const IntersectionBase<dim, fptype, true> &i) {
    /* Compute the coefficients for the determinant
     * We need a precision of 3 times machine precision
     * to compute the coefficients accurately,
     * so set that first
     */
    const unsigned prevPrec =
        mpfr::mpreal::get_default_prec();
    constexpr const unsigned machPrec =
        GenericFP::fpconvert<fptype>::precision;
    /* This is less efficient than we'd like,
     * but is the only way to implement this in a
     * manner which is agnostic floating point type
     */
    constexpr const unsigned coeffPrec = 3 * machPrec;
    mpfr::mpreal::set_default_prec(coeffPrec);
    /* Compute the coefficients */
    Quadric<dim, mpfr::mpreal> q1(q);
    Quadric<dim, mpfr::mpreal> q2(i.q);
    Line<dim, mpfr::mpreal> line(l);
    line.shiftOrigin(q.getOrigin());
    Polynomial<2, mpfr::mpreal> p1 =
        q1.calcLineDistPoly(line);
    Polynomial<2, mpfr::mpreal> p2 =
        q2.calcLineDistPoly(line);
    /* And the determinant
     * The determinant requires 4 times the coefficients
     * precision to avoid multiplicative errors
     */
    /* The fastest way to compute our determinant
     * 1.0(0 0 5 5)+1.0(2 2 3 3)+1.0(1 1 3 5)+1.0(4 4 0 2)+
     * -1.0(1 2 3 4)-1.0(0 1 4 5)-2.0(0 2 3 5)
     * 35 FLOPs
     * (0 5)((0 5)-2.0(2 3)-(1 4))+(2 3)((2 3)-(1 4))+
     * (1 1 3 5)+(4 4 0 2)
     *  3+1+2+1                    +0+0+1+1
     * +3       +3
     * +3
     * 18 FLOPs, reduced by a factor of ~1/2,
     * in exchange for ugly code
     *
     * Start by computing the partial products
     */
    const unsigned partPrec = 2 * coeffPrec;
    constexpr const int numCoeffs = 6;
    mpfr::mpreal coeffs[numCoeffs] = {
        p1.get(2).setPrecision(partPrec, MPFR_RNDN),
        p1.get(1).setPrecision(partPrec, MPFR_RNDN),
        p1.get(0).setPrecision(partPrec, MPFR_RNDN),
        p2.get(2).setPrecision(partPrec, MPFR_RNDN),
        p2.get(1).setPrecision(partPrec, MPFR_RNDN),
        p2.get(0).setPrecision(partPrec, MPFR_RNDN)};

    constexpr const int numPP = 7;
    mpfr::mpreal partialProds[numPP];
    partialProds[0] = coeffs[0] * coeffs[5];
    partialProds[1] = coeffs[2] * coeffs[3];
    partialProds[2] = coeffs[1] * coeffs[4];

    partialProds[3] = sqr(coeffs[1]);
    partialProds[4] = coeffs[3] * coeffs[5];

    partialProds[5] = sqr(coeffs[4]);
    partialProds[6] = coeffs[0] * coeffs[2];

    /* Increase the precision for the multiplication
     * Again, an inefficiency :(
     */
    const unsigned detPrec = 2 * partPrec;
    for(int i = 0; i < numPP; i++) {
      partialProds[i].setPrecision(detPrec, MPFR_RNDN);
    }

    /* Now compute the larger terms */
    constexpr const int numTerms = 4;
    mpfr::mpreal detTerms[numTerms];
    // (0 5)((0 5)-2.0(2 3)-(1 4))
    detTerms[0] = (partialProds[0] - (partialProds[1] << 1) -
                   partialProds[2]) *
                  partialProds[0];
    // (2 3)((2 3)-(1 4))
    detTerms[1] = (partialProds[1] - partialProds[2]) *
                  partialProds[1];
    // (1 1 3 5)
    detTerms[2] = partialProds[3] * partialProds[4];
    // (4 4 0 2)
    detTerms[3] = partialProds[5] * partialProds[6];
    /* There are only 4 terms to sum,
     * which isn't enough for compensated summation to be
     * worthwhile in my experience
     */
    mpfr::mpreal det(0.0);
    for(int i = 0; i < numTerms; i++) {
      det += detTerms[i];
    }
    mpfr::mpreal::set_default_prec(prevPrec);
    return det;
  }

  fptype accurateCompare(
      const IntersectionBase<dim, fptype, true> &i) {
    /* Only works when there are two roots for both quadrics
     * Also requires the differences of roots to be
     * greater than the minimum precision.
     */

    /* Keep track of the number of precision increases
     * required for statistics
     */
    numIP++;
    mpfr::mpreal det = resultantDet(i);
    /* Now compute the signs of the other resultant terms.
     * Since the sign flips with each negative, just
     * use XOR on one bit to keep track of the final sign
     */
    constexpr const int numPartialTerms = 3;
    fptype partResTerms[numPartialTerms] = {
        intPos - i.otherIntPos, otherIntPos - i.intPos,
        otherIntPos - i.otherIntPos};
    int numNeg = 0;
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
    if(det < 0) {
      numNeg ^= 1;
    }
    if(numNeg == 0) {
      return 1.0;
    } else {
      return -1.0;
    }
  }

  fptype compare(
      const IntersectionBase<dim, fptype, true> &i) {
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

  int incPrecCount() { return numIP; }

  static bool cmp(
      IntersectionBase<dim, fptype, true> &lhs,
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
        absErrMargin(absErrMargin) {}
  IntersectionBase(
      const IntersectionBase<dim, fptype, false> &i)
      : q(i.q),
        l(i.l),
        intPos(i.intPos),
        otherIntPos(i.otherIntPos),
        absErrMargin(i.absErrMargin) {}

  IntersectionBase<dim, fptype, false> operator=(
      const IntersectionBase<dim, fptype, false> &i) {
    q = i.q;
    l = i.l;
    intPos = i.intPos;
    otherIntPos = i.otherIntPos;
    absErrMargin = i.absErrMargin;
    return *this;
  }

  fptype compare(
      const IntersectionBase<dim, fptype, false> &i) const {
    fptype delta = intPos - i.intPos;
    return delta;
  }

  int incPrecCount() { return 0; }

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
  inter->sort([](IP &lhs, const IP &rhs) {
    return lhs.compare(rhs) < 0.0;
  });
  /* TODO: Fix this */
  // quicksortInt(inter->begin(), inter->end());
  return inter;
}
}
