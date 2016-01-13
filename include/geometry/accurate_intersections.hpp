
#include "geometry.hpp"
#include "quadrics.hpp"
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
  const Quadric<dim, fptype> q;
  const Line<dim, fptype> l;
  fptype intPos, otherIntPos;
  fptype absErrMargin;
  int numIP;

  IntersectionBase()
      : q(),
        l(Point<dim, fptype>(), Vector<dim, fptype>()),
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

  fptype accurateCompare(
      const IntersectionBase<dim, fptype, true> &i) {
    numIP++;
    /* Only works when there are two roots for both quadrics
     * Also requires the differences of roots to be
     * greater than the minimum precision.
     * Compute the coefficients for the determinant
     * We need a precision of 3 times machine precision
     * to compute the coefficients accurately,
     * so set that first
     */
    constexpr const unsigned machPrec =
        GenericFP::fpconvert<fptype>::precision;
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
    /* And the determinant */
    const unsigned prevPrec =
        mpfr::mpreal::get_default_prec();
    const unsigned detPrec = 4 * coeffPrec;
    mpfr::mpreal::set_default_prec(detPrec);
    mpfr::mpreal coeffs[] = {p1.get(2), p1.get(1),
                             p1.get(0), p2.get(2),
                             p2.get(1), p2.get(0)};
    constexpr const int numTerms = 7;
    constexpr const int numProds = 4;
    int detTermCoeffs[numTerms][numProds] = {
        {0, 0, 5, 5}, {2, 2, 3, 3}, {1, 1, 3, 5},
        {4, 4, 0, 2}, {1, 2, 3, 4}, {0, 1, 4, 5},
        {0, 2, 3, 5}};
    fptype termSigns[numTerms] = {1, 1, 1, 1, -1, -1, -2};
    mpfr::mpreal det(0.0);
    for(int i = 0; i < numTerms; i++) {
      mpfr::mpreal term(termSigns[i]);
      for(int j = 0; j < numProds; j++) {
        int c = detTermCoeffs[i][j];
        term *= coeffs[c];
      }
      det += term;
    }
    /* Now compute the signs of the resultant terms */
    constexpr const int numPartialTerms = 3;
    fptype partResTerms[numPartialTerms] = {
        intPos - i.otherIntPos, otherIntPos - i.intPos,
        otherIntPos - i.otherIntPos};
    fptype partialResultant = 1.0;
    for(int i = 0; i < numPartialTerms; i++)
      partialResultant *= partResTerms[i];
    /* If the sign of this matches the determinants sign,
     * then this intersection occurs after the other one,
     * otherwise, this intersection occurs before */
    fptype tmp = static_cast<fptype>(det);
    tmp *= partialResultant;
    mpfr::mpreal::set_default_prec(prevPrec);
    if(tmp < 0.0) {
      /* Opposite signs! */
      return -1.0;
    } else if(tmp > 0.0) {
      /* Same Signs! */
      return 1.0;
    }
    /* We can't distinguish the two intersections */
    return 0.0;
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
  const Quadric<dim, fptype> q;
  const Line<dim, fptype> l;
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
      : IntersectionBase<dim, fptype, !std::is_same<fptype, mpfr::mpreal>::value>() {}

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
    auto safeMem =
        std::shared_ptr<IntersectionData>(cudaMem, cudaFree);
    err = cudaCopy(cudaMem);
    return safeMem;
  }

  cudaError_t cudaCopy(
      std::shared_ptr<IntersectionData> cudaMem) const {
    cudaError_t err =
        cudaMemcpy(&cudaMem->intPos, &this->intPos,
                   sizeof(this->intPos), cudaMemcpyHostToDevice);
    err =
        cudaMemcpy(&cudaMem->otherIntPos, &this->otherIntPos,
                   sizeof(this->otherIntPos), cudaMemcpyHostToDevice);
    err =
        cudaMemcpy(&cudaMem->absErrMargin, &this->otherIntPos,
                   sizeof(this->absErrMargin), cudaMemcpyHostToDevice);
    err = this->q.cudaCopy(&cudaMem->quad);
		err = this->l.cudaCopy(&cudaMem->line);
    return err;
  }

  cudaError_t cudaCopy(
      IntersectionData *cudaMem) const {
    cudaError_t err =
        cudaMemcpy(&cudaMem->intPos, &this->intPos,
                   sizeof(this->intPos), cudaMemcpyHostToDevice);
    err =
        cudaMemcpy(&cudaMem->otherIntPos, &this->otherIntPos,
                   sizeof(this->otherIntPos), cudaMemcpyHostToDevice);
    err =
        cudaMemcpy(&cudaMem->absErrMargin, &this->absErrMargin,
                   sizeof(this->absErrMargin), cudaMemcpyHostToDevice);
    err = this->q.cudaCopy(&cudaMem->quad);
		err = this->l.cudaCopy(&cudaMem->line);
    return err;
  }

  cudaError_t cudaRetrieve(
      std::shared_ptr<IntersectionData> cudaMem) {
    IntersectionData buffer;
    cudaError_t err = cudaMemcpy(
        &buffer.intPos, &cudaMem->intPos,
        sizeof(buffer.intPos), cudaMemcpyDeviceToHost);
		err = cudaMemcpy(
        &buffer.otherIntPos, &cudaMem->otherIntPos,
        sizeof(buffer.otherIntPos), cudaMemcpyDeviceToHost);
		err = cudaMemcpy(
        &buffer.absErrMargin, &cudaMem->absErrMargin,
        sizeof(buffer.absErrMargin), cudaMemcpyDeviceToHost);
		return err;
  }

  cudaError_t cudaRetrieve(
      IntersectionData *cudaMem) {
    IntersectionData buffer;
    cudaError_t err = cudaMemcpy(
        &buffer.intPos, &cudaMem->intPos,
        sizeof(buffer.intPos), cudaMemcpyDeviceToHost);
		err = cudaMemcpy(
        &buffer.otherIntPos, &cudaMem->otherIntPos,
        sizeof(buffer.otherIntPos), cudaMemcpyDeviceToHost);
		err = cudaMemcpy(
        &buffer.absErrMargin, &cudaMem->absErrMargin,
        sizeof(buffer.absErrMargin), cudaMemcpyDeviceToHost);
		return err;
  }
	#endif

	friend std::ostream &operator<<(std::ostream &os,
																	const Intersection<dim, fptype> &i) {
		os << i.q << "\n" << i.l << "\n(" << i.intPos << ", " << i.otherIntPos << ")";
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
        Intersection<dim, fptype> i(
            q, line, intPos[k], intPos[1 - k],
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
