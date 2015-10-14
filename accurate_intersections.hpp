
#include "geometry.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "polynomial.hpp"
#include "genericfp.hpp"

#include <cmath>
#include <list>
#include <memory>
#include <utility>

#include "mpreal.hpp"

namespace Geometry {

template <int dim, typename fptype>
class Intersection {
 public:
  const Quadric<dim, fptype> *q;
  const Line<dim, fptype> *l;
  fptype intPos, otherIntPos;
  fptype absErrMargin;

  Intersection(const Quadric<dim, fptype> &quad,
               const Line<dim, fptype> &line,
               fptype intPos = NAN,
               fptype otherIntPos = NAN,
               fptype absErrMargin = defAbsPrecision)
      : q(&quad),
        l(&line),
        intPos(intPos),
        otherIntPos(otherIntPos),
        absErrMargin(absErrMargin) {}
  Intersection(const Intersection<dim, fptype> &i)
      : q(i.q),
        intPos(i.intPos),
        otherIntPos(i.otherIntPos),
        absErrMargin(i.absErrMargin) {}

  Intersection<dim, fptype> operator=(
      const Intersection<dim, fptype> &i) {
    q = i.q;
    intPos = i.intPos;
    otherIntPos = i.otherIntPos;
    absErrMargin = i.absErrMargin;
    return *this;
  }

  fptype compare(const Intersection<dim, fptype> &i) const {
    fptype delta = intPos - i.intPos;
    if(std::abs(delta) < absErrMargin) {
      /* Compute the coefficients for the determinant
       * We need a precision of 3 times machine precision
       * to compute the coefficients accurately,
       * so set that first
       */
      const unsigned machPrec =
          GenericFP::fpconvert<fptype>::precision;
      const unsigned coeffPrec = 3 * machPrec;
      mpfr::mpreal::set_default_prec(coeffPrec);
      /* Compute the coefficients */
      Quadric<dim, mpfr::mpreal> q1(*q);
      Quadric<dim, mpfr::mpreal> q2(*i.q);
      Line<dim, mpfr::mpreal> line(*l);
      line.shiftOrigin(q->getOrigin());
      Polynomial<2, mpfr::mpreal> p1 =
          q1.calcLineDistPoly(line);
      Polynomial<2, mpfr::mpreal> p2 =
          q2.calcLineDistPoly(line);
      /* And the determinant */
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
      fptype termSigns[numTerms] = {1, 1, 1, 1, -1, -1, -1};
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
      fptype tmp = det;
      tmp *= partialResultant;
      if(tmp < 0.0) {
        /* Opposite signs! */
        return -1.0;
      }
      else if(tmp > 0.0) {
        /* Same Signs! */
        return 1.0;
      }
      /* We're in trouble */
      assert(tmp != 0.0);
    }
    return delta;
  }
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
  /* Now rearrange the part larger than the pivot */
  quicksortInt(bot, end);
  /* Put the pivot into the correct place,
   * and then rearrange the part smaller than the pivot */
  bot--;
  std::iter_swap(start, bot);
  quicksortInt(start, bot);
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
      if(!std::isnan(intPos[k])) {
        struct Intersection<dim, fptype> i(
            q, line, intPos[k], intPos[!k], absErrMargin);
        inter->push_back(i);
      }
    }
  }
  quicksortInt(inter->begin(), inter->end());
  return inter;
}
}
