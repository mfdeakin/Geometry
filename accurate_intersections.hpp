
#include "geometry.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "polynomial.hpp"

#include <cmath>
#include <list>
#include <memory>
#include <utility>

namespace Geometry {

template <int dim, typename fptype>
class Intersection {
 public:
  const Quadric<dim, fptype> *q;
  fptype intPos, otherIntPos;
  fptype absErrMargin;

  Intersection(const Quadric<dim, fptype> &quad,
               fptype intPos = NAN,
               fptype otherIntPos = NAN,
               fptype absErrMargin = defAbsPrecision)
      : q(&quad),
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
            q, intPos[k], intPos[!k], absErrMargin);
        inter->push_back(i);
      }
    }
  }
  quicksortInt(inter->begin(), inter->end());
  return inter;
}
}
