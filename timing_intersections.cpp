
#include "geometry.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "accurate_intersections.hpp"
#include "timer.hpp"

#include <list>
#include <fstream>
#include <iomanip>

template <int dim, typename fptype>
std::list<Geometry::Quadric<dim, fptype>> parseQuadrics(
    const char *fname) {
  std::ifstream file(fname);
  using Qf = Geometry::Quadric<dim, fptype>;
  std::list<Qf> quads;
  while(!file.eof()) {
    Qf q;
    file >> q;
    quads.push_back(q);
  }
  return quads;
}

bool isSameInt(mpfr::mpreal int1, mpfr::mpreal int2) {
  /* Use a heuristic on truth to determine whether adjacent
   * intersections occur at the same place.
   * This will basically be a threshold on the largest
   * different bit in the mantissa.
   */
  mpfr::mpreal largest =
      mpfr::max(mpfr::fabs(int1), mpfr::fabs(int2));
  mp_exp_t largestExp = largest.get_exp();
  mpfr::mpreal diff = mpfr::fabs(int1 - int2);
  int p1 = int1.getPrecision(), p2 = int2.getPrecision();
  int minPrec = std::min(p1, p2);
  mp_exp_t deltaExp =
      largestExp - minPrec / 4 - diff.get_exp();
  return deltaExp > 0;
}

void validateResults(std::ostream &results, auto &inter,
                     auto &truth) {
  auto j = truth->begin();
  for(auto i = inter->begin();
      i != inter->end() || j != truth->end();) {
    bool printQ = false;
    if(j != truth->end()) {
      /* Create the region boundaries [sameBeg, sameEnd).
       * These are intersections which are probably the same
       */
      auto sameBeg = j;
      auto sameEnd = j;
      int regLen = 0;
      while(sameEnd != truth->end() &&
            isSameInt(j->intPos, sameEnd->intPos)) {
        sameEnd++;
        regLen++;
      }
      /* Increment i until it's not found in the region */
      int numInRegion = 0;
      bool isInRegion = true;
      while(i != inter->end() && isInRegion) {
        j = sameBeg;
        while(j != sameEnd && i->q != j->q) j++;
        if(j == sameEnd) {
          isInRegion = false;
        } else {
          i++;
          numInRegion--;
        }
      }
      /* i isn't in the region.
       * Verify all elements in the region were used */
      if(regLen != numInRegion) {
        printQ = true;
      }
    }
    /* And output the result */
    if(printQ) {
      results << "Incorrect Result\n";
      results << "Expected: " << std::setprecision(20)
              << j->intPos
              << "; Got: " << std::setprecision(20)
              << i->intPos
              << "\nDelta: " << std::setprecision(20)
              << j->intPos - i->intPos << "\n";
      if(i != inter->end()) {
        results << "Estimated: " << i->q << "\n";
        i++;
      } else {
        results << "Computed intersections ended "
                   "prematurely\n";
      }
      if(j != truth->end()) {
        results << "Correct:   " << j->q << "\n";
        j++;
      } else {
        results
            << "Computed intersections ended too late\n";
      }
    } else {
      i++;
      j++;
    }
  }
  results << "\n";
}

template <int dim, typename fptype>
void intersectionTest(
    std::list<Geometry::Quadric<dim, fptype>> &quads,
    const int numTests) {
  /* First build a scene of quadrics.
   * Then generate random lines on a disk centered at the
   * intersection of the cylinders.
   * Then sort the intersections
   * Finally validate them with a higher precision sort
   */
  using Vf = Geometry::Vector<dim, fptype>;
  using Pf = Geometry::Point<dim, fptype>;
  using Lf = Geometry::Line<dim, fptype>;
  using Qm = Geometry::Quadric<dim, mpfr::mpreal>;
  using Lm = Geometry::Line<dim, mpfr::mpreal>;
  constexpr const int machPrec =
      GenericFP::fpconvert<fptype>::precision;
  constexpr const int truthPrec = 48 * machPrec;
  mpfr::mpreal::set_default_prec(truthPrec);
  std::list<Qm> truthQuads;
  /* Generate the quadrics */
  for(auto q : quads) {
    Qm quadMP(q);
    truthQuads.push_back(quadMP);
  }
  std::random_device rd;
  std::mt19937_64 engine(rd());
  constexpr const fptype maxPos = 10000;
  std::uniform_real_distribution<fptype> genPos(-maxPos,
                                                maxPos);
  std::uniform_real_distribution<fptype> genDir(-1.0, 1.0);
  std::ofstream results("results");
  Timer::Timer fp_time, mp_time;
  /* Run the tests */
  for(int t = 0; t < numTests; t++) {
    /* First build the line */
    Vf lineDir;
    Vf lineInt;
    for(int i = 0; i < dim; i++) {
      fptype tmp = genPos(engine);
      lineInt.set(i, tmp);
      tmp = genDir(engine);
      lineDir.set(i, tmp);
    }
    Lf line(Pf(lineInt), lineDir);
    Lm truthLine(line);
    constexpr const fptype eps = 1e-3;
    /* Then sort the intersections */
    fp_time.startTimer();
    auto inter =
        Geometry::sortIntersections(line, quads, eps);
    fp_time.stopTimer();
    mp_time.startTimer();
    auto truth =
        Geometry::sortIntersections<dim, mpfr::mpreal>(
            truthLine, truthQuads, eps);
    mp_time.stopTimer();
    results << "Test " << t << "\n";
    results << "FP Time:\n" << fp_time << "\n";
    results << "MP Time:\n" << mp_time << "\n";
    validateResults(results, inter, truth);
  }
}

int main(int argc, char **argv) {
  using fptype = float;
  constexpr const int dim = 3;
  std::list<Geometry::Quadric<dim, fptype>> quads;
  if(argc > 1)
    quads = parseQuadrics<dim, fptype>(argv[1]);
  else
    quads = parseQuadrics<dim, fptype>("cylinders.csg");
  intersectionTest(quads, 1e5);
  return 0;
}
