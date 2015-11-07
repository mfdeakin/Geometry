
#include "geometry.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "accurate_intersections.hpp"
#include "timer.hpp"

#include <list>
#include <fstream>
#include <iomanip>

#include <signal.h>

struct GlobalVars {
  bool run = true;
} globals;

void sigInt(int signum) { globals.run = false; }

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
      largestExp - minPrec * 3 / 4 - diff.get_exp();
  return deltaExp > 0;
}

bool validateResults(auto &inter, auto &truth) {
  auto j = truth->begin();
  for(auto i = inter->begin();
      i != inter->end() || j != truth->end();) {
    /* First determine the length of the region of equal
     * intersections.
     * Then verify there's an equal number of
     * intersections in the approximately equal range and
     * that they correspond to the correct intersections.
     * If not, then print the intersection
     */
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
      while(i != inter->end() && isInRegion &&
            numInRegion < regLen) {
        j = sameBeg;
        while(j != sameEnd && i->q != j->q) j++;
        if(j == sameEnd) {
          isInRegion = false;
        } else {
          i++;
          numInRegion++;
        }
      }
      /* i isn't in the region.
       * Verify all elements in the region were used */
      if(regLen != numInRegion) return false;
      j = sameEnd;
    } else {
      return false;
      i++;
    }
  }
  return true;
}

template <int dim, typename fptype>
void intersectionTest(
    std::list<Geometry::Quadric<dim, fptype>> &quads,
    std::ostream &results, const int numTests) {
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
  Timer::Timer fp_time, mp_time;
  struct TimeArr {
    int fpns;
    int mpns;
    bool correct;
  } *times = new struct TimeArr[numTests];
  /* Run the tests */
  int t;
  for(t = 0; t < numTests && globals.run; t++) {
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
    times[t].fpns = fp_time.instant_ns();
    times[t].mpns = mp_time.instant_ns();
    times[t].correct = validateResults(inter, truth);
  }
  /* Output all of the results */
  results
      << "Test #, FP Time (ns), MP Time (ns), Correct\n";
  for(int i = 0; i < t; i++)
    results << i + 1 << ", " << times[i].fpns << ", "
            << times[i].mpns << ", " << times[i].correct
            << "\n";
  results << "\n"
          << "Total FP Time (s): " << fp_time.elapsed_s()
          << "." << std::setw(9) << std::setfill('0')
          << fp_time.elapsed_ns() << "\n"
          << "Total MP Time (s): " << mp_time.elapsed_s()
          << "." << std::setw(9) << std::setfill('0')
          << mp_time.elapsed_ns() << "\n";
  delete[] times;
}

int main(int argc, char **argv) {
  using fptype = float;
  constexpr const int dim = 3;
  std::list<Geometry::Quadric<dim, fptype>> quads;
  const char *outFName = "results";
  int numTests = 1e5;
  if(argc > 1) {
    quads = parseQuadrics<dim, fptype>(argv[1]);
    if(argc > 2) {
      outFName = argv[2];
      if(argc > 3) numTests = atoi(argv[3]);
    }
  } else
    quads = parseQuadrics<dim, fptype>("cylinders.csg");
  std::ofstream results(outFName);
  signal(SIGINT, sigInt);
  intersectionTest(quads, results, numTests);
  return 0;
}
