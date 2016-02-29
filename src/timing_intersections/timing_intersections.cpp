
#include "geometry.hpp"
#include "quadric.hpp"
#include "line.hpp"
#include "accurate_intersections.hpp"
#include "quadric_classify.hpp"
#include "timer.hpp"

#include <list>
#include <fstream>
#include <iomanip>

#include <signal.h>
#include <unistd.h>

struct GlobalVars {
  bool run = true;
} globals;

void sigInt(int signum) { globals.run = false; }

template <int dim, typename fptype>
std::list<Geometry::Quadric<dim, fptype>> parseQuadrics(
    const char *fname) {
  std::ifstream file(fname);
  using Qf = Geometry::Quadric<dim, fptype>;
  int qtypeCount[QuadricClassify::QUADT_ERRORINVALID];
  for(int i = 0; i < QuadricClassify::QUADT_ERRORINVALID;
      i++)
    qtypeCount[i] = 0;
  std::list<Qf> quads;
  int imQuads = 0;
  while(!file.eof()) {
    Qf q;
    file >> q;
    QuadricClassify::QuadType type =
        QuadricClassify::classifyQuadric(q);
    if(!QuadricClassify::isImaginary(type) &&
       type != QuadricClassify::QUADT_ERROR) {
      quads.push_back(q);
      qtypeCount[type]++;
    } else {
      imQuads++;
    }
  }
  for(int i = 0; i < QuadricClassify::QUADT_ERRORINVALID;
      i++) {
    std::cout << qtypeCount[i] << " "
              << QuadricClassify::QuadTypeNames[i] << "\n";
  }
  return quads;
}

bool isSameInt(mpfr::mpreal int1, mpfr::mpreal int2) {
  /* Use a heuristic on truth to determine whether adjacent
   * intersections occur at the same place.
   * This will basically be a threshold on the largest
   * different bit in the mantissa.
   */
  /* A relative heuristic does not work when one (or both!)
   * is 0 */
  if(int1 == 0.0 || int2 == 0.0) {
    constexpr const double eps = 0.000001;
    return (mpfr::fabs(int1 - int2) < eps);
  }
  mpfr::mpreal largest =
      mpfr::max(mpfr::fabs(int1), mpfr::fabs(int2));
  mp_exp_t largestExp = largest.get_exp();
  mpfr::mpreal diff = mpfr::fabs(int1 - int2);
  int p1 = int1.getPrecision(), p2 = int2.getPrecision();
  int minPrec = std::min(p1, p2);
  mp_exp_t maxExp = largestExp - minPrec * 1 / 48,
           diffExp = diff.get_exp();
  return maxExp >= diffExp;
}

template <typename ListFP, typename ListMP>
bool validateResults(ListFP &inter, ListMP &truth) {
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
        /* sameEnd is not in the region, so if j reaches it,
         * i isn't in the region */
        if(j == sameEnd) {
          isInRegion = false;
        } else {
          i++;
          numInRegion++;
        }
      }
      /* i isn't in the region.
       * Verify all elements in the region were used */
      if(regLen != numInRegion) {
        return false;
      }
      j = sameEnd;
    } else {
      int numRemaining = 0;
      while(i != inter->end()) {
        i++;
        numRemaining++;
      }
      return false;
    }
  }
  return true;
}

template <typename List>
int countIP(List inter) {
  int numIP = 0;
  for(auto i : *inter) numIP += i.incPrecCount();
  return numIP;
}

template <int dim, typename fptype>
using randLineGen =
    Geometry::Line<dim, fptype> (*)(std::mt19937_64 &rng);

template <int dim, typename fptype>
std::shared_ptr<std::list<Geometry::Intersection<
    dim, fptype>>> __attribute__((noinline))
runTest(std::list<Geometry::Quadric<dim, fptype>> &quads,
        Geometry::Line<dim, fptype> &line,
        Timer::Timer &timer) {
  constexpr const double eps =
      std::numeric_limits<double>::infinity();
  timer.startTimer();
  auto inter =
      Geometry::sortIntersections(line, quads, fptype(eps));
  timer.stopTimer();
  return inter;
}

template <int dim, typename fptype>
void intersectionTest(
    std::list<Geometry::Quadric<dim, fptype>> &quads,
    std::ostream &results, const int numTests,
    randLineGen<dim, fptype> rlgf) {
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
  constexpr const int precMult = 24;
  constexpr const int truthPrec = precMult * machPrec;
  mpfr::mpreal::set_default_prec(truthPrec);
  std::list<Qm> truthQuads;
  /* Generate the quadrics */
  for(auto q : quads) {
    Qm quadMP(q);
    truthQuads.push_back(quadMP);
  }
  std::random_device rd;
  std::mt19937_64 engine(rd());
  Timer::Timer fp_time, mp_time;
  struct TimeArr {
    int fpns;
    int mpns;
    bool correct;
    int numIP;
  } *times = new struct TimeArr[numTests];
  /* Run the tests */
  int t;
  for(t = 0; t < numTests && globals.run; t++) {
    Lf line = rlgf(engine);
    Lm truthLine(line);
    constexpr const fptype eps =
        std::numeric_limits<fptype>::infinity();
    /* Then sort the intersections */
    auto inter = runTest(quads, line, fp_time);
    auto truth = runTest(truthQuads, truthLine, mp_time);
    times[t].fpns = fp_time.instant_ns();
    times[t].mpns = mp_time.instant_ns();
    times[t].correct = validateResults(inter, truth);
    times[t].numIP = countIP(inter);
  }
  /* Output all of the results */
  results << "Test #, FP Time (ns), MP Time (ns), Correct, "
             "Precision Increases\n";
  int totIP = 0;
  for(int i = 0; i < t; i++) {
    results << i + 1 << ", " << times[i].fpns << ", "
            << times[i].mpns << ", " << times[i].correct
            << ", " << times[i].numIP << "\n";
    totIP += times[i].numIP;
  }
  results << "\n"
          << "Total FP Time (s): " << fp_time.elapsed_s()
          << "." << std::setw(9) << std::setfill('0')
          << fp_time.elapsed_ns() << "\n"
          << "Total MP Time (s): " << mp_time.elapsed_s()
          << "." << std::setw(9) << std::setfill('0')
          << mp_time.elapsed_ns() << "\n";
  results << "Total precision increases: " << totIP << "\n";
  delete[] times;
}

void lockCPU() {
  const int numCPUs = sysconf(_SC_NPROCESSORS_ONLN);
  const int cpuSets =
      numCPUs / CPU_SETSIZE + ((numCPUs % CPU_SETSIZE) > 0);
  cpu_set_t *cpus = new cpu_set_t[cpuSets];
  const size_t cpuSize = sizeof(cpu_set_t[cpuSets]);
  sched_getaffinity(0, cpuSize, cpus);
  for(int i = 1; i < numCPUs; i++) CPU_CLR(i, cpus);
  CPU_SET(0, cpus);
  sched_setaffinity(0, cpuSize, cpus);
  delete[] cpus;
}

template <int dim, typename fptype>
Geometry::Line<dim, fptype> defRandLine(
    std::mt19937_64 &rng) {
  constexpr const fptype minPos = 0, maxPos = 2;
  std::uniform_real_distribution<fptype> genPos(minPos,
                                                maxPos);
  std::uniform_real_distribution<fptype> genDir(-1.0, 1.0);
  /* First build the line */
  Geometry::Vector<dim, fptype> lineDir;
  Geometry::Vector<dim, fptype> lineInt;
  for(int i = 0; i < dim; i++) {
    fptype tmp = genPos(rng);
    lineInt.set(i, tmp);
    tmp = genDir(rng);
    lineDir.set(i, tmp);
  }
  return Geometry::Line<dim, fptype>(
      Geometry::Point<dim, fptype>(lineInt), lineInt);
}

template <int dim, typename fptype>
Geometry::Line<dim, fptype> cylRandLine(
    std::mt19937_64 &rng) {
  constexpr const fptype minPos = 0.375, maxPos = 0.625;
  std::uniform_real_distribution<fptype> genPos(minPos,
                                                maxPos);
  std::uniform_real_distribution<fptype> genDir(-1.0, 1.0);
  /* First build the line */
  Geometry::Vector<dim, fptype> lineDir;
  Geometry::Vector<dim, fptype> lineInt;
  lineInt.set(1, genPos(rng));
  for(int i = 0; i < dim; i++) {
    lineInt.set(i, lineInt.get(1));
    lineDir.set(i, genDir(rng));
  }
  lineInt.set(0, genPos(rng));
  return Geometry::Line<dim, fptype>(
      Geometry::Point<dim, fptype>(lineInt), lineInt);
}

int main(int argc, char **argv) {
  using fptype = float;
  constexpr const int dim = 3;
  lockCPU();
  std::list<Geometry::Quadric<dim, fptype>> quads;
  const char *outFName = "results";
  int numTests = 1e4;
  randLineGen<dim, fptype> rlg = cylRandLine<dim, fptype>;
  if(argc > 1) {
    quads = parseQuadrics<dim, fptype>(argv[1]);
    rlg = defRandLine<dim, fptype>;
    if(argc > 2) {
      outFName = argv[2];
      if(argc > 3) numTests = atoi(argv[3]);
    }
  } else {
    quads = parseQuadrics<dim, fptype>("cylinders.csg");
  }
  std::ofstream results(outFName);
  signal(SIGINT, sigInt);
  intersectionTest(quads, results, numTests, defRandLine);
  return 0;
}
