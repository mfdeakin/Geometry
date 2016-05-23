
#include "accurate_intersections_approx.hpp"
#include "accurate_intersections_incprec.hpp"
#include "accurate_intersections_resultant.hpp"
#include "geometry.hpp"
#include "line.hpp"
#include "quadric.hpp"
#include "quadric_classify.hpp"
#include "timer.hpp"

#include <array>
#include <fstream>
#include <iomanip>
#include <list>
#include <random>

#include <signal.h>
#include <unistd.h>

struct GlobalVars {
  bool run = true;
} globals;

void sigInt(int signum) { globals.run = false; }

template <int dim, typename fptype>
std::list<Geometry::Quadric<dim, fptype>> parseQuadrics(
    const char *fname, int *minExp, int *maxExp) {
  std::ifstream file(fname);
  if(!file.is_open()) {
    return std::list<Geometry::Quadric<dim, fptype>>();
  }
  int buf = 0;
  file >> buf;
  if(minExp != NULL) {
    *minExp = buf;
  }
  file >> buf;
  if(maxExp != NULL) {
    *maxExp = buf;
  }
  if(*minExp > 0) {
    std::cout
        << "Cannot generate a line with these exponents\n"
        << *minExp << ", " << *maxExp << "\n";
    exit(1);
  }
  std::cout << "Using " << *minExp << ", " << *maxExp
            << " as the range of exponents\n";

  using Qf = Geometry::Quadric<dim, fptype>;
  int qtypeCount[QuadricClassify::QUADT_ERRORINVALID];
  for(int i = 0; i < QuadricClassify::QUADT_ERRORINVALID;
      i++)
    qtypeCount[i] = 0;
  std::list<Qf> quads;
  int imQuads = 0;
  int numQuads = 0;
  while(!file.eof()) {
    Qf q;
    file >> q;
    QuadricClassify::QuadType type =
        QuadricClassify::classifyQuadric(q);
    if(!QuadricClassify::isImaginary(type) &&
       type != QuadricClassify::QUADT_ERROR &&
       type != QuadricClassify::QUADT_DEGENERATE &&
       type != QuadricClassify::QUADT_ERRORINVALID) {
      quads.push_back(q);
      qtypeCount[type]++;
    } else {
      std::cout << "Quadric " << numQuads
                << " is invalid, returned type "
                << QuadricClassify::QuadTypeNames[type]
                << "\n";
      imQuads++;
    }
    numQuads++;
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

template <typename ListTest, typename ListTrue>
bool validateResults(ListTest &inter, ListTrue &truth) {
  /* Note that because the same method is used to determine
   * if an intersection is in the list, these lists will be
   * the same length
   */
  if(inter->size() != truth->size()) {
    return false;
  }
  auto j = truth->begin();
  for(auto i = inter->begin();
      i != inter->end() || j != truth->end(); i++, j++) {
    if(i->q != j->q || i->intPos != j->intPos) {
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

using rngAlg = std::mt19937_64;

template <int dim, typename fptype>
using randLineGen = Geometry::Line<dim, fptype> (*)(
    rngAlg &rng, int minExp, int maxExp);

template <int dim, typename fptype, typename cmpAlg>
std::shared_ptr<std::list<cmpAlg>> __attribute__((noinline))
runTest(
    std::list<Geometry::Quadric<dim, fptype>> &quads,
    Geometry::Line<dim, fptype> &line, Timer::Timer &timer,
    double eps = std::numeric_limits<fptype>::infinity()) {
  timer.startTimer();
  auto inter =
      Geometry::sortIntersections<dim, fptype, cmpAlg>(
          line, quads, fptype(eps));
  timer.stopTimer();
  return inter;
}

struct Results {
  long long ns;
  int numIP;
  bool correct;
};
enum class Tests {
  RESULTANT = 0,
  SINGLE = 1,
  DOUBLE = 2,
  MP = 3,
};

template <int dim, typename fptype>
void intersectionTest(
    std::list<Geometry::Quadric<dim, fptype>> &quads,
    std::ostream &output, const int numTests,
    randLineGen<dim, fptype> rlgf, int minExp, int maxExp,
    bool fixedSeed) {
  /* First build a scene of quadrics.
   * Then generate random lines on a disk centered at the
   * intersection of the cylinders.
   * Then sort the intersections
   * Finally validate them with a higher precision sort
   */
  using Vf = Geometry::Vector<dim, fptype>;
  using Pf = Geometry::Point<dim, fptype>;
  using Lf = Geometry::Line<dim, fptype>;

  rngAlg engine;
  if(fixedSeed) {
    engine.seed(682185716);
  }

  const double eps = 1.52587890625e-5;

  constexpr const int numTestTypes = 4;

  Timer::Timer testTimes[numTestTypes];

  std::vector<std::array<Results, numTestTypes>> results;
  /* Run the tests */
  int t;
  for(t = 0; t < numTests && globals.run; t++) {
    Lf line = rlgf(engine, minExp, maxExp);
    /* Then sort the intersections */
    auto resultant = runTest<
        dim, fptype,
        Geometry::IntersectionResultant<dim, fptype>>(
        quads, line,
        testTimes[static_cast<int>(Tests::RESULTANT)], eps);

    auto sng = runTest<dim, fptype,
                       Geometry::IntersectionApproximate<
                           dim, fptype, float>>(
        quads, line,
        testTimes[static_cast<int>(Tests::SINGLE)], eps);

    auto dbl = runTest<dim, fptype,
                       Geometry::IntersectionApproximate<
                           dim, fptype, double>>(
        quads, line,
        testTimes[static_cast<int>(Tests::DOUBLE)], eps);

    auto mp = runTest<
        dim, fptype,
        Geometry::IntersectionIncreasedPrec<dim, fptype>>(
        quads, line, testTimes[static_cast<int>(Tests::MP)],
        eps);

    std::array<Results, numTestTypes> current;
    for(int i = 0; i < numTestTypes; i++) {
      current[i].ns = testTimes[i].instant_ns();
    }

    current[static_cast<int>(Tests::RESULTANT)].correct =
        validateResults(resultant, resultant);
    current[static_cast<int>(Tests::RESULTANT)].numIP =
        countIP(resultant);

    current[static_cast<int>(Tests::SINGLE)].correct =
        validateResults(sng, resultant);
    current[static_cast<int>(Tests::SINGLE)].numIP =
        countIP(sng);

    current[static_cast<int>(Tests::DOUBLE)].correct =
        validateResults(dbl, resultant);
    current[static_cast<int>(Tests::DOUBLE)].numIP =
        countIP(dbl);

    current[static_cast<int>(Tests::MP)].correct =
        validateResults(mp, resultant);
    current[static_cast<int>(Tests::MP)].numIP =
        countIP(mp);

    results.push_back(current);
  }
  /* Output all of the results */
  output << "Test #, Singles, Single Times (ns), Single "
            "Correct, Doubles, Double Times (ns), Double "
            "Correct, Increased Precs, MP Time (ns), MP "
            "Correct, Resultants, Resultant Time (ns)\n";
  int totIP[numTestTypes];
  int totIncorrect[numTestTypes];
  for(int i = 0; i < numTestTypes; i++) {
    totIP[i] = 0;
    totIncorrect[i] = 0;
  }
  for(int i = 0; i < t; i++) {
    output << i + 1 << ", ";
    for(int j = 0; j < numTestTypes; j++) {
      output << results[i][j].numIP << ", "
             << results[i][j].ns;
      if(j != numTestTypes - 1) {
        output << ", " << results[i][j].correct << ", ";
      }
      totIP[j] += results[i][j].numIP;
      totIncorrect[j] += 1 - results[i][j].correct;
    }
    output << "\n";
  }
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

template <typename fptype>
bool checkExp(fptype fpVal, int minExp, int maxExp) {
  if(fpVal == fptype(0.0)) {
    return true;
  }
  GenericFP::fpconvert<fptype> *fpbits =
      reinterpret_cast<GenericFP::fpconvert<fptype> *>(
          &fpVal);
  int exponent = fpbits->exponent;
  exponent -= fpbits->centralExp;
  if(exponent < minExp || exponent > maxExp) {
    return false;
  } else {
    return true;
  }
}

template <int dim, typename fptype>
Geometry::Line<dim, fptype> defRandLine(rngAlg &rng,
                                        int minExp,
                                        int maxExp) {
  constexpr const fptype minPos = 0, maxPos = 1;
  std::uniform_real_distribution<fptype> genPos(minPos,
                                                maxPos);
  std::uniform_real_distribution<fptype> genDir(-1.0, 1.0);
  /* First build the line */
  Geometry::Vector<dim, fptype> lineDir;
  Geometry::Vector<dim, fptype> lineInt;
  for(int i = 0; i < dim; i++) {
    do {
      lineInt.set(i, genPos(rng));
    } while(checkExp(lineInt.get(i), minExp, maxExp) ==
            false);
    do {
      lineDir.set(i, genDir(rng));
    } while(checkExp(lineDir.get(i), minExp, maxExp) ==
            false);
  }
  return Geometry::Line<dim, fptype>(
      Geometry::Point<dim, fptype>(lineInt), lineInt);
}

template <int dim, typename fptype>
Geometry::Line<dim, fptype> nestedEllRandLine(rngAlg &rng,
                                              int minExp,
                                              int maxExp) {
  constexpr const fptype minPos = 0.0, maxPos = 1.0;
  std::uniform_real_distribution<fptype> genPos(minPos,
                                                maxPos);
  Geometry::Vector<dim, fptype> lineInt;
  constexpr const fptype maxDeltaMag = 0.0625;
  std::uniform_real_distribution<fptype> genDelta(
      -maxDeltaMag, maxDeltaMag);
  Geometry::Vector<dim, fptype> delta;
  for(int i = 0; i < dim; i++) {
    do {
      lineInt.set(i, genPos(rng));
    } while(checkExp(lineInt.get(i), minExp, maxExp) ==
            false);
    do {
      delta.set(i, genDelta(rng));
    } while(checkExp(delta.get(i), minExp, maxExp) ==
            false);
  }
  /* Direct the line towards (1.0, 0.5, 0.5) + delta */
  const Geometry::Vector<dim, fptype> destination =
      Geometry::Vector<dim, fptype>({1.0, 0.5, 0.5}) +
      delta;
  Geometry::Vector<dim, fptype> lineDir =
      destination - lineInt;
  for(int i = 0; i < dim; i++) {
    if(checkExp(lineDir.get(i), minExp, maxExp) == false) {
      if(lineDir.get(i) > 1.0) {
        lineDir.set(i, 1.0);
      } else {
        lineDir.set(i, 0);
      }
    }
  }
  return Geometry::Line<dim, fptype>(
      Geometry::Point<dim, fptype>(lineInt), lineDir);
}

template <int dim, typename fptype>
Geometry::Line<dim, fptype> cylRandLine(rngAlg &rng,
                                        int minExp,
                                        int maxExp) {
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
  int minExp, maxExp;
  if(argc > 1) {
    quads = parseQuadrics<dim, fptype>(argv[1], &minExp,
                                       &maxExp);
  } else {
    quads = parseQuadrics<dim, fptype>("cylinders.csg",
                                       &minExp, &maxExp);
  }
  if(argc > 2) {
    outFName = argv[2];
  }
  if(argc > 3) {
    numTests = atoi(argv[3]);
  }
  randLineGen<dim, fptype> rlg = defRandLine<dim, fptype>;
  if(argc > 4) {
    int lineGenAlg = atoi(argv[4]);
    switch(lineGenAlg) {
      case 1:
        rlg = nestedEllRandLine<dim, fptype>;
        std::cout
            << "Using the nested spheres random lines\n";
        break;
      case 2:
        std::cout << "Using the cylinders random lines\n";
        rlg = cylRandLine<dim, fptype>;
        break;
      default:
        std::cout << "Using the default random lines\n";
        rlg = defRandLine<dim, fptype>;
    }
  }
  bool fixedSeed = false;
  if(argc > 5) {
    if(strcmp(argv[5], "fixed") == 0) {
      fixedSeed = true;
    }
  }
  std::ofstream results(outFName);
  signal(SIGINT, sigInt);
  intersectionTest(quads, results, numTests, rlg, minExp,
                   maxExp, fixedSeed);
  return 0;
}
