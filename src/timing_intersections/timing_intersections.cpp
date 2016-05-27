
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

#include "tester.hpp"

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

using rngAlg = std::mt19937_64;

template <int dim, typename fptype>
using randLineGen = Geometry::Line<dim, fptype> (*)(
    rngAlg &rng, int minExp, int maxExp);

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

  Tester<dim, fptype,
         Geometry::IntersectionResultant<dim, fptype>>
      resultantTest(numTests, &quads, eps);

  Tester<dim, fptype, Geometry::IntersectionApproximate<
                          dim, fptype, float>>
      singleTest(numTests, &quads, eps);

  Tester<dim, fptype, Geometry::IntersectionApproximate<
                          dim, fptype, float>>
      doubleTest(numTests, &quads, eps);

  Tester<dim, fptype,
         Geometry::IntersectionIncreasedPrec<dim, fptype>>
      mpTest(numTests, &quads, eps);

  std::vector<Lf> lines;
  /* Run the tests */
  int t;
  for(t = 0; t < numTests && globals.run; t++) {
    Lf line = rlgf(engine, minExp, maxExp);
    /* Then sort the intersections */
    resultantTest.runTest(line);
    singleTest.runTest(line);
    doubleTest.runTest(line);
    mpTest.runTest(line);
    lines.push_back(line);
    singleTest.updateResults(resultantTest.getResults());
    doubleTest.updateResults(resultantTest.getResults());
    mpTest.updateResults(resultantTest.getResults());
    resultantTest.updateResults(resultantTest.getResults());
  }
  /* Output the aggregate of the results */
  output << "Singles Incorrect: "
         << singleTest.getTotalIncorrect() << "\n"
         << "Single Total Time (ns): "
         << singleTest.getTotalTime_ns() << "\n\n";
  output << "Doubles Incorrect: "
         << doubleTest.getTotalIncorrect() << "\n"
         << "Double Total Time (ns): "
         << doubleTest.getTotalTime_ns() << "\n\n";
  output << "Increased Precision Incorrect: "
         << mpTest.getTotalIncorrect() << "\n"
         << "Increased Precision Total Time (ns): "
         << mpTest.getTotalTime_ns() << "\n\n";
  output << "Resultant Total Time (ns): "
         << resultantTest.getTotalTime_ns() << "\n\n\n";
  /* Output all of the results */
  output << "Test #, Singles, Single Times (ns), Single "
            "Correct, Doubles, Double Times (ns), Double "
            "Correct, Increased Precs, MP Time (ns), MP "
            "Correct, Resultants, Resultant Time (ns), "
            "Line Dir X, Line Dir Y, Line Dir Z, "
            "Line Int X, Line Int Y, Line Int Z\n";
  for(int i = 0; i < t; i++) {
    output << i + 1 << ", ";
    printResult(singleTest.getResult(i), output);
    output << ", ";
    printResult(doubleTest.getResult(i), output);
    output << ", ";
    printResult(mpTest.getResult(i), output);
    output << ", ";
    printResult(resultantTest.getResult(i), output);
    output << ", ";
    output << lines[i].getDirection().get(0);
    output << ", ";
    output << lines[i].getDirection().get(1);
    output << ", ";
    output << lines[i].getDirection().get(2);
    output << ", ";
    output << lines[i].getIntercept().getOffset().get(0);
    output << ", ";
    output << lines[i].getIntercept().getOffset().get(1);
    output << ", ";
    output << lines[i].getIntercept().getOffset().get(2);
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
      union {
        fptype tmp;
        GenericFP::fpconvert<fptype> bitmanipulate;
      };
      tmp = genPos(rng);
      bitmanipulate.mantissa &= ~((1 << 12) - 1);
      lineInt.set(i, tmp);
    } while(checkExp(lineInt.get(i), minExp, maxExp) ==
            false);
    do {
      union {
        fptype tmp;
        GenericFP::fpconvert<fptype> bitmanipulate;
      };
      tmp = genDir(rng);
      bitmanipulate.mantissa &= ~((1 << 12) - 1);
      lineDir.set(i, tmp);
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
  constexpr const fptype maxDeltaMag = 0.0078125;
  std::uniform_real_distribution<fptype> genDelta(
      -maxDeltaMag, maxDeltaMag);
  Geometry::Vector<dim, fptype> delta;
  for(int i = 0; i < dim; i++) {
    do {
      union {
        fptype tmp;
        GenericFP::fpconvert<fptype> bitmanipulate;
      };
      tmp = genPos(rng);
      bitmanipulate.mantissa &= ~((1 << 12) - 1);
      lineInt.set(i, tmp);
    } while(checkExp(lineInt.get(i), minExp, maxExp) ==
            false);
    delta.set(i, genDelta(rng));
  }
  /* Direct the line towards (1.0, 0.5, 0.5) + delta */
  const Geometry::Vector<dim, fptype> destination =
      Geometry::Vector<dim, fptype>({1.0, 0.5, 0.5}) +
      delta;
  Geometry::Vector<dim, fptype> lineDir =
      destination - lineInt;
  /* Fix the exponent and the mantissa to fit within our
   * bounds for correctness
   */
  int maximum = std::numeric_limits<int>::min();
  for(int i = 0; i < dim; i++) {
    int exp = 0;
    MathFuncs::MathFuncs<fptype>::frexp(lineDir.get(i),
                                        &exp);
    if(exp > maximum) {
      maximum = exp;
    }
  }
  fptype multiplier = MathFuncs::MathFuncs<fptype>::ldexp(
      1.0, maxExp - maximum);
  for(int i = 0; i < dim; i++) {
    union {
      fptype tmp;
      GenericFP::fpconvert<fptype> bitmanipulate;
    };

    tmp = lineDir.get(i) * multiplier;
    bitmanipulate.mantissa &= ~((1 << 12) - 1);
    if(MathFuncs::MathFuncs<fptype>::fabs(tmp) <
       MathFuncs::MathFuncs<fptype>::ldexp(1.0, minExp)) {
      if(MathFuncs::MathFuncs<fptype>::fabs(tmp) >
         MathFuncs::MathFuncs<fptype>::ldexp(1.0,
                                             minExp - 1)) {
        tmp = MathFuncs::MathFuncs<fptype>::ldexp(1.0,
                                                  minExp);
      } else {
        tmp = 0.0;
      }
    }
    lineDir.set(i, tmp);
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
  union {
    fptype tmp;
    GenericFP::fpconvert<fptype> bitmanipulate;
  };
  tmp = genPos(rng);
  bitmanipulate.mantissa &= ~((1 << 12) - 1);
  lineInt.set(1, tmp);
  for(int i = 0; i < dim; i++) {
    lineInt.set(i, lineInt.get(1));
    union {
      fptype tmp;
      GenericFP::fpconvert<fptype> bitmanipulate;
    };
    tmp = genDir(rng);
    bitmanipulate.mantissa &= ~((1 << 12) - 1);
    lineDir.set(i, tmp);
  }
  tmp = genPos(rng);
  bitmanipulate.mantissa &= ~((1 << 12) - 1);
  lineInt.set(0, tmp);
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
