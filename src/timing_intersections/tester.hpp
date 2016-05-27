
#ifndef _TESTER_HPP_
#define _TESTER_HPP_

#include "accurate_intersections_approx.hpp"
#include "accurate_intersections_incprec.hpp"
#include "accurate_intersections_resultant.hpp"
#include "timer.hpp"

#include <iostream>

struct TestResults {
  long long time_ns;
  int numAccurateComp;
  bool correct;
};

void printResult(const TestResults &result,
                 std::ostream &output) {
  output << result.numAccurateComp << ", " << result.time_ns
         << ", " << result.correct;
}

template <int dim, typename fptype, typename Comparator>
class Tester {
 public:
  Tester(
      int numTests,
      std::list<Geometry::Quadric<dim, fptype>> *quads,
      double eps = std::numeric_limits<fptype>::infinity())
      : timer(),
        quads(quads),
        prevResults(NULL),
        eps(eps),
        numTests(numTests),
        curTest(0),
        totalTime(0),
        totalIncorrect(0),
        results(new TestResults[numTests]) {}

  virtual ~Tester() { delete[] results; }

  void runTest(const Geometry::Line<dim, fptype> &line) {
    timer.startTimer();
    prevResults = Geometry::sortIntersections<dim, fptype,
                                              Comparator>(
        line, *quads, fptype(eps));
    timer.stopTimer();
  }

  template <typename TruthComparator>
  void updateResults(
      std::shared_ptr<std::list<TruthComparator>> truth) {
    assert(curTest < numTests);
    results[curTest].time_ns = timer.instant_ns();
    results[curTest].numAccurateComp = countCompares();
    results[curTest].correct = validateResults(truth);
    totalTime += results[curTest].time_ns;
    totalIncorrect += !results[curTest].correct;
    curTest++;
  }

  int countCompares() {
    int compares = 0;
    for(auto intersect : *prevResults) {
      compares += intersect.incPrecCount();
    }
    return compares;
  }

  template <typename TruthComparator>
  bool validateResults(
      std::shared_ptr<std::list<TruthComparator>> truth) {
    if(prevResults->size() != truth->size()) {
      return false;
    }
    auto t = truth->begin();
    for(auto i = prevResults->begin();
        i != prevResults->end(); i++, t++) {
      if(i->q != t->q || (i->intPos < i->otherIntPos &&
                          t->intPos >= t->otherIntPos) ||
         (i->intPos == i->otherIntPos &&
          t->intPos != t->otherIntPos) ||
         (i->intPos > i->otherIntPos &&
          t->intPos <= t->otherIntPos)) {
        return false;
      }
      if(i->intPos != t->intPos) {
        std::cout << "Something went terribly wrong\n";
        return false;
      }
    }
    return true;
  }

  const TestResults &getResult(int testNum) {
    assert(testNum < curTest);
    return results[testNum];
  }

  std::shared_ptr<std::list<Comparator>> getResults() {
    return prevResults;
  }

  long long getTotalTime_ns() { return totalTime; }

  int getTotalIncorrect() { return totalIncorrect; }

 private:
  Timer::Timer timer;
  std::list<Geometry::Quadric<dim, fptype>> *quads;
  std::shared_ptr<std::list<Comparator>> prevResults;
  double eps;
  int numTests;
  int curTest;
  long long totalTime;
  int totalIncorrect;
  TestResults *results;
};

#endif
