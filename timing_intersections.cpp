
#include "geometry.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "accurate_intersections.hpp"
#include "timer.hpp"

#include <list>
#include <fstream>

void intersectionTest(const int numTests) {
  constexpr const int dim = 3;
  using fptype = float;
  using Qf = Geometry::Quadric<dim, fptype>;
  using Vf = Geometry::Vector<dim, fptype>;
  using Pf = Geometry::Point<dim, fptype>;
  using Lf = Geometry::Line<dim, fptype>;
  using Qm = Geometry::Quadric<dim, mpfr::mpreal>;
  using Lm = Geometry::Line<dim, mpfr::mpreal>;
  constexpr const int machPrec =
      GenericFP::fpconvert<fptype>::precision;
  constexpr const int truthPrec = 96 * machPrec;
  mpfr::mpreal::set_default_prec(truthPrec);
  constexpr const fptype maxMag = 100000.0;
  constexpr const fptype radius = 50.0;
  constexpr const fptype cyl1XOff = 50.0;
  constexpr const fptype cyl1YOff = 50.0;
  constexpr const fptype cyl2XOff = cyl1XOff;
  constexpr const fptype cyl2ZOff = 50.0;
  constexpr const fptype quadCoeffs[][Qf::numCoeffs] = {
      /* Exterior planes */
      {0, 0, 0, maxMag, 0, 0, 1.0, 0, 0, 0},
      {0, 0, 0, -maxMag, 0, 0, 1.0, 0, 0, 0},
      {0, 0, 0, maxMag, 0, 0, 0, 0, 1.0, 0},
      {0, 0, 0, -maxMag, 0, 0, 0, 0, 1.0, 0},
      {0, 0, 0, maxMag, 0, 0, 0, 0, 0, 1.0},
      {0, 0, 0, -maxMag, 0, 0, 0, 0, 0, 1.0},
      /* Cylinder 1 */
      {1.0, 1.0, 0, -radius * radius + cyl1XOff * cyl1XOff +
                        cyl1YOff * cyl1YOff,
       0, 0, -2 * cyl1XOff, 0, -2 * cyl1YOff, 0},
      /* Cylinder 2 */
      {1.0, 0, 1.0, -radius * radius + cyl2XOff * cyl2XOff +
                        cyl2ZOff * cyl2ZOff,
       0, 0, -2 * cyl2XOff, 0, 0, -2 * cyl2ZOff},
      /* Intersection Plane 1 */
      {0, 0, 0, cyl1XOff + cyl2ZOff, 0, 0, 0, 0, 1, -1},
      /* Intersection Plane 1 */
      {0, 0, 0, -cyl1XOff - cyl2ZOff, 0, 0, 0, 0, -1, 1}};
  constexpr const int numQuads =
      sizeof(quadCoeffs) / sizeof(quadCoeffs[0]);
  std::list<Qf> quads;
  std::list<Qm> truthQuads;
  for(int i = 0; i < numQuads; i++) {
    Qf quadFP;
    Qm quadMP;
    for(int j = 0; j < Qf::numCoeffs; j++) {
      quadFP.setCoeff(j, quadCoeffs[i][j]);
      quadMP.setCoeff(j, quadCoeffs[i][j]);
    }
    quads.push_back(quadFP);
    truthQuads.push_back(quadMP);
  }
  std::random_device rd;
  std::mt19937_64 engine(rd());
  constexpr const fptype maxTheta = M_PI;
  std::uniform_real_distribution<fptype> genTheta(-maxTheta,
                                                  maxTheta);
  constexpr fptype maxDisp = std::max(
      std::max(std::fabs(cyl1XOff), std::fabs(cyl1YOff)),
      std::fabs(cyl2ZOff));
  std::uniform_real_distribution<fptype> genDist(
      -maxMag - maxDisp, maxMag - maxDisp);
  std::ofstream results("results");
  Timer::Timer fp_time, mp_time;
  for(int t = 0; t < numTests; t++) {
    fptype theta = genTheta(engine);
    fptype dist = genDist(engine);
    fptype xDisp = dist * std::sin(theta);
    fptype yDisp = dist * std::cos(theta);
    fptype zDisp = yDisp;
    Vf lineDir({xDisp, yDisp, zDisp});
    Vf lineInt({xDisp - cyl1XOff, yDisp - cyl1YOff,
                zDisp - cyl2ZOff});
    Pf intercept(lineInt);
    Lf line(intercept, lineDir);
    Lm truthLine(line);
    constexpr const fptype eps = 1e-3;
    fp_time.startTimer();
    auto inter =
        Geometry::sortIntersections(line, quads, eps);
    fp_time.stopTimer();
    mp_time.startTimer();
    auto truth =
        Geometry::sortIntersections<dim, mpfr::mpreal>(
            truthLine, truthQuads, eps);
    mp_time.stopTimer();
    auto j = truth->begin();
    results << "Test " << t << "\n";
    results << "FP Time:\n" << fp_time << "\n";
    results << "MP Time:\n" << mp_time << "\n";
    for(auto i = inter->begin();
        i != inter->end() || j != truth->end();) {
      bool printQ = false;
      if(i != inter->end()) {
        if(j != truth->end()) {
          if(i->q != j->q) {
            results << "Incorrect Result\n";
            printQ = true;
          }
        } else {
          printQ = true;
        }
      }
      if(printQ) {
        if(i != inter->end()) {
          results << "Estimated: " << i->q << "\n";
          i++;
        } else {
          results << "Computed intersections ended "
                     "prematurely\n";
        }
        if(j != truth->end()) {
          results << "Correct: " << j->q << "\n";
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
}

int main(int argc, char **argv) {
  intersectionTest(5e6);
  return 0;
}
