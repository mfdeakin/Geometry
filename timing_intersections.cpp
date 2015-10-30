
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
  ifstream file(fname);
  using Qf = Geometry::Quadric<dim, fptype>;
  std::list<Qf> quads;
  return quads;
}

void intersectionTest(const int numTests) {
  /* TODO: Split this up into readable pieces */
  /* First build a scene of quadrics.
   * Then generate random lines on a disk centered at the
   * intersection of the cylinders.
   * Then sort the intersections
   * Finally validate them with a higher precision sort
   */
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
  /* The default scene: Two intersecting cylinders in a box
   */
  constexpr const fptype maxMag = 1000.0;
  constexpr const fptype radius = 50.0;
  constexpr const fptype cylXOff = 0.0;
  constexpr const fptype cylYOff = 50.0;
  constexpr const fptype cylZOff = 50.0;
  constexpr const fptype quadCoeffs[][Qf::numCoeffs] = {
      /* Exterior planes */
      {0, 0, 0, maxMag, 0, 0, 1.0, 0, 0, 0},
      {0, 0, 0, -maxMag, 0, 0, 1.0, 0, 0, 0},
      {0, 0, 0, maxMag, 0, 0, 0, 0, 1.0, 0},
      {0, 0, 0, -maxMag, 0, 0, 0, 0, 1.0, 0},
      {0, 0, 0, maxMag, 0, 0, 0, 0, 0, 1.0},
      {0, 0, 0, -maxMag, 0, 0, 0, 0, 0, 1.0},
      /* Cylinder 1 */
      {1.0, 1.0, 0, -radius * radius + cylXOff * cylXOff +
                        cylYOff * cylYOff,
       0, 0, -2 * cylXOff, 0, -2 * cylYOff, 0},
      /* Cylinder 2 */
      {1.0, 0, 1.0, -radius * radius + cylXOff * cylXOff +
                        cylZOff * cylZOff,
       0, 0, -2 * cylXOff, 0, 0, -2 * cylZOff},
      /* Intersection Plane 1 */
      {0, 0, 0, -cylXOff + cylZOff, 0, 0, 0, 0, 1, -1},
      /* Intersection Plane 1 */
      {0, 0, 0, cylXOff - cylZOff, 0, 0, 0, 0, -1, 1}};
  constexpr const int numQuads =
      sizeof(quadCoeffs) / sizeof(quadCoeffs[0]);
  std::list<Qf> quads;
  std::list<Qm> truthQuads;
  /* Generate the quadrics */
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
      std::max(std::fabs(cylXOff), std::fabs(cylYOff)),
      std::fabs(cylZOff));
  std::uniform_real_distribution<fptype> genDist(
      -maxMag + maxDisp, maxMag - maxDisp);
  std::ofstream results("results");
  Timer::Timer fp_time, mp_time;
  /* Run the tests */
  for(int t = 0; t < numTests; t++) {
    /* First build the line */
    fptype theta = genTheta(engine);
    fptype dist = genDist(engine);
    fptype xDisp = dist * std::sin(theta);
    fptype yDisp = dist * std::cos(theta);
    fptype zDisp = yDisp;
    Vf lineDir({-xDisp, -yDisp, -zDisp});
    Vf lineInt({xDisp - cylXOff, yDisp - cylYOff,
                zDisp - cylZOff});
    Pf intercept(lineInt);
    Lf line(intercept, lineDir);
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
    auto j = truth->begin();
    results << "Test " << t << "\n";
    results << "FP Time:\n" << fp_time << "\n";
    results << "MP Time:\n" << mp_time << "\n";
    /* Now validate them */
    for(auto i = inter->begin();
        i != inter->end() || j != truth->end();) {
      bool printQ = false;
      if(i != inter->end()) {
        if(j != truth->end()) {
          if(i->q != j->q) {
            results << "Incorrect Result\n";
            results << "Expected: " << std::setprecision(20)
                    << j->intPos
                    << "; Got: " << std::setprecision(20)
                    << i->intPos
                    << "\nDelta: " << std::setprecision(20)
                    << j->intPos - i->intPos << "\n";
            printQ = true;
          }
        } else {
          printQ = true;
        }
      }
      /* And output the result */
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
  intersectionTest(1e5);
  return 0;
}
