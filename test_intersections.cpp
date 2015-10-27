
#include <array>
#include <list>
#include <fstream>

#include <gtest/gtest.h>

#include "quadrics.hpp"
#include "origin.hpp"
#include "accurate_intersections.hpp"
#include "accurate_math.hpp"
#include "genericfp.hpp"

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
  constexpr const int truthPrec = 48 * machPrec;
  mpfr::mpreal::set_default_prec(truthPrec);
  constexpr const fptype maxMag = 100000.0;
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
      {
          1.0, 1.0, 0.0, -100.0, 0.0, 0.0, 2 * cyl1XOff,
          2 * cyl1YOff, 0.0,
          cyl1XOff * cyl1XOff + cyl1YOff * cyl1YOff,
      },
      /* Cylinder 2 */
      {
          1.0, 0.0, 1.0, -100.0, 0.0, 0.0, 2 * cyl2XOff,
          0.0, 2 * cyl2ZOff,
          cyl2XOff * cyl2XOff + cyl2ZOff * cyl2ZOff,
      },
      /* Intersection Plane 1 */
      {0, 0, 0, cyl1XOff + cyl2ZOff, 0, 0, 1, 0, 0, 1},
      /* Intersection Plane 1 */
      {0, 0, 0, -cyl1XOff - cyl2ZOff, 0, 0, -1, 0, 0, -1}};
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
  std::uniform_real_distribution<fptype> rgenf(-maxMag,
                                               maxMag);
  std::ofstream results("results");
  for(int t = 0; t < numTests; t++) {
    Vf lineDir;
    Vf lineInt;
    for(int i = 0; i < dim; i++) {
      lineDir.set(i, rgenf(engine));
      lineInt.set(i, rgenf(engine));
    }
    Pf intercept(lineInt);
    Lf line(intercept, lineDir);
    Lm truthLine(line);
    constexpr const fptype eps = 1e-3;
    auto inter =
        Geometry::sortIntersections(line, quads, eps);
    auto truth =
        Geometry::sortIntersections<dim, mpfr::mpreal>(
            truthLine, truthQuads, eps);
    auto j = truth->begin();
    results << "Test " << t << "\n";
    for(auto i = inter->begin();
        i != inter->end() || j != truth->end();) {
      if(i != inter->end()) {
        if(j != truth->end()) {
          if(i->q == j->q)
            results << "Correct Result\n";
          else
            results << "Incorrect Result\n";
        }
        results << "Estimated: " << i->q
                << "\n";
        i++;
      }
      if(j != truth->end()) {
        results << "Correct: " << j->q << "\n";
        j++;
      }
    }
    results << "\n\n";
  }
}

TEST(QLIntersect, LineIntersection) {
  using fptype = float;
  constexpr const int dim = 3;
  constexpr const fptype eps = 1e-3;
  using Q = Geometry::Quadric<dim, fptype>;
  using V = Geometry::Vector<dim, fptype>;
  using P = Geometry::Point<dim, fptype>;
  using L = Geometry::Line<dim, fptype>;
  P intercept(V({1.0, 0.0, -1.0}));
  L l(intercept, V({1.0, 1.0, 1.0}));
  fptype quadCoeffs[][Q::numCoeffs] = {
      {1.0, 1.0, 1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0, 1.0, 1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0, 0.0, 0.25, -3.0, -2.0, 0.0, 0.0, 0.0, 2.0, 0.0},
      {1.0, 0.0, 0.25, -3.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0},
      {1.0, 1.0, 1.0, -3.00003, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0}};
  constexpr const fptype expected[] = {
      -3.4055080413818359375,  -1.0000149011611938477,
      -0.99999994039535522461, -0.99999994039535522461,
      0.47434464097023010254,  0.99999988079071044922,
      0.99999988079071044922,  1.0000150203704833984,
  };
  constexpr const int numQuadrics =
      sizeof(quadCoeffs) / sizeof(quadCoeffs[0]);
  std::list<Q> quadrics;
  for(int i = 0; i < numQuadrics; i++) {
    Q q;
    for(int j = 0; j < Q::numCoeffs; j++) {
      q.setCoeff(j, quadCoeffs[i][j]);
    }
    quadrics.push_back(q);
  }
  auto inter =
      Geometry::sortIntersections(l, quadrics, eps);
  int i = 0;
  for(auto intersects : *inter) {
    std::cout << "Intersections: (" << std::setprecision(20)
              << intersects.intPos << ", "
              << std::setprecision(20)
              << intersects.otherIntPos
              << "), vs expected: "
              << intersects.intPos - expected[i] << "\n";
    EXPECT_EQ(intersects.intPos, expected[i]);
    i++;
  }
}
