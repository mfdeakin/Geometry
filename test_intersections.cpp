
#include <array>
#include <list>

#include <gtest/gtest.h>

#include "quadrics.hpp"
#include "origin.hpp"
#include "accurate_intersections.hpp"
#include "accurate_math.hpp"

TEST(QLIntersect, LineIntersection) {
  using fptype = float;
  constexpr const int dim = 3;
  constexpr const int numCoeffs = (dim + 2) * (dim + 1) / 2;
  constexpr const int numQuadrics = 2;
  constexpr const fptype eps = 1e-3;
  using Q = Geometry::Quadric<dim, fptype>;
  using V = Geometry::Vector<dim, fptype>;
  using P = Geometry::Point<dim, fptype>;
  using L = Geometry::Line<dim, fptype>;
  P intercept(V({1.0, 0.0, -1.0}));
  L l(intercept, V({1.0, 1.0, 1.0}));
  fptype quadCoeffs[numQuadrics][Q::numCoeffs] = {
      {1.0, 1.0, 1.0, -3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
      {1.0, 1.0, 1.0, -3.00003, 0.0, 0.0, 0.0, 0.0, 0.0,
       0.0}};
  std::list<Q> quadrics;
  for(int i = 0; i < numQuadrics; i++) {
    Q q;
    for(int j = 0; j < numCoeffs; j++) {
      q.setCoeff(j, quadCoeffs[i][j]);
    }
    quadrics.push_back(q);
  }
  auto inter =
      Geometry::sortIntersections(l, quadrics, eps);
  for(auto intersects : *inter) {
    std::cout << "Intersections: (" << intersects.intPos
              << ", " << intersects.otherIntPos << ")\n";
  }
}
// Intersections: (-1, 1.00001)
// Intersections: (-1, 1)
// Intersections: (1, -1)
// Intersections: (1.00001, -1)
// Intersections: (-1.00005, 1.00005)
// Intersections: (-1, 1)
// Intersections: (1, -1)
// Intersections: (1.00005, -1.00005)
