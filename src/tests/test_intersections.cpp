
#include <array>
#include <list>
#include <fstream>

#include <gtest/gtest.h>

#include "quadric.hpp"
#include "origin.hpp"
#include "accurate_intersections.hpp"
#include "accurate_math.hpp"
#include "timer.hpp"
#include "genericfp.hpp"

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
    EXPECT_EQ(intersects.intPos, expected[i]);
    i++;
  }
}
