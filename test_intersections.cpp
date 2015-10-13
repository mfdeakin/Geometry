
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
  constexpr const int numQuadrics = 10;
  constexpr const fptype eps = 1e-3;
  using Q = Geometry::Quadric<dim, fptype>;
  using V = Geometry::Vector<dim, fptype>;
  using O = Geometry::Origin<dim, fptype>;
  using P = Geometry::Point<dim, fptype>;
  using L = Geometry::Line<dim, fptype>;
  O o;
  P intercept(o, V({1.0, 0.0, -1.0}));
  L l(intercept, V({1.0, 1.0, 1.0}));
  std::list<Q> quadrics;
  for(int i = 0; i < numQuadrics; i++) {
    Q q(o);
    for(int j = 0; j < numCoeffs; j++) {
      q.coeff(j) = 0.0;
    }
    quadrics.push_back(q);
  }
  auto inter =
      Geometry::sortIntersections(l, quadrics, eps);
  for(auto intersects : *inter) {
  }
}
