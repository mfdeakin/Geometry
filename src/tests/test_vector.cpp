
#include <gtest/gtest.h>

#include <iostream>

#include "geometry.hpp"
#include "vector.hpp"
#include "array.hpp"

TEST(Vector, OrthogonalBasis) {
  constexpr const int dim = 3;
  using fptype = float;
  const fptype eps = 1e-4;
  Array<fptype, dim> tests[] = {
      {{1.0, 0.0, 0.0}},     {{0.0, 1.0, 0.0}},
      {{0.0, 0.0, 1.0}},     {{1.0, 1.0, 0.0}},
      {{1.0, 0.0, 1.0}},     {{0.0, 1.0, 1.0}},
      {{1.0, 1.0, 1.0}},     {{1.0, 1.0, 2.0}},
      {{1.0, 1.0e3, 1.0e-3}}};
  for(auto t : tests) {
    Geometry::Vector<dim, fptype> v(t);
    for(unsigned i = 0; i < dim - 1; i++) {
      auto o1 = v.getOrthogonal(i);
      fptype dp = o1.dot(v);
      EXPECT_NEAR(dp, 0.0, eps);
      for(unsigned j = i + 1; j < dim - 1; j++) {
        auto o2 = v.getOrthogonal(j);
        dp = o1.dot(o2);
        EXPECT_NEAR(dp, 0.0, eps);
      }
    }
  }
}
