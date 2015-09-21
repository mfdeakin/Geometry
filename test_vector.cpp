
#include <gtest/gtest.h>

#include <iostream>
#include <array>

#include "geometry.hpp"
#include "vector.hpp"

TEST(Vector, OrthogonalBasis) {
  constexpr const int dim = 3;
  typedef float fptype;
  const fptype eps = 1e-7;
  std::array<fptype, dim> tests[] = {{{1.0, 0.0, 0.0}},
                                     {{0.0, 1.0, 0.0}},
                                     {{0.0, 0.0, 1.0}},
                                     {{1.0, 1.0, 0.0}},
                                     {{1.0, 0.0, 1.0}},
                                     {{0.0, 1.0, 1.0}},
                                     {{1.0, 1.0, 1.0}}};
  for(auto t : tests) {
    Geometry::Vector<dim, fptype> v(t);
    auto orthogs = v.calcOrthogonals();
    for(unsigned i = 0; i < orthogs.size(); i++) {
      auto o1 = orthogs[i];
      fptype dp = o1.dot(v);
      EXPECT_NEAR(dp, 0.0, eps);
      for(unsigned j = i + 1; j < orthogs.size(); j++) {
        auto o2 = orthogs[j];
        dp = o1.dot(o2);
        EXPECT_NEAR(dp, 0.0, eps);
      }
    }
  }
}