
#include <gtest/gtest.h>

#include <iostream>
#include <array>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "vector.hpp"

TEST(Line, ShiftOrigin) {
  static constexpr const int dim = 3;
  using fptype = float;
  const fptype eps = 1e-4;
  struct TestCase {
    std::array<fptype, dim> newOrigin;
    std::array<fptype, dim> origOffset;
    std::array<fptype, dim> dir;
    std::array<fptype, dim> expectedOffset;
  };
  struct TestCase tests[] = {
      {{{0.0, 0.0, 0.0}},
       {{0.0, 0.0, 0.0}},
       {{0.0, 0.0, 1.0}},
       {{0.0, 0.0, 0.0}}},

      {{{0.0, 0.0, 0.0}},
       {{1.0, 0.0, 0.0}},
       {{0.0, 1.0, 1.0}},
       {{1.0, 0.0, 0.0}}},

      {{{2.0, 0.0, 0.0}},
       {{1.0, 0.0, 0.0}},
       {{0.0, 1.0, 1.0}},
       {{-1.0, 0.0, 0.0}}},

      {{{1.0, 0.0, 0.0}},
       {{1.0, 0.0, 0.0}},
       {{0.0, 1.0, 1.0}},
       {{0.0, 0.0, 0.0}}},

      {{{0.0, 0.0, 0.0}},
       {{1.0, 1.0, 0.0}},
       {{0.0, 1.0, 1.0}},
       {{1.0, 0.5, -0.5}}},
  };
  Geometry::Origin<dim, fptype> defOrigin;
  for(auto t : tests) {
    Geometry::Point<dim, fptype> intersect(
        defOrigin,
        Geometry::Vector<dim, fptype>(t.origOffset));
    Geometry::Vector<dim, fptype> dir(
        Geometry::Vector<dim, fptype>(t.dir).normalize());
    Geometry::Line<dim, fptype> line(intersect, dir);
    Geometry::Origin<dim, fptype> o(t.newOrigin);
    line.shiftOrigin(o);
    Geometry::Vector<dim, fptype> newOff(
        line.getIntercept().getOffset());
    for(unsigned i = 0; i < dim; i++) {
      EXPECT_NEAR(newOff(i), t.expectedOffset[i], eps);
    }
  }
}
