
#include <gtest/gtest.h>

#include <iostream>

#include "geometry.hpp"
#include "origin.hpp"
#include "point.hpp"
#include "vector.hpp"
#include "array.hpp"

TEST(Line, ShiftOrigin) {
  static constexpr const int dim = 3;
  using fptype = float;
  const fptype eps = 1e-4;
  struct TestCase {
    Array<fptype, dim> newOrigin;
    Array<fptype, dim> origOffset;
    Array<fptype, dim> dir;
    Array<fptype, dim> expectedOffset;
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

TEST(Line, CalcPointAtDistance) {
  static constexpr const int dim = 3;
  using fptype = float;
  const fptype eps = 1e-6;
  struct TestCase {
    fptype distance;
    Array<fptype, dim> offset;
    Array<fptype, dim> dir;
    Array<fptype, dim> expected;
  };
  static const fptype pi =
      (fptype)std::atan((fptype)1) * 4.0;
  struct TestCase tests[] = {
      {0.0,
       {{0.0, 0.0, 0.0}},
       {{0.0, 0.0, 1.0}},
       {{0.0, 0.0, 0.0}}},
      {1.0,
       {{0.0, 0.0, 0.0}},
       {{0.0, 0.0, 1.0}},
       {{0.0, 0.0, 1.0}}},
      {std::sqrt((fptype)2),
       {{0.0, 0.0, 0.0}},
       {{1.0, 0.0, 1.0}},
       {{1.0, 0.0, 1.0}}},
      {3 * std::sqrt((fptype)2),
       {{0.0, 0.0, 0.0}},
       {{1.0, 0.0, 1.0}},
       {{3.0, 0.0, 3.0}}},
      {pi * std::sqrt((fptype)3),
       {{0.0, 0.0, 0.0}},
       {{1.0, 1.0, 1.0}},
       {{pi, pi, pi}}},
      {-1.0,
       {{0.0, 0.0, 0.0}},
       {{0.0, 0.0, 1.0}},
       {{0.0, 0.0, -1.0}}},
      {-std::sqrt((fptype)2),
       {{0.0, 0.0, 0.0}},
       {{1.0, 0.0, 1.0}},
       {{-1.0, 0.0, -1.0}}},
      {-3 * std::sqrt((fptype)2),
       {{0.0, 0.0, 0.0}},
       {{1.0, 0.0, 1.0}},
       {{-3.0, 0.0, -3.0}}},
      {-pi * std::sqrt((fptype)3),
       {{0.0, 0.0, 0.0}},
       {{1.0, 1.0, 1.0}},
       {{-pi, -pi, -pi}}},
  };
  Geometry::Origin<dim, fptype> defOrigin;
  for(auto t : tests) {
    Geometry::Vector<dim, fptype> offset(t.offset);
    Geometry::Point<dim, fptype> intersect(defOrigin,
                                           offset);
    Geometry::Vector<dim, fptype> dir(t.dir);
    Geometry::Line<dim, fptype> line(intersect, dir.normalize());
    Geometry::Point<dim, fptype> pos =
        line.getPosAtDist(t.distance);
    Geometry::Vector<dim, fptype> off = pos.getOffset();
    for(unsigned i = 0; i < dim; i++) {
      EXPECT_NEAR(off(i), t.expected[i], eps);
    }
  }
}
