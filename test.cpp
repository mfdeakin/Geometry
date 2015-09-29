
#include "origin.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "bsgtree.hpp"
#include "accurate_math.hpp"

#include <iostream>

#include <array>
#include <cmath>

#include <gtest/gtest.h>

int main(int argc, char **argv) {
  static constexpr const int dim = 3;
  using fptype = float;
  Geometry::Origin<dim, fptype> o;
  Geometry::Vector<dim, fptype> ptOffset;
  Geometry::Point<dim, fptype> pt(o, ptOffset);
  Geometry::Vector<dim, fptype> dir(
      std::array<fptype, dim>({{1.0, 1.0, 1.0}}));
  Geometry::Line<dim, fptype> line(pt, dir);
  Geometry::Quadric<dim, fptype> q(o);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
