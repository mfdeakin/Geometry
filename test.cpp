
#include "origin.hpp"
#include "quadrics.hpp"
#include "line.hpp"
#include "bsgtree.hpp"
#include "accurate_math.hpp"
#include <array>
#include <math.h>
#include <stdio.h>
#include <gtest/gtest.h>

int main(int argc, char **argv) {
  constexpr int dim = 3;
  typedef float fptype;
  Geometry::Origin<dim, fptype> o;
  Geometry::Vector<dim, fptype> ptOffset;
  Geometry::Point<dim, fptype> pt(o, ptOffset);
  Geometry::Vector<dim, fptype> dir;
  Geometry::Line<dim, fptype> line(pt, dir);
  Geometry::Quadric<dim, fptype> q(o);
  auto orthogs = dir.calcOrthogonals();
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
