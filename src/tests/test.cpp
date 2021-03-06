
#include <cstddef>

#include "origin.hpp"
#include "quadric.hpp"
#include "line.hpp"
#include "bsgtree.hpp"
#include "accurate_math.hpp"
#include "polynomial.hpp"
#include "array.hpp"

#include "mpreal.hpp"

#include <iostream>

#include <list>
#include <cmath>

#include <gtest/gtest.h>

int main(int argc, char **argv) {
  static constexpr const int dim = 3;
  using fptype = float;
  Geometry::Origin<dim, fptype> o;
  Geometry::Vector<dim, fptype> ptOffset;
  Geometry::Point<dim, fptype> pt(o, ptOffset);
  Geometry::Vector<dim, fptype> dir(
      Array<fptype, dim>({1.0, 1.0, 1.0}));
  Geometry::Line<dim, fptype> line(pt, dir);
  Geometry::Quadric<dim, fptype> q(o);
  Geometry::Polynomial<dim, fptype> p1, p2;
  Geometry::Polynomial<2 * dim, fptype> prod =
      p1.product(p2);
  mpfr::mpreal::set_default_prec(72);
  mpfr::mpreal mt;
  Geometry::Quadric<dim, mpfr::mpreal> qmp(o);
  QuadricClassify::classifyQuadric(qmp);
  std::list<Geometry::Quadric<dim, mpfr::mpreal>> qmplist;
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
