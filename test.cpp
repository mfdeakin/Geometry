
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
  q.getCoeff(0, 0) = 2.0;
  q.getCoeff(1, 1) = 3.0;
  q.getCoeff(2, 2) = 5.0;
  q.getCoeff(3, 3) = 7.0;

  q.getCoeff(0, 1) = 11.0;
  q.getCoeff(0, 2) = 13.0;
  q.getCoeff(0, 3) = 17.0;

  q.getCoeff(1, 2) = 19.0;
  q.getCoeff(1, 3) = 23.0;

  q.getCoeff(2, 3) = 29.0;
  auto quadtype = AccurateMath::classifyQuadric(q);
  assert(quadtype < AccurateMath::QUADT_ERRORINVALID);
  printf("Quadric Type: %d, %s\n", quadtype,
         AccurateMath::QuadTypeNames[quadtype]);
  return 0;
}
