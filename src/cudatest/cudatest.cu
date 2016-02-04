
#include "cudadef.h"

#include "line.hpp"
#include "quadric.hpp"
#include "accurate_intersections.hpp"

#include <gtest/gtest.h>

constexpr static const int dim = 3;
using fptype = float;

using O = Geometry::Origin<dim, fptype>;
using V = Geometry::Vector<dim, fptype>;
using P = Geometry::Point<dim, fptype>;
using L = Geometry::Line<dim, fptype>;
using Q = Geometry::Quadric<dim, fptype>;
using I = Geometry::Intersection<dim, fptype>;

__global__ void CompileTest() {
  V v1, v2;
  v1.dot(v2);
  Q q;
  for(int i = 0; i < dim + 1; i++)
    for(int j = 0; j < dim + 1; j++)
      q.setCoeff(i, j, (3.0 * i * i - i) / 1.5);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  int retVal = RUN_ALL_TESTS();
  constexpr const int numThreads = 1;
  constexpr const int numBlocks = 1;
  CompileTest<<<numBlocks, numThreads>>>();
  return retVal;
}
