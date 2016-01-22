
#include "cudadef.h"

#include "line.hpp"
#include "quadrics.hpp"
#include "accurate_intersections.hpp"

#include <curand.h>
#include <curand_kernel.h>

constexpr static const int dim = 3;
using fptype = float;

using O = Geometry::Origin<dim, fptype>;
using V = Geometry::Vector<dim, fptype>;
using P = Geometry::Point<dim, fptype>;
using L = Geometry::Line<dim, fptype>;
using Q = Geometry::Quadric<dim, fptype>;
using I = Geometry::Intersection<dim, fptype>;

__global__ void setupKernel(curandState *state,
                            unsigned long long seed) {
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  curand_init(seed, idx, 0, &state[idx]);
}

template <L (*genLine)(curandState *), int numQuads,
          int numLines>
__global__ void computeIntersections(
    curandState *devRNGState,
    const Q::QuadricData *quadData,
    const I::IntersectionData *intResults) {
  /* The number of intersections needs to be the number of
   * quadrics times the number of lines to be computed */
  const int idx = threadIdx.x + blockDim.x * blockIdx.x;
  const int numThreads = blockDim.x * gridDim.x;
  const int linesPerThread = numLines / numThreads + 1;
  int lineIdx = linesPerThread * idx;
  curandState *rngState = &devRNGState[idx];
  Q quads[numQuads];
  for(int i = 0; i < numQuads; i++)
    for(int j = 0; j < dim; j++)
      quads[i].setCoeff(j, quadData[i].coeffs[j]);
  for(int i = 0; i < linesPerThread && lineIdx < numLines;
      i++, lineIdx++) {
    L line = genLine(rngState);
    for(int j = 0; j < numQuads; j++) {
      auto q = quads[j];
      auto intParams = q.calcLineDistToIntersect(line);
      I::IntersectionData i1;
      i1.quad = quadData[j];
      line.copyData(i1.line);
      i1.intPos = intParams[0];
      i1.intPos = intParams[1];
      I::IntersectionData i2;
      i2.quad = quadData[j];
      line.copyData(i2.line);
      i2.intPos = intParams[1];
      i2.intPos = intParams[0];
    }
  }
}

CUDA_DEVICE L uDistLine(curandState *rngState) {
  V intercept;
  V dir;
  fptype dSum = 0.0;
  for(int i = 0; i < dim; i++) {
    fptype p_i = curand_uniform_double(rngState);
    intercept.set(i, p_i);
    fptype d_i = curand_uniform_double(rngState);
    dir.set(i, d_i);
    dSum = fma(d_i, d_i, dSum);
  }
  if(dSum == 0.0) dir.set(0, 1.0);
  return L(P(intercept), dir);
}

int main(int argc, char **argv) {
  constexpr const int numBlocks = 128;
  constexpr const int numTPB = 32;
  constexpr const int numThreads = numBlocks * numTPB;
  curandState *curng;
  cudaError_t err =
      cudaMalloc(&curng, sizeof(*curng) * numThreads);
  int seed = 0;
  setupKernel<<<numBlocks, numTPB>>>(curng, seed);
  computeIntersections<uDistLine, 1,
                       1><<<numBlocks, numTPB>>>(
      curng, NULL, NULL);
  return 0;
}
