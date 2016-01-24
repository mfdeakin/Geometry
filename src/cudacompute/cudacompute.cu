
#include "cudadef.h"

#include "line.hpp"
#include "quadrics.hpp"
#include "accurate_intersections.hpp"

#include <curand.h>
#include <curand_kernel.h>

#include <fstream>
#include <iostream>
#include <list>

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
      I::IntersectionData inter[2];
      for(int k = 0; k < 2; k++) {
        inter[k].quad = quadData[j];
        line.copyData(inter[k].line);
        inter[k].intPos = intParams[k];
        inter[k].intPos = intParams[(k + 1) % 2];
      }
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

template <int dim, typename fptype>
std::list<Geometry::Quadric<dim, fptype>> parseQuadrics(
    const char *fname) {
  std::ifstream file(fname);
  using Qf = Geometry::Quadric<dim, fptype>;
  std::list<Qf> quads;
  int imQuads = 0;
  while(!file.eof()) {
    Qf q;
    file >> q;
    QuadricClassify::QuadType type =
        QuadricClassify::classifyQuadric(q);
    if(!QuadricClassify::isImaginary(type) &&
       type != QuadricClassify::QUADT_ERROR)
      quads.push_back(q);
    else
      imQuads++;
  }
  std::cout << imQuads << " imaginary quadrics, "
            << quads.size() << " real quadrics.\n";
  return quads;
}

int main(int argc, char **argv) {
  if(argc < 2) {
    std::cout << "Not enough parameters. Usage:\n"
              << argv[0] << " quadricfile\n";
  }
  const char *fname = argv[1];
  constexpr const int numBlocks = 128;
  constexpr const int numTPB = 32;
  constexpr const int numThreads = numBlocks * numTPB;
  curandState *curng;
  cudaError_t err =
      cudaMalloc(&curng, sizeof(*curng) * numThreads);
  if(err != cudaSuccess || curng == NULL) {
    std::cout << "Error allocated GPU memory 1\n";
    return 1;
  }
  int seed = 0;
  setupKernel<<<numBlocks, numTPB>>>(curng, seed);

  auto quads = parseQuadrics<dim, fptype>(fname);
  constexpr const int numLines = 1024;
  Q::QuadricData *cuQuadMem;
  err = cudaMalloc(&cuQuadMem,
                   sizeof(*cuQuadMem) * quads.size());
  if(err != cudaSuccess || cuQuadMem == NULL) {
    cudaFree(curng);
    std::cout << "Error allocated GPU memory 2\n";
    return 1;
  }
  int i = 0;
  for(auto q : quads) {
    q.cudaCopy(&cuQuadMem[i]);
    i++;
  }

  I::IntersectionData *cuIntMem;
  int maxNumInts = quads.size() * numLines;
  err =
      cudaMalloc(&cuIntMem, sizeof(*cuIntMem) * maxNumInts);
  if(err != cudaSuccess || cuQuadMem == NULL) {
    cudaFree(curng);
    cudaFree(cuQuadMem);
    std::cout << "Error allocated GPU memory 3\n";
    return 1;
  }

  computeIntersections<uDistLine, 1,
                       1><<<numBlocks, numTPB>>>(
      curng, cuQuadMem, cuIntMem);

  cudaFree(curng);
  cudaFree(cuQuadMem);
  cudaFree(cuIntMem);
  return 0;
}
