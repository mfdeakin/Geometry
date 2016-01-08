
#include "cudadef.h"

#include "line.hpp"
#include "quadrics.hpp"

template <typename Q>
__global__ void CMTest(typename Q::QuadricData *qIn,
                       typename Q::QuadricData *qOut,
                       int numQuads) {
  int base = blockDim.x * blockIdx.x;
  int idx = base + threadIdx.x;
  if(idx < Q::numCoeffs)
    qOut->coeffs[idx] = qIn->coeffs[idx];
}

int main(int argc, char **argv) {
  constexpr const int dim = 3;
  using fptype = float;
  using Q = Geometry::Quadric<dim, fptype>;
  Geometry::Vector<dim, fptype> v;
  Geometry::Point<dim, fptype> p(v);
  Q src;
  for(int i = 0; i < src.numCoeffs; i++)
    src.setCoeff(i, (2 * i + 1) * (2 * i + 2) / 2);
  auto spCudaMem = src.cudaCopy();
  Q dest;
  auto newMem = dest.cudaCopy();

  constexpr const int numBlocks =
      Geometry::Quadric<dim, fptype>::numCoeffs;
  constexpr const int thrdPerBlk = 1;

  cudaDeviceProp cudev;
  constexpr const int devNum = 0;
  cudaSetDevice(devNum);
  cudaGetDeviceProperties(&cudev, 0);

  CMTest<Q><<<numBlocks, thrdPerBlk>>>(spCudaMem.get(),
                                       newMem.get(), 1);
  dest.cudaRetrieve(newMem);
  std::cout << src << "\n" << dest << "\n";
  return 0;
}
