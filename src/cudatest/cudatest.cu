
#include "cudadef.h"

#include "line.hpp"
#include "quadrics.hpp"
#include "accurate_intersections.hpp"

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
	using O = Geometry::Origin<dim, fptype>;
	using V = Geometry::Vector<dim, fptype>;
	using P = Geometry::Point<dim, fptype>;
	using L = Geometry::Line<dim, fptype>;
  using Q = Geometry::Quadric<dim, fptype>;
	using I = Geometry::Intersection<dim, fptype>;
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
  std::cout << "Copied Quadrics:\n" << src << "\n" << dest << "\n";
	V lineOff;
	P lineInt(lineOff);
	V lineDir;
	for(int i = 0; i < dim; i++)
		lineDir.set(i, MathFuncs::MathFuncs<fptype>::sin(static_cast<fptype>(i + 0.2345)));
	lineDir.normalize();
	std::shared_ptr<V::VectorData> vcudamem = lineDir.cudaCopy();
	V lineDirDest;
	lineDirDest.cudaRetrieve(vcudamem);
	std::cout << "\nCopied Vectors:\n" << lineDir << "\n\n" << lineDirDest << "\n";
	L line(lineInt, lineDir);
	I iSrc(src, line, 0.23426, 539.234);
	auto icudamem = iSrc.cudaCopy();
	I iDest;
	iDest.cudaRetrieve(icudamem);
	std::cout << "\nCopied Intersections:\n" << iSrc << "\n\n" << iDest << "\n";
  return 0;
}
