
#include "cudadef.h"

#include "line.hpp"
#include "quadrics.hpp"

int main(int argc, char **argv) {
  constexpr const int dim = 3;
  using fptype = float;
  Geometry::Quadric<dim, fptype> q;
  auto spCudaMem = q.cudaCopy();
  return 0;
}
