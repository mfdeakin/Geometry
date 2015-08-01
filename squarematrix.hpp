
#ifndef _SQUAREMATRIX_HPP_
#define _SQUAREMATRIX_HPP_

#include "matrix.hpp"
#include "vector.hpp"

template <unsigned _dim, typename fptype>
class SquareMatrix : public Matrix<_dim, _dim, fptype>;

template <typename fptype>
class SquareMatrix<1, fptype> : public Matrix<1, 1, fptype> {
public:
  SquareMatrix() {}
}

#endif
