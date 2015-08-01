
#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_

#include "matrix.hpp"
#include "cassert.h"

template <unsigned _dim, typename fptype>
class Vector : public Matrix<_dim, 1, fptype> {
public:
  Vector() {}
  Vector(const Vector &src) : Matrix<_dim, 1, fptype>(src) {}

  fptype norm() {
    fptype sumSquared;
    fptype extra;
    for(unsigned i = 0; i < _dim; i++) {
      fptype normalized = this->array[i] * this->array[i] - extra;
      fptype tmp = sumSquared + normalized;
      extra = (tmp - sumSquared) - this->array[i];
      sumSquared = tmp;
    }
    return sqrt(sumSquared);
  }
};

#endif
