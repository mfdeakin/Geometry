
#ifndef _ARRAY_HPP_
#define _ARRAY_HPP_

#include "cudadef.h"

template <typename T, int sz>
struct Array {
  T data[sz];

  CUDA_CALLABLE Array() {}

  CUDA_CALLABLE Array(const Array<T, sz> &src) {
    for(int i = 0; i < sz; i++) data[i] = src[i];
  }

  CUDA_CALLABLE Array(const T src[sz]) {
    for(int i = 0; i < sz; i++) data[i] = src[i];
  }

  CUDA_CALLABLE Array(std::initializer_list<T> src) {
    std::copy(src.begin(), src.end(), data);
  }

  CUDA_CALLABLE ~Array() {}

  CUDA_CALLABLE const T &operator[](int idx) const {
    assert(idx >= 0);
    assert(idx < sz);
    return data[idx];
  }

  CUDA_CALLABLE T &operator[](int idx) {
    assert(idx >= 0);
    assert(idx < sz);
    return data[idx];
  }

  CUDA_CALLABLE Array<T, sz> operator=(
      const Array<T, sz> &src) {
    for(int i = 0; i < sz; i++) data[i] = src[i];
    return *this;
  }

  CUDA_CALLABLE static constexpr int size() { return sz; }
};

#endif
