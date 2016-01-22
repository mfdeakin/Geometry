
#ifndef _POLYNOMIAL_HPP_
#define _POLYNOMIAL_HPP_

#include <memory>
#include <cmath>

#include "array.hpp"
#include "cudadef.h"
#include "typecast.hpp"

namespace Geometry {

template <int order, typename fptype>
class Polynomial;

template <int order, typename fptype>
class PolynomialBase {
 public:
  CUDA_CALLABLE PolynomialBase() {
    for(int i = 0; i < order; i++) set(i, NAN);
  }

  CUDA_CALLABLE PolynomialBase(
      const Array<fptype, order> &coeffs) {
    for(int i = 0; i < order; i++) set(i, coeffs[i]);
  }

  template <typename srctype>
  CUDA_CALLABLE PolynomialBase(
      const PolynomialBase<order, srctype> &p) {
    for(int i = 0; i < order; i++) set(i, p.get(i));
  }

  CUDA_CALLABLE fptype get(int coeff) const {
    assert(coeff >= 0);
    assert(coeff < numCoeffs);
    return coeffs[coeff];
  }

  CUDA_CALLABLE fptype &set(int coeff, fptype val) {
    assert(coeff >= 0);
    assert(coeff < numCoeffs);
    coeffs[coeff] = val;
    return coeffs[coeff];
  }

  template <typename valType>
  CUDA_CALLABLE fptype eval(valType v) {
    fptype sum = 0.0;
    for(int i = numCoeffs - 1; i >= 0; i--) {
      sum = std::fma(v, sum, get(i));
    }
    return sum;
  }

  static constexpr const int _order = order;
  static constexpr const int numCoeffs = order + 1;
  using _fptype = fptype;

 protected:
  Array<fptype, numCoeffs> coeffs;
};

template <int order, typename fptype, bool isSolvable>
class PolyIsSolvable;

template <int order, typename fptype>
class PolyIsSolvable<order, fptype, false>
    : public PolynomialBase<order, fptype> {};

template <typename fptype>
class PolyIsSolvable<0, fptype, true>
    : public PolynomialBase<0, fptype> {
  Array<fptype, 1> calcRoots() { return {}; }
};

template <typename fptype>
class PolyIsSolvable<1, fptype, true>
    : public PolynomialBase<1, fptype> {
 public:
  Array<fptype, 1> calcRoots() {
    return {-PolynomialBase<1, fptype>::get(0) /
            PolynomialBase<1, fptype>::get(1)};
  }
};

template <typename fptype>
class PolyIsSolvable<2, fptype, true>
    : public PolynomialBase<2, fptype> {
 public:
  Array<fptype, 2> calcRoots() {
    return AccurateMath::kahanQuadratic(
        PolynomialBase<2, fptype>::get(2),
        PolynomialBase<2, fptype>::get(1),
        PolynomialBase<2, fptype>::get(0));
  }
};

template <int order, typename fptype>
    class Polynomial : public PolyIsSolvable < order,
    fptype, order<3> {
 public:
  template <typename p2Type>
  Polynomial<
      order + p2Type::_order,
      typename Typecast::typecast<
          fptype, typename p2Type::_fptype>::higherPrec>
  product(p2Type rhs) {
    Polynomial<
        order + p2Type::_order,
        typename Typecast::typecast<
            fptype, typename p2Type::_fptype>::higherPrec>
        res;
    for(int i = 0; i < this->numCoeffs; i++) {
      for(int j = 0; j < p2Type::numCoeffs; j++) {
        fptype sum = std::fma(this->get(i), rhs.get(j),
                              res.get(i + j));
        res.set(i + j, sum);
      }
    }
    return res;
  }
};
}

#endif
