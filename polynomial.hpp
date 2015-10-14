
#ifndef _POLYNOMIAL_HPP_
#define _POLYNOMIAL_HPP_

#include <memory>
#include <array>
#include <cmath>

#include "typecast.hpp"

namespace Geometry {

template <int order, typename fptype>
class Polynomial;

template <int order, typename fptype>
class PolynomialBase {
 public:
  PolynomialBase()
      : coeffs(new fptype[numCoeffs],
               std::default_delete<fptype[]>()) {
    for(int i = 0; i < order; i++) coeffs.get()[i] = 0.0;
  }

  PolynomialBase(const std::array<fptype, order> &coeffs)
      : coeffs(new fptype[numCoeffs],
               std::default_delete<fptype[]>()) {
    for(int i = 0; i < order; i++)
      this->coeffs.get()[i] = coeffs[i];
  }

  PolynomialBase(const PolynomialBase<order, fptype> &p)
      : coeffs(p.coeffs) {}

  fptype get(int coeff) const {
    assert(coeff >= 0);
    assert(coeff < numCoeffs);
    return coeffs.get()[coeff];
  }

  fptype &set(int coeff, fptype val) {
    assert(coeff >= 0);
    assert(coeff < numCoeffs);
    if(coeffs.unique() == false) {
      fptype *newCoeffs = new fptype[numCoeffs];
      for(int i = 0; i < numCoeffs; i++)
        newCoeffs[i] = get(i);
      coeffs.reset(newCoeffs,
                   std::default_delete<fptype[]>());
    }
    coeffs.get()[coeff] = val;
    return coeffs.get()[coeff];
  }

  template <typename valType>
  fptype eval(valType v) {
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
  std::shared_ptr<fptype> coeffs;
};

template <int order, typename fptype, bool isSolvable>
class PolyIsSolvable;

template <int order, typename fptype>
class PolyIsSolvable<order, fptype, false>
    : public PolynomialBase<order, fptype> {};

template <typename fptype>
class PolyIsSolvable<0, fptype, true>
    : public PolynomialBase<0, fptype> {
  std::array<fptype, 1> calcRoots() { return {}; }
};

template <typename fptype>
class PolyIsSolvable<1, fptype, true>
    : public PolynomialBase<1, fptype> {
 public:
  std::array<fptype, 1> calcRoots() {
    return {-PolynomialBase<1, fptype>::get(0) /
            PolynomialBase<1, fptype>::get(1)};
  }
};

template <typename fptype>
class PolyIsSolvable<2, fptype, true>
    : public PolynomialBase<2, fptype> {
 public:
  std::array<fptype, 2> calcRoots() {
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
