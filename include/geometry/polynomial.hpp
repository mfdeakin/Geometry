
#ifndef _POLYNOMIAL_HPP_
#define _POLYNOMIAL_HPP_

#include <memory>
#include <cmath>
#include <ostream>
#include <istream>
#include <iostream>

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
    return coeffs[coeff];
  }

  CUDA_CALLABLE fptype &set(int coeff, fptype val) {
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

  friend std::ostream &operator<<(
      std::ostream &os,
      const Polynomial<order, fptype> &p) {
    bool once = false;
    for(int i = p.numCoeffs - 1; i >= 0; i--) {
      if(once) {
        os << " + ";
      }
      once = true;
      os << p.get(i) << " t^" << i;
    }
    return os;
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
 public:
  PolyIsSolvable() : PolynomialBase<0, fptype>() {}
  PolyIsSolvable(const PolyIsSolvable<0, fptype, true> &p)
      : PolynomialBase<0, fptype>(p) {}
  Array<fptype, 1> calcRoots() { return {}; }
};

template <typename fptype>
class PolyIsSolvable<1, fptype, true>
    : public PolynomialBase<1, fptype> {
 public:
  PolyIsSolvable() : PolynomialBase<1, fptype>() {}
  PolyIsSolvable(const PolyIsSolvable<1, fptype, true> &p)
      : PolynomialBase<1, fptype>(p) {}
  Array<fptype, 1> calcRoots() {
    return {-PolynomialBase<1, fptype>::get(0) /
            PolynomialBase<1, fptype>::get(1)};
  }
};

template <typename fptype>
class PolyIsSolvable<2, fptype, true>
    : public PolynomialBase<2, fptype> {
 public:
  PolyIsSolvable() : PolynomialBase<2, fptype>() {}
  PolyIsSolvable(const PolyIsSolvable<2, fptype, true> &p)
      : PolynomialBase<2, fptype>(p) {}
  Array<fptype, 2> calcRoots() {
    fptype sqCoeff = this->get(2);
    fptype linCoeff = this->get(1);
    fptype constant = this->get(0);
    if(sqCoeff == 0.0) {
      return Array<fptype, 2>(
          {{-constant / linCoeff, fptype(NAN)}});
    }
    fptype disc = fptype(linCoeff * linCoeff / 4.0) -
                  fptype(sqCoeff * constant);
    if(disc < 0) {
      return Array<fptype, 2>({{fptype(NAN), fptype(NAN)}});
    }
    fptype fracPart =
        -MathFuncs::MathFuncs<fptype>::copysign(
            MathFuncs::MathFuncs<fptype>::fabs(linCoeff /
                                               2.0) +
                MathFuncs::MathFuncs<fptype>::sqrt(disc),
            linCoeff);
    Array<fptype, 2> roots = {
        {fracPart / sqCoeff, constant / fracPart}};
    return roots;
  }
};

template <int order, typename fptype>
class Polynomial
    : public PolyIsSolvable<order, fptype, (order < 3)> {
 public:
  Polynomial()
      : PolyIsSolvable<order, fptype, (order < 3)>() {}
  Polynomial(const Polynomial<order, fptype> &p)
      : PolyIsSolvable<order, fptype, (order < 3)>(p) {}
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
