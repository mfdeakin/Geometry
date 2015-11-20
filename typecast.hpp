
#ifndef _TYPECAST_HPP_
#define _TYPECAST_HPP_

#include "mpreal.hpp"

namespace Typecast {

template <typename t1, typename t2>
struct typecast {};

template <>
struct typecast<float, float> {
  typedef float higherPrec;
  typedef float lowerPrec;
};

template <>
struct typecast<float, double> {
  typedef double higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<double, float> {
  typedef double higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<float, long double> {
  typedef long doublehigherPrec;
  typedef float lowerPrec;
};

template <>
struct typecast<long double, float> {
  typedef long double higherPrec;
  typedef float lowerPrec;
};

template <>
struct typecast<float, mpfr::mpreal> {
  typedef mpfr::mpreal higherPrec;
  typedef float lowerPrec;
};

template <>
struct typecast<mpfr::mpreal, float> {
  typedef mpfr::mpreal higherPrec;
  typedef float lowerPrec;
};

template <>
struct typecast<double, double> {
  typedef double higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<double, long double> {
  typedef long double higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<long double, double> {
  typedef long double higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<double, mpfr::mpreal> {
  typedef mpfr::mpreal higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<mpfr::mpreal, double> {
  typedef mpfr::mpreal higherPrec;
  typedef double lowerPrec;
};

template <>
struct typecast<long double, long double> {
  typedef long double higherPrec;
  typedef long double lowerPrec;
};

template <>
struct typecast<long double, mpfr::mpreal> {
  typedef mpfr::mpreal higherPrec;
  typedef long double lowerPrec;
};

template <>
struct typecast<mpfr::mpreal, long double> {
  typedef mpfr::mpreal higherPrec;
  typedef long double lowerPrec;
};

template <>
struct typecast<mpfr::mpreal, mpfr::mpreal> {
  typedef mpfr::mpreal higherPrec;
  typedef mpfr::mpreal lowerPrec;
};

template <>
struct typecast<short int, short int> {
  typedef short int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<short int, int> {
  typedef int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<short int, long int> {
  typedef long int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<short int, long long int> {
  typedef long long int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<int, short int> {
  typedef int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<long int, short int> {
  typedef long int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<long long int, short int> {
  typedef long long int higherPrec;
  typedef short int lowerPrec;
};

template <>
struct typecast<int, int> {
  typedef int higherPrec;
  typedef int lowerPrec;
};

template <>
struct typecast<int, long int> {
  typedef long int higherPrec;
  typedef int lowerPrec;
};

template <>
struct typecast<int, long long int> {
  typedef long long int higherPrec;
  typedef int lowerPrec;
};

template <>
struct typecast<long int, int> {
  typedef long int higherPrec;
  typedef int lowerPrec;
};

template <>
struct typecast<long long int, int> {
  typedef long long int higherPrec;
  typedef int lowerPrec;
};

template <>
struct typecast<long int, long int> {
  typedef long int higherPrec;
  typedef long int lowerPrec;
};

template <>
struct typecast<long int, long long int> {
  typedef long long int higherPrec;
  typedef long int lowerPrec;
};

template <>
struct typecast<long long int, long int> {
  typedef long long int higherPrec;
  typedef long int lowerPrec;
};

template <>
struct typecast<short unsigned int, short unsigned int> {
  typedef short unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<short unsigned int, unsigned int> {
  typedef unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<short unsigned int, long unsigned int> {
  typedef long unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<short unsigned int,
                long long unsigned int> {
  typedef long long unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<unsigned int, short unsigned int> {
  typedef unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<long unsigned int, short unsigned int> {
  typedef long unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<long long unsigned int,
                short unsigned int> {
  typedef long long unsigned int higherPrec;
  typedef short unsigned int lowerPrec;
};

template <>
struct typecast<unsigned int, unsigned int> {
  typedef unsigned int higherPrec;
  typedef unsigned int lowerPrec;
};

template <>
struct typecast<unsigned int, long unsigned int> {
  typedef long unsigned int higherPrec;
  typedef unsigned int lowerPrec;
};

template <>
struct typecast<unsigned int, long long unsigned int> {
  typedef long long unsigned int higherPrec;
  typedef unsigned int lowerPrec;
};

template <>
struct typecast<long unsigned int, unsigned int> {
  typedef long unsigned int higherPrec;
  typedef unsigned int lowerPrec;
};

template <>
struct typecast<long long unsigned int, unsigned int> {
  typedef long long unsigned int higherPrec;
  typedef unsigned int lowerPrec;
};

template <>
struct typecast<long unsigned int, long unsigned int> {
  typedef long unsigned int higherPrec;
  typedef long unsigned int lowerPrec;
};

template <>
struct typecast<long unsigned int, long long unsigned int> {
  typedef long long unsigned int higherPrec;
  typedef long unsigned int lowerPrec;
};

template <>
struct typecast<long long unsigned int, long unsigned int> {
  typedef long long unsigned int higherPrec;
  typedef long unsigned int lowerPrec;
};
}

#endif
