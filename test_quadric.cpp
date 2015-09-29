
#include <array>

#include <gtest/gtest.h>

#include "quadrics.hpp"
#include "origin.hpp"
#include "quadric_classify.hpp"
#include "accurate_math.hpp"

TEST(Quadric, LineIntersection) {
  using fptype = float;
  constexpr const int dim = 3;
  constexpr const unsigned numCoeffs =
      (dim + 2) * (dim + 1) / 2;
  constexpr const fptype eps = 1e-5;
  struct teststruct {
    std::array<fptype, numCoeffs> coeffs;
    std::array<fptype, dim> lineDir;
    std::array<fptype, dim> lineInt;
    std::array<fptype, dim> roots[2];
  } tests[] = {
      {{1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}},
      {{1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{NAN, NAN, NAN}, {NAN, NAN, NAN}}},
      {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{0.0, 0.0, 0.0}, {NAN, NAN, NAN}}},
      {{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0},
       {1.0, 1.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{-11.0, -1.0, 0.0}, {0.0, 10.0, 0.0}}},
      {{1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, 60.0},
       {{-1.0, 0.0, 60.0}, {1.0, 0.0, 60.0}}},
      {{1.0, 1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, 0.0},
       {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}},
  };
  Geometry::Origin<dim, fptype> o;
  for(auto t : tests) {
    Geometry::Vector<dim, fptype> offset(t.lineInt);
    Geometry::Point<dim, fptype> intercept(o, offset);
    Geometry::Vector<dim, fptype> dir(t.lineDir);
    Geometry::Line<dim, fptype> l(intercept, dir);
    Geometry::Quadric<dim, fptype> q(o);
    for(unsigned i = 0; i < numCoeffs; i++)
      q.coeff(i) = t.coeffs[i];
    std::cout << "Quadric: " << q << "\n";
    std::cout << "Line: " << l << "\n";
    auto intersects = q.calcLineIntersect(l);
    for(unsigned i = 0; i < 2; i++) {
      Geometry::Point<dim, fptype> expected(
          o, Geometry::Vector<dim, fptype>(t.roots[i]));
      Geometry::Point<dim, fptype> p(intersects[i]);
      std::cout << "Expected: " << expected << " vs. " << p
                << "\n";
      /* Either none or all coordinates are NAN */
      if(std::isnan(t.roots[i][0])) {
        for(unsigned j = 0; j < dim; j++)
          EXPECT_EQ(std::isnan(p.getOffset()(j)), true);
      } else {
        Geometry::Vector<dim, fptype> delta =
            p.ptDiff(expected);
        fptype diffMag = delta.norm();
        EXPECT_NEAR(diffMag, 0.0, eps);
      }
    }
  }
}
