

#include <gtest/gtest.h>

#include "quadrics.hpp"
#include "origin.hpp"
#include "quadric_classify.hpp"
#include "accurate_math.hpp"
#include "array.hpp"

TEST(Quadric, LineIntersection) {
  using fptype = float;
  constexpr const int dim = 3;
  constexpr const unsigned numCoeffs =
      Geometry::Quadric<dim, fptype>::numCoeffs;
  constexpr const int numRoots = 2;
  constexpr const fptype eps = 1e-3;
  struct teststruct {
    Array<fptype, numCoeffs> coeffs;
    Array<fptype, dim> lineDir;
    Array<fptype, dim> lineInt;
    Array<fptype, dim> roots[numRoots];
  } tests[] = {
      /* Hyperboloid of One Sheet */
      {{1.0, 1.0, -3.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, 1.0},
       {{static_cast<fptype>(-std::sqrt(19)), 0.0, 1.0},
        {static_cast<fptype>(std::sqrt(19)), 0.0, 1.0}}},
      /* Hyperboloid of Two Sheets */
      {{1.0, 1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, 1.0},
       {{0.0, 0.0, 1.0}, {NAN, NAN, NAN}}},
      /* Elliptic Paraboloid */
      {{1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, -1.0},
       {{-1.0, 0.0, -1.0}, {1.0, 0.0, -1.0}}},
      /* Hyperbolic Paraboloid */
      {{1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, -1.0},
       {{-1.0, 0.0, -1.0}, {1.0, 0.0, -1.0}}},
      /* Unit Sphere */
      {{1.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}},
      /* Imaginary Unit Sphere */
      {{1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{NAN, NAN, NAN}, {NAN, NAN, NAN}}},
      /* Imaginary Cone */
      {{1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{0.0, 0.0, 0.0}, {NAN, NAN, NAN}}},
      /* Real Cone */
      {{1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 4.0},
       {-10.0, 0.0, 0.0},
       {{-40.0 / 3.0, 0.0, -40.0 / 3.0}, {-8.0, 0.0, 8.0}}},
      /* Parabolic Cylinder */
      {{0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 4.0},
       {{-15.0, 0.0, 4.0}, {NAN, NAN, NAN}}},
      /* Hyperbolic Cylinder */
      {{1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {0.0, 1.0, 4.0},
       {-10.0, 0.0, 0.0},
       {{-10.0, static_cast<fptype>(std::sqrt(101.0)),
         static_cast<fptype>(4.0 * std::sqrt(101.0))},
        {-10.0, static_cast<fptype>(-std::sqrt(101.0)),
         static_cast<fptype>(-4.0 * std::sqrt(101.0))}}},
      /* Imaginary Elliptic Cylinder */
      {{1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{NAN, NAN, NAN}, {NAN, NAN, NAN}}},
      /* Real Elliptic Cylinder */
      {{1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}},
      /* Two Imaginary Parallel Planes, x=(+/-)i */
      {{1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{NAN, NAN, NAN}, {NAN, NAN, NAN}}},
      /* Two Real Parallel Planes, x=(+/-)1 */
      {{1.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {10.0, 0.0, 60.0},
       {{-1.0, 0.0, 60.0}, {1.0, 0.0, 60.0}}},
      /* Two Real Intersecting Planes/The line x=0, y=-1 */
      {{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0},
       {1.0, 1.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{-11.0, -1.0, 0.0}, {0.0, 10.0, 0.0}}},
      /* Two Imaginary Intersecting Planes */
      {{1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 1.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{NAN, NAN, NAN}, {NAN, NAN, NAN}}},
      /* Coincident Planes at x=0 */
      {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
       {1.0, 0.0, 0.0},
       {-10.0, 0.0, 0.0},
       {{0.0, 0.0, 0.0}, {NAN, NAN, NAN}}},
      /* Two Sheeted Hyperboloid;
       * a more stringent test with error ~1e-3 */
      {{static_cast<fptype>(std::sqrt(2)),
        static_cast<fptype>(std::exp(20)),
        static_cast<fptype>(-std::atan(1) * 4 / 1e3),
        static_cast<fptype>(std::log(7)),
        static_cast<fptype>(std::sin(2)),
        static_cast<fptype>(std::cos(2)), 2.0, -1.0,
        static_cast<fptype>(-std::sqrt(3)),
        static_cast<fptype>(-std::exp(10))},
       {0.0, 1.0, 4.0},
       {0.0, 0.0, 50.0},
       {{0.0, 0.0477354, 50.1909},
        {0.0, -0.0475539, 49.8098}}},
  };
  Geometry::Origin<dim, fptype> o;
  for(auto t : tests) {
    Geometry::Vector<dim, fptype> offset(t.lineInt);
    Geometry::Point<dim, fptype> intercept(o, offset);
    Geometry::Vector<dim, fptype> dir(t.lineDir);
    Geometry::Line<dim, fptype> l(intercept, dir);
    Geometry::Quadric<dim, fptype> q(o);
    for(unsigned i = 0; i < numCoeffs; i++)
      q.setCoeff(i, t.coeffs[i]);
    auto intersects = q.calcLineIntersect(l);
    for(unsigned i = 0; i < 2; i++) {
      Geometry::Point<dim, fptype> expected(
          o, Geometry::Vector<dim, fptype>(t.roots[i]));
      Geometry::Point<dim, fptype> p(intersects[i]);
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
