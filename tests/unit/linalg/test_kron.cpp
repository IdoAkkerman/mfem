// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "../../../linalg/kronecker.cpp"
#include "mfem.hpp"
#include "unit_tests.hpp"

#include <algorithm>
#include <random>
#include <type_traits>

using namespace mfem;

namespace {

constexpr real_t kTol =
    std::is_same_v<real_t, float> ? real_t(1e-4) : real_t(1e-12);

DenseMatrix MakeRandomMatrix(int rows, int cols, std::mt19937 &rng) {
  std::uniform_real_distribution<real_t> dist(-1.0, 1.0);
  DenseMatrix M(rows, cols);
  for (int j = 0; j < cols; ++j) {
    for (int i = 0; i < rows; ++i) {
      M(i, j) = dist(rng);
    }
  }
  return M;
}

Vector MakeRandomVector(int size, std::mt19937 &rng) {
  std::uniform_real_distribution<real_t> dist(-1.0, 1.0);
  Vector v(size);
  for (int i = 0; i < size; ++i) {
    v(i) = dist(rng);
  }
  return v;
}

DenseMatrix MakeIdentity(int n) {
  DenseMatrix I(n, n);
  I = 0.0;
  for (int i = 0; i < n; ++i) {
    I(i, i) = 1.0;
  }
  return I;
}

// Explicit A ⊗ B (B is the fastest-varying block)
DenseMatrix ExplicitKroneckerPair(const DenseMatrix &A, const DenseMatrix &B) {
  const int mA = A.Height(), nA = A.Width();
  const int mB = B.Height(), nB = B.Width();
  DenseMatrix K(mA * mB, nA * nB);

  for (int iA = 0; iA < mA; ++iA) {
    for (int jA = 0; jA < nA; ++jA) {
      const real_t a_val = A(iA, jA);
      for (int iB = 0; iB < mB; ++iB) {
        for (int jB = 0; jB < nB; ++jB) {
          K(iA * mB + iB, jA * nB + jB) = a_val * B(iB, jB);
        }
      }
    }
  }
  return K;
}

// Matches KroneckerProductOperator(D_1, ..., D_N) -> D_N ⊗ ... ⊗ D_1
DenseMatrix ExplicitKronecker(const DenseMatrix &D1) { return D1; }

template <typename... Rest>
DenseMatrix ExplicitKronecker(const DenseMatrix &D1, const Rest &...rest) {
  return ExplicitKroneckerPair(ExplicitKronecker(rest...), D1);
}

// Vector Kronecker product a ⊗ b (b is the fastest-varying index)
Vector KronVec(const Vector &a, const Vector &b) {
  Vector out(a.Size() * b.Size());
  for (int i = 0; i < a.Size(); ++i) {
    for (int j = 0; j < b.Size(); ++j) {
      out(i * b.Size() + j) = a(i) * b(j);
    }
  }
  return out;
}

real_t RelativeError(const Vector &actual, const Vector &expected) {
  Vector diff(expected);
  diff -= actual;
  const real_t denom = std::max(real_t(1.0), expected.Norml2());
  return diff.Norml2() / denom;
}

template <typename... Matrices>
void CheckAgainstExplicit(const Matrices &...mats) {
  KroneckerProductOperator op(mats...); // Also exercises CTAD
  DenseMatrix K = ExplicitKronecker(mats...);

  REQUIRE(op.Height() == K.Height());
  REQUIRE(op.Width() == K.Width());

  std::mt19937 rng(12345);
  Vector x = MakeRandomVector(K.Width(), rng);
  Vector y_op(K.Height()), y_ref(K.Height());
  y_op = 999.0; // Ensure Mult overwrites rather than accumulates

  op.Mult(x, y_op);
  K.Mult(x, y_ref);

  REQUIRE(RelativeError(y_op, y_ref) < kTol);
}

} // namespace

TEST_CASE("KroneckerProductOperator hand-calculated ordering",
          "[Linalg][Kronecker]") {
  // D_1 = [[1, 2], [3, 4]], D_2 = [[0, 5], [6, 7]]
  // (D_2 ⊗ D_1) * [1, 2, 3, 4]^T = [55, 125, 107, 241]^T
  // (Checks both arithmetic and non-commuting argument order independently.)
  DenseMatrix D1(2, 2), D2(2, 2);
  D1(0, 0) = 1.0;
  D1(0, 1) = 2.0;
  D1(1, 0) = 3.0;
  D1(1, 1) = 4.0;

  D2(0, 0) = 0.0;
  D2(0, 1) = 5.0;
  D2(1, 0) = 6.0;
  D2(1, 1) = 7.0;

  KroneckerProductOperator op(D1, D2);
  REQUIRE(op.Height() == 4);
  REQUIRE(op.Width() == 4);

  Vector x({1.0, 2.0, 3.0, 4.0});
  Vector y(4);
  op.Mult(x, y);

  Vector expected({55.0, 125.0, 107.0, 241.0});
  REQUIRE(RelativeError(y, expected) < kTol);
}

TEST_CASE("KroneckerProductOperator matches explicit DenseMatrix",
          "[Linalg][Kronecker]") {
  std::mt19937 rng(2026);

  SECTION("N = 1 (square and rectangular)") {
    CheckAgainstExplicit(MakeRandomMatrix(4, 4, rng));
    CheckAgainstExplicit(MakeRandomMatrix(5, 2, rng));
    CheckAgainstExplicit(MakeRandomMatrix(2, 5, rng));
  }

  SECTION("N = 2 (square, tall, wide, and mixed)") {
    CheckAgainstExplicit(MakeRandomMatrix(3, 3, rng),
                         MakeRandomMatrix(4, 4, rng));
    CheckAgainstExplicit(MakeRandomMatrix(4, 2, rng),
                         MakeRandomMatrix(5, 3, rng));
    CheckAgainstExplicit(MakeRandomMatrix(2, 4, rng),
                         MakeRandomMatrix(3, 5, rng));
    CheckAgainstExplicit(MakeRandomMatrix(5, 2, rng),
                         MakeRandomMatrix(2, 6, rng));
  }

  SECTION("N = 3 (exercises interior mode with L > 1 and R > 1)") {
    CheckAgainstExplicit(MakeRandomMatrix(3, 3, rng),
                         MakeRandomMatrix(4, 4, rng),
                         MakeRandomMatrix(5, 5, rng));
    CheckAgainstExplicit(MakeRandomMatrix(4, 2, rng),
                         MakeRandomMatrix(3, 5, rng),
                         MakeRandomMatrix(2, 3, rng));
    CheckAgainstExplicit(MakeRandomMatrix(2, 5, rng),
                         MakeRandomMatrix(6, 2, rng),
                         MakeRandomMatrix(3, 4, rng));
  }

  SECTION("N = 4 (higher-dimensional mixed rectangular)") {
    CheckAgainstExplicit(
        MakeRandomMatrix(2, 3, rng), MakeRandomMatrix(3, 2, rng),
        MakeRandomMatrix(4, 2, rng), MakeRandomMatrix(2, 3, rng));
  }

  SECTION("1x1 and degenerate dimensions") {
    CheckAgainstExplicit(MakeRandomMatrix(1, 1, rng));
    CheckAgainstExplicit(MakeRandomMatrix(1, 4, rng),
                         MakeRandomMatrix(3, 1, rng),
                         MakeRandomMatrix(2, 2, rng));
    CheckAgainstExplicit(MakeRandomMatrix(3, 2, rng),
                         MakeRandomMatrix(1, 1, rng),
                         MakeRandomMatrix(4, 3, rng));
  }
}

TEST_CASE("KroneckerProductOperator edge cases and statefulness",
          "[Linalg][Kronecker]") {
  std::mt19937 rng(777);

  SECTION("Identity factors preserve input vector") {
    KroneckerProductOperator op(MakeIdentity(3), MakeIdentity(4),
                                MakeIdentity(2));
    Vector x = MakeRandomVector(24, rng);
    Vector y(24);
    op.Mult(x, y);
    REQUIRE(RelativeError(y, x) < kTol);
  }

  SECTION("Zero matrix factor produces zero output") {
    DenseMatrix Z(3, 4);
    Z = 0.0;
    KroneckerProductOperator op(MakeRandomMatrix(2, 3, rng), Z,
                                MakeRandomMatrix(4, 2, rng));
    Vector x = MakeRandomVector(op.Width(), rng);
    Vector y(op.Height());
    y = 42.0;
    op.Mult(x, y);
    REQUIRE(y.Norml2() == real_t(0.0));
  }

  SECTION(
      "Repeated Mult calls do not corrupt internal buffers (odd & even N)") {
    auto D1 = MakeRandomMatrix(3, 4, rng);
    auto D2 = MakeRandomMatrix(5, 2, rng);
    auto D3 = MakeRandomMatrix(2, 3, rng);

    KroneckerProductOperator op_even(D1, D2);
    DenseMatrix K_even = ExplicitKronecker(D1, D2);

    KroneckerProductOperator op_odd(D1, D2, D3);
    DenseMatrix K_odd = ExplicitKronecker(D1, D2, D3);

    for (int call = 0; call < 4; ++call) {
      Vector x_e = MakeRandomVector(op_even.Width(), rng);
      Vector y_e(op_even.Height()), ref_e(op_even.Height());
      op_even.Mult(x_e, y_e);
      K_even.Mult(x_e, ref_e);
      REQUIRE(RelativeError(y_e, ref_e) < kTol);

      Vector x_o = MakeRandomVector(op_odd.Width(), rng);
      Vector y_o(op_odd.Height()), ref_o(op_odd.Height());
      op_odd.Mult(x_o, y_o);
      K_odd.Mult(x_o, ref_o);
      REQUIRE(RelativeError(y_o, ref_o) < kTol);
    }
  }

  SECTION("In-place aliasing Mult(x, x) for square operator") {
    auto D1 = MakeRandomMatrix(3, 3, rng);
    auto D2 = MakeRandomMatrix(4, 4, rng);
    auto D3 = MakeRandomMatrix(2, 2, rng);

    KroneckerProductOperator op(D1, D2, D3);
    DenseMatrix K = ExplicitKronecker(D1, D2, D3);

    Vector x = MakeRandomVector(op.Width(), rng);
    Vector ref(op.Height());
    K.Mult(x, ref);

    op.Mult(x, x);
    REQUIRE(RelativeError(x, ref) < kTol);
  }
}

TEST_CASE("KroneckerProductOperator algebraic identities (larger sizes)",
          "[Linalg][Kronecker]") {
  std::mt19937 rng(9001);

  SECTION("Separable vector action: (D3 ⊗ D2 ⊗ D1)(v3 ⊗ v2 ⊗ v1)") {
    auto D1 = MakeRandomMatrix(10, 8, rng);
    auto D2 = MakeRandomMatrix(7, 9, rng);
    auto D3 = MakeRandomMatrix(6, 11, rng);

    Vector v1 = MakeRandomVector(8, rng);
    Vector v2 = MakeRandomVector(9, rng);
    Vector v3 = MakeRandomVector(11, rng);

    Vector x = KronVec(v3, KronVec(v2, v1));
    Vector y(10 * 7 * 6);

    KroneckerProductOperator op(D1, D2, D3);
    op.Mult(x, y);

    Vector w1(10), w2(7), w3(6);
    D1.Mult(v1, w1);
    D2.Mult(v2, w2);
    D3.Mult(v3, w3);
    Vector expected = KronVec(w3, KronVec(w2, w1));

    REQUIRE(RelativeError(y, expected) < kTol);
  }

  SECTION("Mixed-product property: (A3 B3 ⊗ A2 B2 ⊗ A1 B1)x") {
    auto A1 = MakeRandomMatrix(6, 5, rng), B1 = MakeRandomMatrix(5, 7, rng);
    auto A2 = MakeRandomMatrix(4, 6, rng), B2 = MakeRandomMatrix(6, 3, rng);
    auto A3 = MakeRandomMatrix(5, 4, rng), B3 = MakeRandomMatrix(4, 6, rng);

    DenseMatrix C1(6, 7), C2(4, 3), C3(5, 6);
    mfem::Mult(A1, B1, C1);
    mfem::Mult(A2, B2, C2);
    mfem::Mult(A3, B3, C3);

    KroneckerProductOperator opA(A1, A2, A3);
    KroneckerProductOperator opB(B1, B2, B3);
    KroneckerProductOperator opC(C1, C2, C3);

    Vector x = MakeRandomVector(opB.Width(), rng);
    Vector tmp(opB.Height()), y_composed(opA.Height()), y_direct(opC.Height());

    opB.Mult(x, tmp);
    opA.Mult(tmp, y_composed);
    opC.Mult(x, y_direct);

    REQUIRE(RelativeError(y_composed, y_direct) < kTol);
  }

  SECTION("Linearity: K(alpha * x + beta * z) == alpha * Kx + beta * Kz") {
    KroneckerProductOperator op(MakeRandomMatrix(5, 4, rng),
                                MakeRandomMatrix(3, 6, rng),
                                MakeRandomMatrix(4, 3, rng));

    const real_t alpha = 2.75, beta = -1.5;
    Vector x = MakeRandomVector(op.Width(), rng);
    Vector z = MakeRandomVector(op.Width(), rng);

    Vector xz(x.Size());
    xz = 0.0;
    xz.Add(alpha, x);
    xz.Add(beta, z);

    Vector y_lhs(op.Height()), yx(op.Height()), yz(op.Height());
    op.Mult(xz, y_lhs);
    op.Mult(x, yx);
    op.Mult(z, yz);

    Vector y_rhs(op.Height());
    y_rhs = 0.0;
    y_rhs.Add(alpha, yx);
    y_rhs.Add(beta, yz);

    REQUIRE(RelativeError(y_lhs, y_rhs) < kTol);
  }
}
