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

#include "../../../linalg/kronecker.hpp"
#include "mfem.hpp"
#include "unit_tests.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace mfem {

namespace kronecker_test {

template <size_t N>
const char *
BackendName(typename KroneckerProductOperator<N>::ContractionBackend b) {
  using Backend = typename KroneckerProductOperator<N>::ContractionBackend;
  switch (b) {
  case Backend::Auto:
    return "Auto";
  case Backend::BLAS:
    return "BLAS";
  case Backend::BlockedCPP:
    return "BlockedCPP";
  }
  return "Unknown";
}

template <size_t N>
const char *SymmetryName(typename KroneckerProductOperator<N>::SymmetryType s) {
  using Sym = typename KroneckerProductOperator<N>::SymmetryType;
  switch (s) {
  case Sym::None:
    return "None";
  case Sym::RowEvenOdd:
    return "RowEvenOdd";
  case Sym::ColEvenOdd:
    return "ColEvenOdd";
  case Sym::Centrosymmetric:
    return "Centrosymmetric";
  }
  return "Unknown";
}

real_t ContractionTol(int width, int height) {
  const real_t eps = std::numeric_limits<real_t>::epsilon();
  const real_t scale = std::sqrt(static_cast<real_t>(std::max(width, height)));
  return real_t(64.0) * eps * std::max(real_t(1.0), scale);
}

Vector MakeRandomVector(int size, int seed) {
  Vector v(size);
  v.Randomize(seed);
  for (int i = 0; i < size; ++i) {
    v(i) = real_t(2.0) * v(i) - real_t(1.0);
  }
  return v;
}

DenseMatrix MakeRandomMatrix(int rows, int cols, int seed) {
  const Vector vals = MakeRandomVector(rows * cols, seed);
  DenseMatrix M(rows, cols);
  std::copy_n(vals.HostRead(), rows * cols, M.HostWrite());
  return M;
}

DenseMatrix MakeIdentity(int n) {
  DenseMatrix I(n, n);
  I = real_t(0.0);
  for (int i = 0; i < n; ++i) {
    I(i, i) = real_t(1.0);
  }
  return I;
}

/// Exchange (anti-diagonal reversal) matrix J_n: (J_n)_{i, n - 1 - i} = 1.
DenseMatrix MakeExchange(int n) {
  DenseMatrix J(n, n);
  J = real_t(0.0);
  for (int i = 0; i < n; ++i) {
    J(i, n - 1 - i) = real_t(1.0);
  }
  return J;
}

/// Algebraic row even-odd projection: rows from R*(I+J)/2 or R*(I-J)/2.
/// mode: 0 = mixed even/odd rows, 1 = all even rows, 2 = all odd rows.
DenseMatrix MakeRowEvenOddMatrix(int m, int n, int seed, int mode) {
  const DenseMatrix R = MakeRandomMatrix(m, n, seed);
  const DenseMatrix I = MakeIdentity(n);
  const DenseMatrix J = MakeExchange(n);

  DenseMatrix P_even(I), P_odd(I);
  P_even += J;
  P_even *= real_t(0.5);
  P_odd -= J;
  P_odd *= real_t(0.5);

  DenseMatrix R_even(m, n), R_odd(m, n);
  mfem::Mult(R, P_even, R_even);
  mfem::Mult(R, P_odd, R_odd);

  if (mode == 1) {
    return R_even;
  }
  if (mode == 2) {
    return R_odd;
  }

  DenseMatrix A(m, n);
  for (int i = 0; i < m; ++i) {
    const DenseMatrix &src = (i % 3 != 0) ? R_even : R_odd;
    for (int j = 0; j < n; ++j) {
      A(i, j) = src(i, j);
    }
  }
  return A;
}

/// Algebraic column even-odd matrix: transpose of a row even-odd matrix.
DenseMatrix MakeColEvenOddMatrix(int m, int n, int seed, int mode) {
  const DenseMatrix B = MakeRowEvenOddMatrix(n, m, seed, mode);
  DenseMatrix A;
  A.Transpose(B);
  return A;
}

/// Algebraic centrosymmetric projection: 0.5 * (R + J_m * R * J_n).
DenseMatrix MakeCentrosymmetricMatrix(int m, int n, int seed) {
  const DenseMatrix R = MakeRandomMatrix(m, n, seed);
  DenseMatrix tmp(m, n), rev(m, n);
  mfem::Mult(MakeExchange(m), R, tmp);
  mfem::Mult(tmp, MakeExchange(n), rev);

  DenseMatrix A(R);
  A += rev;
  A *= real_t(0.5);
  return A;
}

/// Algebraic anti-centrosymmetric projection: 0.5 * (R - J_m * R * J_n).
DenseMatrix MakeAntiCentrosymmetricMatrix(int m, int n, int seed) {
  const DenseMatrix R = MakeRandomMatrix(m, n, seed);
  DenseMatrix tmp(m, n), rev(m, n);
  mfem::Mult(MakeExchange(m), R, tmp);
  mfem::Mult(tmp, MakeExchange(n), rev);

  DenseMatrix A(R);
  A -= rev;
  A *= real_t(0.5);
  return A;
}

/// Exact orthogonal DST-I matrix (rows and columns alternate even/odd).
DenseMatrix MakeDSTMatrix(int n) {
  DenseMatrix V(n, n);
  const double pi = std::acos(-1.0);
  const double scale = std::sqrt(2.0 / static_cast<double>(n + 1));
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const double arg = pi * static_cast<double>((i + 1) * (j + 1)) /
                         static_cast<double>(n + 1);
      V(i, j) = static_cast<real_t>(scale * std::sin(arg));
    }
  }
  return V;
}

/// Explicit A \otimes B (B is the fastest-varying block).
DenseMatrix ExplicitKroneckerPair(const DenseMatrix &A, const DenseMatrix &B) {
  const int mA = A.Height(), nA = A.Width();
  const int mB = B.Height(), nB = B.Width();
  DenseMatrix K(mA * mB, nA * nB);

  for (int jA = 0; jA < nA; ++jA) {
    for (int jB = 0; jB < nB; ++jB) {
      const int col = jA * nB + jB;
      for (int iA = 0; iA < mA; ++iA) {
        const real_t a_val = A(iA, jA);
        for (int iB = 0; iB < mB; ++iB) {
          K(iA * mB + iB, col) = a_val * B(iB, jB);
        }
      }
    }
  }
  return K;
}

/// Explicit (D_{N-1} \otimes ... \otimes D_0).
DenseMatrix ExplicitKronecker(const DenseMatrix &D0) { return D0; }

template <typename... Rest>
DenseMatrix ExplicitKronecker(const DenseMatrix &D0, const Rest &...rest) {
  return ExplicitKroneckerPair(ExplicitKronecker(rest...), D0);
}

/// Vector Kronecker product a \otimes b (b is the fastest-varying index).
Vector KronVec(const Vector &a, const Vector &b) {
  Vector out(a.Size() * b.Size());
  for (int i = 0; i < a.Size(); ++i) {
    const real_t ai = a(i);
    const int off = i * b.Size();
    for (int j = 0; j < b.Size(); ++j) {
      out(off + j) = ai * b(j);
    }
  }
  return out;
}

template <size_t N> Vector KronVecArray(const std::array<Vector, N> &vecs) {
  Vector res = vecs[0];
  for (size_t k = 1; k < N; ++k) {
    res = KronVec(vecs[k], res);
  }
  return res;
}

real_t RelativeError(const Vector &actual, const Vector &expected) {
  Vector diff(expected);
  diff -= actual;
  const real_t denom = std::max(real_t(1.0), expected.Norml2());
  return diff.Norml2() / denom;
}

real_t StrictRelativeError(const Vector &actual, const Vector &expected) {
  Vector diff(expected);
  diff -= actual;
  const real_t ref_norm = expected.Norml2();
  REQUIRE(ref_norm > real_t(0.0));
  return diff.Norml2() / ref_norm;
}

/// Check Mult, MultTranspose, and adjointness against an explicit DenseMatrix
/// for one caller-selected (backend, use_eo) configuration.
template <size_t N>
void VerifyExplicitConfig(
    KroneckerProductOperator<N> &op, const DenseMatrix &K_ref,
    typename KroneckerProductOperator<N>::ContractionBackend backend,
    bool use_eo, int seed = 101) {
  const char *backend_name = BackendName<N>(backend);
  const int h = op.Height();
  const int w = op.Width();
  CAPTURE(N, h, w, backend_name, use_eo);

  REQUIRE(h == K_ref.Height());
  REQUIRE(w == K_ref.Width());

  op.SetBackend(backend);
  op.SetUseEvenOdd(use_eo);
  REQUIRE(op.GetBackend() == backend);
  REQUIRE(op.GetUseEvenOdd() == use_eo);

  const Vector x = MakeRandomVector(w, seed);
  const Vector z = MakeRandomVector(h, seed + 1);

  Vector y_ref(h), w_ref(w);
  K_ref.Mult(x, y_ref);
  K_ref.MultTranspose(z, w_ref);

  Vector y(h), adj(w);
  y = real_t(999.0);
  adj = real_t(-999.0);

  op.Mult(x, y);
  op.MultTranspose(z, adj);

  const real_t tol = ContractionTol(w, h);
  {
    INFO("Mult vs explicit DenseMatrix");
    REQUIRE(RelativeError(y, y_ref) <= tol);
  }
  {
    INFO("MultTranspose vs explicit DenseMatrix");
    REQUIRE(RelativeError(adj, w_ref) <= tol);
  }
  {
    INFO("Adjoint identity <K x, z> == <x, K^T z>");
    const real_t lhs = y * z;
    const real_t rhs = x * adj;
    const real_t scale = std::max(real_t(1.0), std::abs(rhs));
    REQUIRE(std::abs(lhs - rhs) <= tol * scale);
  }
}

/// Verify large operators using:
/// 1) exact analytical rank-1 Kronecker sums (independent ground truth), and
/// 2) dense random adjointness + cross-backend differential check.
template <size_t N>
void VerifyLargeConfig(
    KroneckerProductOperator<N> &op, const std::array<DenseMatrix, N> &mats,
    typename KroneckerProductOperator<N>::ContractionBackend backend,
    bool use_eo, int seed = 201) {
  using Backend = typename KroneckerProductOperator<N>::ContractionBackend;
  const char *backend_name = BackendName<N>(backend);
  const int h = op.Height();
  const int w = op.Width();
  CAPTURE(N, h, w, backend_name, use_eo);

  // 1. Independent rank-1 Kronecker sum oracle.
  Vector x_r1(w), y_r1_ref(h), z_r1(h), w_r1_ref(w);
  x_r1 = real_t(0.0);
  y_r1_ref = real_t(0.0);
  z_r1 = real_t(0.0);
  w_r1_ref = real_t(0.0);

  int s = seed;
  for (int term = 0; term < 3; ++term) {
    std::array<Vector, N> vx, vy, vz, vw;
    for (size_t k = 0; k < N; ++k) {
      vx[k] = MakeRandomVector(mats[k].Width(), s++);
      vy[k].SetSize(mats[k].Height());
      mats[k].Mult(vx[k], vy[k]);

      vz[k] = MakeRandomVector(mats[k].Height(), s++);
      vw[k].SetSize(mats[k].Width());
      mats[k].MultTranspose(vz[k], vw[k]);
    }
    x_r1 += KronVecArray(vx);
    y_r1_ref += KronVecArray(vy);
    z_r1 += KronVecArray(vz);
    w_r1_ref += KronVecArray(vw);
  }

  // 2. Cross-backend dense baseline (BlockedCPP, even-odd disabled).
  const Vector x_dense = MakeRandomVector(w, s++);
  const Vector z_dense = MakeRandomVector(h, s++);

  op.SetBackend(Backend::BlockedCPP);
  op.SetUseEvenOdd(false);
  Vector y_base(h), w_base(w);
  y_base = real_t(42.0);
  w_base = real_t(-42.0);
  op.Mult(x_dense, y_base);
  op.MultTranspose(z_dense, w_base);

  // Evaluate requested configuration.
  op.SetBackend(backend);
  op.SetUseEvenOdd(use_eo);

  const real_t tol = ContractionTol(w, h);
  Vector y(h), adj(w);

  y = real_t(777.0);
  adj = real_t(-777.0);
  op.Mult(x_r1, y);
  op.MultTranspose(z_r1, adj);
  {
    INFO("Mult vs exact rank-1 Kronecker sum");
    REQUIRE(RelativeError(y, y_r1_ref) <= tol);
  }
  {
    INFO("MultTranspose vs exact rank-1 Kronecker sum");
    REQUIRE(RelativeError(adj, w_r1_ref) <= tol);
  }

  y = real_t(555.0);
  adj = real_t(-555.0);
  op.Mult(x_dense, y);
  op.MultTranspose(z_dense, adj);
  {
    INFO("Dense Mult cross-backend agreement");
    REQUIRE(RelativeError(y, y_base) <= tol);
  }
  {
    INFO("Dense MultTranspose cross-backend agreement");
    REQUIRE(RelativeError(adj, w_base) <= tol);
  }
  {
    INFO("Dense adjoint identity <K x, z> == <x, K^T z>");
    const real_t lhs = y * z_dense;
    const real_t rhs = x_dense * adj;
    const real_t scale = std::max(real_t(1.0), std::abs(rhs));
    REQUIRE(std::abs(lhs - rhs) <= tol * scale);
  }
}

} // namespace kronecker_test

TEST_CASE("KroneckerProductOperator hand-calculated ordering",
          "[Linalg][Kronecker]") {
  using namespace kronecker_test;
  using Backend = KroneckerProductOperator<2>::ContractionBackend;

  const auto backend =
      GENERATE(Backend::Auto, Backend::BLAS, Backend::BlockedCPP);
  const bool use_eo = GENERATE(false, true);
  const char *backend_name = BackendName<2>(backend);
  CAPTURE(backend_name, use_eo);

  // D1 = [[1, 2], [3, 4]], D2 = [[0, 5], [6, 7]]
  // (D2 \otimes D1)   * [1, 2, 3, 4]^T = [55, 125, 107, 241]^T
  // (D2 \otimes D1)^T * [1, 2, 3, 4]^T = [90, 132, 140, 204]^T
  DenseMatrix D1(2, 2), D2(2, 2);
  D1(0, 0) = real_t(1.0);
  D1(0, 1) = real_t(2.0);
  D1(1, 0) = real_t(3.0);
  D1(1, 1) = real_t(4.0);

  D2(0, 0) = real_t(0.0);
  D2(0, 1) = real_t(5.0);
  D2(1, 0) = real_t(6.0);
  D2(1, 1) = real_t(7.0);

  KroneckerProductOperator op(D1, D2);
  op.SetBackend(backend);
  op.SetUseEvenOdd(use_eo);

  REQUIRE(op.Height() == 4);
  REQUIRE(op.Width() == 4);

  const Vector x({real_t(1.0), real_t(2.0), real_t(3.0), real_t(4.0)});
  const Vector exp_fwd(
      {real_t(55.0), real_t(125.0), real_t(107.0), real_t(241.0)});
  const Vector exp_adj(
      {real_t(90.0), real_t(132.0), real_t(140.0), real_t(204.0)});
  const real_t tol = ContractionTol(4, 4);

  Vector y(4);
  y = real_t(999.0);
  op.Mult(x, y);
  REQUIRE(RelativeError(y, exp_fwd) <= tol);

  y = real_t(-999.0);
  op.MultTranspose(x, y);
  REQUIRE(RelativeError(y, exp_adj) <= tol);

  // Non-commutativity: swapping factors produces a distinct operator.
  KroneckerProductOperator op_rev(D2, D1);
  Vector y_rev(4);
  op_rev.Mult(x, y_rev);
  REQUIRE(RelativeError(y_rev, exp_fwd) > real_t(1e-2));
}

TEST_CASE("KroneckerProductOperator matches explicit DenseMatrix",
          "[Linalg][Kronecker]") {
  using namespace kronecker_test;

  const int b_idx = GENERATE(0, 1, 2);
  const bool use_eo = GENERATE(false, true);

  auto get_backend = [&](auto tag) {
    constexpr size_t D = decltype(tag)::value;
    using B = typename KroneckerProductOperator<D>::ContractionBackend;
    return b_idx == 0 ? B::Auto : (b_idx == 1 ? B::BLAS : B::BlockedCPP);
  };

  SECTION("N = 1 (square, rectangular, and degenerate)") {
    const auto backend = get_backend(std::integral_constant<size_t, 1>{});
    const std::vector<std::pair<int, int>> shapes = {
        {1, 1}, {4, 4}, {5, 2}, {2, 5}, {1, 7}, {7, 1}, {13, 9}};
    int seed = 100;
    for (const auto &sh : shapes) {
      CAPTURE(sh.first, sh.second);
      const DenseMatrix D0 = MakeRandomMatrix(sh.first, sh.second, seed++);
      KroneckerProductOperator op(D0);
      VerifyExplicitConfig(op, D0, backend, use_eo, seed++);
    }
  }

  SECTION("N = 2 (square, tall, wide, mixed, and aggressive reordering)") {
    const auto backend = get_backend(std::integral_constant<size_t, 2>{});
    struct Case2 {
      int m0, n0, m1, n1;
    };
    const std::vector<Case2> cases = {
        {3, 3, 4, 4},  {4, 2, 5, 3},  {2, 4, 3, 5}, {5, 2, 2, 6},
        {2, 9, 11, 3}, {11, 3, 2, 9}, {1, 7, 9, 1}};
    int seed = 200;
    for (const Case2 &c : cases) {
      CAPTURE(c.m0, c.n0, c.m1, c.n1);
      const DenseMatrix D0 = MakeRandomMatrix(c.m0, c.n0, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(c.m1, c.n1, seed++);
      KroneckerProductOperator op(D0, D1);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1), backend, use_eo,
                           seed++);
    }
  }

  SECTION("N = 3 (interior modes and degenerate dimensions)") {
    const auto backend = get_backend(std::integral_constant<size_t, 3>{});
    struct Case3 {
      int m0, n0, m1, n1, m2, n2;
    };
    const std::vector<Case3> cases = {{1, 1, 1, 1, 1, 1}, {1, 4, 3, 1, 2, 2},
                                      {3, 2, 1, 1, 4, 3}, {3, 3, 4, 4, 5, 5},
                                      {4, 2, 3, 5, 2, 3}, {2, 5, 6, 2, 3, 4},
                                      {3, 7, 8, 4, 5, 2}, {7, 3, 4, 8, 2, 5}};
    int seed = 300;
    for (const Case3 &c : cases) {
      CAPTURE(c.m0, c.n0, c.m1, c.n1, c.m2, c.n2);
      const DenseMatrix D0 = MakeRandomMatrix(c.m0, c.n0, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(c.m1, c.n1, seed++);
      const DenseMatrix D2 = MakeRandomMatrix(c.m2, c.n2, seed++);
      KroneckerProductOperator op(D0, D1, D2);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1, D2), backend, use_eo,
                           seed++);
    }
  }

  SECTION("N = 4 (multi-swap ping-pong and asymmetric reordering)") {
    const auto backend = get_backend(std::integral_constant<size_t, 4>{});
    int seed = 400;
    {
      const DenseMatrix D0 = MakeRandomMatrix(2, 3, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(3, 2, seed++);
      const DenseMatrix D2 = MakeRandomMatrix(4, 2, seed++);
      const DenseMatrix D3 = MakeRandomMatrix(2, 3, seed++);
      KroneckerProductOperator op(D0, D1, D2, D3);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1, D2, D3), backend,
                           use_eo, seed++);
    }
    {
      const DenseMatrix D0 = MakeRandomMatrix(6, 2, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(2, 5, seed++);
      const DenseMatrix D2 = MakeRandomMatrix(3, 3, seed++);
      const DenseMatrix D3 = MakeRandomMatrix(1, 4, seed++);
      KroneckerProductOperator op(D0, D1, D2, D3);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1, D2, D3), backend,
                           use_eo, seed++);
    }
  }

  SECTION("N = 5 and N = 6 (deep ping-pong buffer alternation)") {
    int seed = 500;
    {
      const auto backend = get_backend(std::integral_constant<size_t, 5>{});
      const DenseMatrix D0 = MakeRandomMatrix(2, 3, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(3, 2, seed++);
      const DenseMatrix D2 = MakeRandomMatrix(2, 2, seed++);
      const DenseMatrix D3 = MakeRandomMatrix(4, 2, seed++);
      const DenseMatrix D4 = MakeRandomMatrix(2, 3, seed++);
      KroneckerProductOperator op(D0, D1, D2, D3, D4);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1, D2, D3, D4), backend,
                           use_eo, seed++);
    }
    {
      const auto backend = get_backend(std::integral_constant<size_t, 6>{});
      const DenseMatrix D0 = MakeRandomMatrix(2, 2, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(3, 2, seed++);
      const DenseMatrix D2 = MakeRandomMatrix(2, 3, seed++);
      const DenseMatrix D3 = MakeRandomMatrix(2, 2, seed++);
      const DenseMatrix D4 = MakeRandomMatrix(3, 2, seed++);
      const DenseMatrix D5 = MakeRandomMatrix(2, 3, seed++);
      KroneckerProductOperator op(D0, D1, D2, D3, D4, D5);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1, D2, D3, D4, D5),
                           backend, use_eo, seed++);
    }
  }
}

TEST_CASE("KroneckerProductOperator power-of-two boundary stress",
          "[Linalg][Kronecker]") {
  using namespace kronecker_test;

  const int b_idx = GENERATE(0, 1, 2);
  const bool use_eo = GENERATE(false, true);

  SECTION("1D boundary neighborhoods around 4, 8, 32, and 64") {
    using B1 = KroneckerProductOperator<1>::ContractionBackend;
    const B1 backend =
        b_idx == 0 ? B1::Auto : (b_idx == 1 ? B1::BLAS : B1::BlockedCPP);

    const std::vector<std::pair<int, int>> shapes = {
        {3, 35},  {4, 32},  {5, 33},  {7, 8},   {8, 9},
        {31, 63}, {32, 64}, {33, 65}, {65, 33}, {71, 37}};
    int seed = 600;
    for (const auto &sh : shapes) {
      CAPTURE(sh.first, sh.second);
      const DenseMatrix D0 = MakeRandomMatrix(sh.first, sh.second, seed++);
      KroneckerProductOperator op(D0);
      VerifyExplicitConfig(op, D0, backend, use_eo, seed++);
    }
  }

  SECTION("2D boundary neighborhoods around 32, 64, and 128") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    const B2 backend =
        b_idx == 0 ? B2::Auto : (b_idx == 1 ? B2::BLAS : B2::BlockedCPP);

    struct Case2 {
      int m0, n0, m1, n1;
    };
    const std::vector<Case2> cases = {{63, 31, 15, 17},
                                      {64, 32, 32, 32},
                                      {65, 33, 33, 35},
                                      {69, 37, 19, 23},
                                      {9, 11, 127, 129}};
    int seed = 700;
    for (const Case2 &c : cases) {
      CAPTURE(c.m0, c.n0, c.m1, c.n1);
      const DenseMatrix D0 = MakeRandomMatrix(c.m0, c.n0, seed++);
      const DenseMatrix D1 = MakeRandomMatrix(c.m1, c.n1, seed++);
      KroneckerProductOperator op(D0, D1);
      VerifyLargeConfig<2>(op, {D0, D1}, backend, use_eo, seed++);
    }
  }

  SECTION("3D large rectangular factors") {
    using B3 = KroneckerProductOperator<3>::ContractionBackend;
    const B3 backend =
        b_idx == 0 ? B3::Auto : (b_idx == 1 ? B3::BLAS : B3::BlockedCPP);

    const DenseMatrix D0 = MakeRandomMatrix(17, 13, 801);
    const DenseMatrix D1 = MakeRandomMatrix(19, 15, 802);
    const DenseMatrix D2 = MakeRandomMatrix(11, 14, 803);
    KroneckerProductOperator op(D0, D1, D2);
    VerifyLargeConfig<3>(op, {D0, D1, D2}, backend, use_eo, 804);
  }
}

TEST_CASE("KroneckerProductOperator edge cases, aliasing, and statefulness",
          "[Linalg][Kronecker]") {
  using namespace kronecker_test;

  SECTION("Identity factors preserve input and trigger symmetry for n >= 8") {
    using B3 = KroneckerProductOperator<3>::ContractionBackend;
    using Sym3 = KroneckerProductOperator<3>::SymmetryType;

    const auto backend = GENERATE(B3::Auto, B3::BLAS, B3::BlockedCPP);
    const bool use_eo = GENERATE(false, true);
    const char *backend_name = BackendName<3>(backend);
    CAPTURE(backend_name, use_eo);

    KroneckerProductOperator op(MakeIdentity(3), MakeIdentity(9),
                                MakeIdentity(10));
    REQUIRE(op.GetSymmetryType(0, false) == Sym3::None);
    REQUIRE(op.GetSymmetryType(1, false) == Sym3::Centrosymmetric);
    REQUIRE(op.GetSymmetryType(2, false) == Sym3::Centrosymmetric);

    op.SetBackend(backend);
    op.SetUseEvenOdd(use_eo);

    const Vector x = MakeRandomVector(op.Width(), 901);
    Vector y(op.Height()), w(op.Width());
    y = real_t(123.0);
    w = real_t(-123.0);

    op.Mult(x, y);
    op.MultTranspose(x, w);

    const real_t tol = ContractionTol(op.Width(), op.Height());
    REQUIRE(RelativeError(y, x) <= tol);
    REQUIRE(RelativeError(w, x) <= tol);
  }

  SECTION("Zero matrix factor produces exact zero and no symmetry trigger") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;

    // Small zero factor inside 3D product.
    {
      DenseMatrix Z(3, 4);
      Z = real_t(0.0);
      KroneckerProductOperator op(MakeRandomMatrix(2, 3, 910), Z,
                                  MakeRandomMatrix(4, 2, 911));
      const Vector x = MakeRandomVector(op.Width(), 912);
      const Vector z = MakeRandomVector(op.Height(), 913);
      Vector y(op.Height()), w(op.Width());
      y = real_t(42.0);
      w = real_t(-42.0);
      op.Mult(x, y);
      op.MultTranspose(z, w);
      REQUIRE(y.Norml2() == real_t(0.0));
      REQUIRE(w.Norml2() == real_t(0.0));
    }

    // Zero factor >= 8x8 exercises max_abs == 0 early exit in BuildEvenOdd.
    {
      const auto backend = GENERATE(B2::Auto, B2::BLAS, B2::BlockedCPP);
      const bool use_eo = GENERATE(false, true);

      DenseMatrix Z(10, 11);
      Z = real_t(0.0);
      const DenseMatrix D1 = MakeCentrosymmetricMatrix(9, 9, 914);
      KroneckerProductOperator op(Z, D1);

      REQUIRE(op.GetSymmetryType(0, false) == Sym2::None);
      REQUIRE(op.GetSymmetryType(0, true) == Sym2::None);
      REQUIRE(op.GetSymmetryType(1, false) == Sym2::Centrosymmetric);

      VerifyExplicitConfig(op, ExplicitKronecker(Z, D1), backend, use_eo, 915);
    }
  }

  SECTION("Deep copy: mutating source matrices does not affect operator") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    const auto backend = GENERATE(B2::Auto, B2::BLAS, B2::BlockedCPP);
    const bool use_eo = GENERATE(false, true);

    DenseMatrix D0 = MakeCentrosymmetricMatrix(10, 10, 920);
    DenseMatrix D1 = MakeRowEvenOddMatrix(9, 10, 921, 0);
    const DenseMatrix K_orig = ExplicitKronecker(D0, D1);

    KroneckerProductOperator op(D0, D1);

    D0 = real_t(999.0);
    D1 = real_t(-888.0);

    VerifyExplicitConfig(op, K_orig, backend, use_eo, 922);
  }

  SECTION("Repeated interleaved Mult and MultTranspose (odd and even N)") {
    const DenseMatrix D1 = MakeRandomMatrix(3, 4, 930);
    const DenseMatrix D2 = MakeRandomMatrix(5, 2, 931);
    const DenseMatrix D3 = MakeRandomMatrix(2, 3, 932);

    KroneckerProductOperator op_even(D1, D2);
    const DenseMatrix K_even = ExplicitKronecker(D1, D2);

    KroneckerProductOperator op_odd(D1, D2, D3);
    const DenseMatrix K_odd = ExplicitKronecker(D1, D2, D3);

    const real_t tol_e = ContractionTol(op_even.Width(), op_even.Height());
    const real_t tol_o = ContractionTol(op_odd.Width(), op_odd.Height());

    for (int call = 0; call < 5; ++call) {
      CAPTURE(call);
      const Vector xe = MakeRandomVector(op_even.Width(), 940 + 4 * call);
      const Vector ze = MakeRandomVector(op_even.Height(), 941 + 4 * call);
      Vector ye(op_even.Height()), ref_ye(op_even.Height());
      Vector we(op_even.Width()), ref_we(op_even.Width());
      op_even.Mult(xe, ye);
      K_even.Mult(xe, ref_ye);
      op_even.MultTranspose(ze, we);
      K_even.MultTranspose(ze, ref_we);
      REQUIRE(RelativeError(ye, ref_ye) <= tol_e);
      REQUIRE(RelativeError(we, ref_we) <= tol_e);

      const Vector xo = MakeRandomVector(op_odd.Width(), 942 + 4 * call);
      const Vector zo = MakeRandomVector(op_odd.Height(), 943 + 4 * call);
      Vector yo(op_odd.Height()), ref_yo(op_odd.Height());
      Vector wo(op_odd.Width()), ref_wo(op_odd.Width());
      op_odd.Mult(xo, yo);
      K_odd.Mult(xo, ref_yo);
      op_odd.MultTranspose(zo, wo);
      K_odd.MultTranspose(zo, ref_wo);
      REQUIRE(RelativeError(yo, ref_yo) <= tol_o);
      REQUIRE(RelativeError(wo, ref_wo) <= tol_o);
    }
  }

  SECTION("Const instance produces bit-identical repeated results") {
    const DenseMatrix D0 = MakeCentrosymmetricMatrix(9, 10, 960);
    const DenseMatrix D1 = MakeRowEvenOddMatrix(10, 8, 961, 0);
    const DenseMatrix D2 = MakeRandomMatrix(7, 9, 962);
    const KroneckerProductOperator op(D0, D1, D2);

    const Vector x1 = MakeRandomVector(op.Width(), 963);
    const Vector x2 = MakeRandomVector(op.Width(), 964);
    const Vector z1 = MakeRandomVector(op.Height(), 965);
    const Vector z2 = MakeRandomVector(op.Height(), 966);

    Vector y1_a(op.Height()), y1_b(op.Height()), y2(op.Height());
    Vector w1_a(op.Width()), w1_b(op.Width()), w2(op.Width());

    op.Mult(x1, y1_a);
    op.MultTranspose(z1, w1_a);
    op.Mult(x2, y2);
    op.MultTranspose(z2, w2);

    y1_b = real_t(-314.15);
    w1_b = real_t(271.82);
    op.Mult(x1, y1_b);
    op.MultTranspose(z1, w1_b);

    REQUIRE(RelativeError(y1_b, y1_a) == real_t(0.0));
    REQUIRE(RelativeError(w1_b, w1_a) == real_t(0.0));
  }

  SECTION("In-place aliasing (&x == &y) for N = 1, 2, 3, 4") {
    {
      using B1 = KroneckerProductOperator<1>::ContractionBackend;
      const auto backend = GENERATE(B1::Auto, B1::BLAS, B1::BlockedCPP);
      CAPTURE(BackendName<1>(backend));

      const DenseMatrix D0 = MakeRandomMatrix(19, 19, 970);
      KroneckerProductOperator op(D0);
      op.SetBackend(backend);

      const Vector x = MakeRandomVector(19, 971);
      Vector y_ref(19), w_ref(19);
      D0.Mult(x, y_ref);
      D0.MultTranspose(x, w_ref);

      Vector v(x), w(x);
      op.Mult(v, v);
      op.MultTranspose(w, w);
      const real_t tol = ContractionTol(19, 19);
      REQUIRE(RelativeError(v, y_ref) <= tol);
      REQUIRE(RelativeError(w, w_ref) <= tol);
    }

    {
      const DenseMatrix D0 = MakeRandomMatrix(17, 17, 972);
      const DenseMatrix D1 = MakeCentrosymmetricMatrix(13, 13, 973);
      KroneckerProductOperator op(D0, D1);
      const Vector x = MakeRandomVector(op.Width(), 974);
      Vector y_ref(op.Height()), w_ref(op.Width());
      op.Mult(x, y_ref);
      op.MultTranspose(x, w_ref);

      Vector v(x), w(x);
      op.Mult(v, v);
      op.MultTranspose(w, w);
      const real_t tol = ContractionTol(op.Width(), op.Height());
      REQUIRE(RelativeError(v, y_ref) <= tol);
      REQUIRE(RelativeError(w, w_ref) <= tol);
    }

    {
      const DenseMatrix D0 = MakeRandomMatrix(9, 9, 975);
      const DenseMatrix D1 = MakeRandomMatrix(10, 10, 976);
      const DenseMatrix D2 = MakeCentrosymmetricMatrix(11, 11, 977);
      KroneckerProductOperator op(D0, D1, D2);
      const Vector x = MakeRandomVector(op.Width(), 978);
      Vector y_ref(op.Height()), w_ref(op.Width());
      op.Mult(x, y_ref);
      op.MultTranspose(x, w_ref);

      Vector v(x), w(x);
      op.Mult(v, v);
      op.MultTranspose(w, w);
      const real_t tol = ContractionTol(op.Width(), op.Height());
      REQUIRE(RelativeError(v, y_ref) <= tol);
      REQUIRE(RelativeError(w, w_ref) <= tol);
    }

    {
      const DenseMatrix D0 = MakeRandomMatrix(3, 3, 979);
      const DenseMatrix D1 = MakeRandomMatrix(4, 4, 980);
      const DenseMatrix D2 = MakeRandomMatrix(3, 3, 981);
      const DenseMatrix D3 = MakeRandomMatrix(2, 2, 982);
      KroneckerProductOperator op(D0, D1, D2, D3);
      const Vector x = MakeRandomVector(op.Width(), 983);
      Vector y_ref(op.Height()), w_ref(op.Width());
      op.Mult(x, y_ref);
      op.MultTranspose(x, w_ref);

      Vector v(x), w(x);
      op.Mult(v, v);
      op.MultTranspose(w, w);
      const real_t tol = ContractionTol(op.Width(), op.Height());
      REQUIRE(RelativeError(v, y_ref) <= tol);
      REQUIRE(RelativeError(w, w_ref) <= tol);
    }
  }

  SECTION("Partial overlapping vector views (N = 1 and N = 2)") {
    {
      const int m = 13, n = 11;
      const DenseMatrix D0 = MakeRandomMatrix(m, n, 985);
      KroneckerProductOperator op(D0);
      const Vector x_orig = MakeRandomVector(n, 986);
      const Vector z_orig = MakeRandomVector(m, 987);
      Vector y_ref(m), w_ref(n);
      D0.Mult(x_orig, y_ref);
      D0.MultTranspose(z_orig, w_ref);

      Vector backing(m + 4);
      Vector x_view(backing.GetData(), n);
      Vector y_view(backing.GetData() + 2, m);
      x_view = x_orig;
      op.Mult(x_view, y_view);
      REQUIRE(RelativeError(y_view, y_ref) <= ContractionTol(n, m));

      Vector z_view(backing.GetData() + 2, m);
      Vector w_view(backing.GetData(), n);
      z_view = z_orig;
      op.MultTranspose(z_view, w_view);
      REQUIRE(RelativeError(w_view, w_ref) <= ContractionTol(n, m));
    }

    {
      const DenseMatrix D0 = MakeRandomMatrix(5, 4, 988);
      const DenseMatrix D1 = MakeRandomMatrix(6, 5, 989);
      KroneckerProductOperator op(D0, D1);
      const int w_sz = op.Width();
      const int h_sz = op.Height();
      const Vector x_orig = MakeRandomVector(w_sz, 990);
      const Vector z_orig = MakeRandomVector(h_sz, 991);
      Vector y_ref(h_sz), w_ref(w_sz);
      op.Mult(x_orig, y_ref);
      op.MultTranspose(z_orig, w_ref);

      Vector backing(std::max(w_sz, h_sz) + 7);
      Vector x_view(backing.GetData(), w_sz);
      Vector y_view(backing.GetData() + 5, h_sz);
      x_view = x_orig;
      op.Mult(x_view, y_view);
      REQUIRE(RelativeError(y_view, y_ref) <= ContractionTol(w_sz, h_sz));

      Vector z_view(backing.GetData() + 5, h_sz);
      Vector w_view(backing.GetData(), w_sz);
      z_view = z_orig;
      op.MultTranspose(z_view, w_view);
      REQUIRE(RelativeError(w_view, w_ref) <= ContractionTol(w_sz, h_sz));
    }
  }
}

TEST_CASE("KroneckerProductOperator algebraic and spectral invariants",
          "[Linalg][Kronecker]") {
  using namespace kronecker_test;

  SECTION("Separable vector action for Mult and MultTranspose") {
    const DenseMatrix D1 = MakeRandomMatrix(10, 8, 1001);
    const DenseMatrix D2 = MakeRandomMatrix(7, 9, 1002);
    const DenseMatrix D3 = MakeRandomMatrix(6, 11, 1003);

    const Vector v1 = MakeRandomVector(8, 1004);
    const Vector v2 = MakeRandomVector(9, 1005);
    const Vector v3 = MakeRandomVector(11, 1006);
    const Vector x = KronVec(v3, KronVec(v2, v1));

    Vector w1(10), w2(7), w3(6);
    D1.Mult(v1, w1);
    D2.Mult(v2, w2);
    D3.Mult(v3, w3);
    const Vector exp_fwd = KronVec(w3, KronVec(w2, w1));

    const Vector u1 = MakeRandomVector(10, 1007);
    const Vector u2 = MakeRandomVector(7, 1008);
    const Vector u3 = MakeRandomVector(6, 1009);
    const Vector z = KronVec(u3, KronVec(u2, u1));

    Vector q1(8), q2(9), q3(11);
    D1.MultTranspose(u1, q1);
    D2.MultTranspose(u2, q2);
    D3.MultTranspose(u3, q3);
    const Vector exp_adj = KronVec(q3, KronVec(q2, q1));

    KroneckerProductOperator op(D1, D2, D3);
    Vector y(op.Height()), w(op.Width());
    op.Mult(x, y);
    op.MultTranspose(z, w);

    const real_t tol = ContractionTol(op.Width(), op.Height());
    REQUIRE(RelativeError(y, exp_fwd) <= tol);
    REQUIRE(RelativeError(w, exp_adj) <= tol);
  }

  SECTION("Mixed-product property for Mult and MultTranspose") {
    const DenseMatrix A1 = MakeRandomMatrix(9, 8, 1010);
    const DenseMatrix B1 = MakeRandomMatrix(8, 10, 1011);
    const DenseMatrix A2 = MakeRandomMatrix(8, 9, 1012);
    const DenseMatrix B2 = MakeRandomMatrix(9, 7, 1013);
    const DenseMatrix A3 = MakeRandomMatrix(10, 8, 1014);
    const DenseMatrix B3 = MakeRandomMatrix(8, 9, 1015);

    DenseMatrix C1(9, 10), C2(8, 7), C3(10, 9);
    mfem::Mult(A1, B1, C1);
    mfem::Mult(A2, B2, C2);
    mfem::Mult(A3, B3, C3);

    KroneckerProductOperator opA(A1, A2, A3);
    KroneckerProductOperator opB(B1, B2, B3);
    KroneckerProductOperator opC(C1, C2, C3);
    const real_t tol = ContractionTol(opC.Width(), opC.Height());

    const Vector x = MakeRandomVector(opB.Width(), 1016);
    Vector tmp(opB.Height()), y_comp(opA.Height()), y_dir(opC.Height());
    opB.Mult(x, tmp);
    opA.Mult(tmp, y_comp);
    opC.Mult(x, y_dir);
    REQUIRE(RelativeError(y_comp, y_dir) <= tol);

    const Vector z = MakeRandomVector(opA.Height(), 1017);
    Vector tmp_t(opA.Width()), w_comp(opB.Width()), w_dir(opC.Width());
    opA.MultTranspose(z, tmp_t);
    opB.MultTranspose(tmp_t, w_comp);
    opC.MultTranspose(z, w_dir);
    REQUIRE(RelativeError(w_comp, w_dir) <= tol);
  }

  SECTION("Linearity for Mult and MultTranspose") {
    KroneckerProductOperator op(MakeRandomMatrix(5, 4, 1020),
                                MakeRandomMatrix(3, 6, 1021),
                                MakeRandomMatrix(4, 3, 1022));
    const real_t tol = ContractionTol(op.Width(), op.Height());
    const real_t alpha = real_t(2.75);
    const real_t beta = real_t(-1.5);

    const Vector x = MakeRandomVector(op.Width(), 1023);
    const Vector z = MakeRandomVector(op.Width(), 1024);
    Vector xz(x.Size());
    xz = real_t(0.0);
    xz.Add(alpha, x);
    xz.Add(beta, z);

    Vector y_lhs(op.Height()), yx(op.Height()), yz(op.Height());
    op.Mult(xz, y_lhs);
    op.Mult(x, yx);
    op.Mult(z, yz);
    Vector y_rhs(op.Height());
    y_rhs = real_t(0.0);
    y_rhs.Add(alpha, yx);
    y_rhs.Add(beta, yz);
    REQUIRE(RelativeError(y_lhs, y_rhs) <= tol);

    const Vector u = MakeRandomVector(op.Height(), 1025);
    const Vector v = MakeRandomVector(op.Height(), 1026);
    Vector uv(u.Size());
    uv = real_t(0.0);
    uv.Add(alpha, u);
    uv.Add(beta, v);

    Vector w_lhs(op.Width()), wu(op.Width()), wv(op.Width());
    op.MultTranspose(uv, w_lhs);
    op.MultTranspose(u, wu);
    op.MultTranspose(v, wv);
    Vector w_rhs(op.Width());
    w_rhs = real_t(0.0);
    w_rhs.Add(alpha, wu);
    w_rhs.Add(beta, wv);
    REQUIRE(RelativeError(w_lhs, w_rhs) <= tol);
  }

  SECTION("Transpose-of-factors equivalence via base const Operator &") {
    const DenseMatrix D0 = MakeRowEvenOddMatrix(9, 11, 1030, 0);
    const DenseMatrix D1 = MakeCentrosymmetricMatrix(10, 8, 1031);
    const DenseMatrix D2 = MakeRandomMatrix(7, 9, 1032);

    DenseMatrix D0_t, D1_t, D2_t;
    D0_t.Transpose(D0);
    D1_t.Transpose(D1);
    D2_t.Transpose(D2);

    const KroneckerProductOperator op(D0, D1, D2);
    const KroneckerProductOperator op_t(D0_t, D1_t, D2_t);
    const Operator &base_op = op;
    const Operator &base_op_t = op_t;

    REQUIRE(base_op.Height() == base_op_t.Width());
    REQUIRE(base_op.Width() == base_op_t.Height());

    const real_t tol = ContractionTol(base_op.Width(), base_op.Height());
    const Vector x = MakeRandomVector(base_op.Width(), 1033);
    const Vector z = MakeRandomVector(base_op.Height(), 1034);

    Vector y1(base_op.Height()), y2(base_op.Height());
    base_op.Mult(x, y1);
    base_op_t.MultTranspose(x, y2);
    REQUIRE(RelativeError(y1, y2) <= tol);

    Vector w1(base_op.Width()), w2(base_op.Width());
    base_op.MultTranspose(z, w1);
    base_op_t.Mult(z, w2);
    REQUIRE(RelativeError(w1, w2) <= tol);
  }

  SECTION("2D orthogonal DST-I eigenvectors V^T V = I and V V^T = I") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;

    const auto backend = GENERATE(B2::Auto, B2::BLAS, B2::BlockedCPP);
    const bool use_eo = GENERATE(false, true);
    CAPTURE(BackendName<2>(backend), use_eo);

    KroneckerProductOperator op(MakeDSTMatrix(16), MakeDSTMatrix(17));
    REQUIRE(op.GetSymmetryType(0, false) == Sym2::RowEvenOdd);
    REQUIRE(op.GetSymmetryType(1, false) == Sym2::RowEvenOdd);

    op.SetBackend(backend);
    op.SetUseEvenOdd(use_eo);

    const Vector x = MakeRandomVector(op.Width(), 1040);
    Vector tmp(op.Height()), y(op.Width());
    const real_t tol = ContractionTol(op.Width(), op.Height());

    op.Mult(x, tmp);
    op.MultTranspose(tmp, y);
    REQUIRE(RelativeError(y, x) <= tol);

    op.MultTranspose(x, tmp);
    op.Mult(tmp, y);
    REQUIRE(RelativeError(y, x) <= tol);
  }

  SECTION("3D column-scaled DST-I (strict ColEvenOdd / RowEvenOdd)") {
    using B3 = KroneckerProductOperator<3>::ContractionBackend;
    using Sym3 = KroneckerProductOperator<3>::SymmetryType;

    const auto backend = GENERATE(B3::Auto, B3::BLAS, B3::BlockedCPP);
    const bool use_eo = GENERATE(false, true);
    CAPTURE(BackendName<3>(backend), use_eo);

    auto make_scaled_dst = [](int n, Vector &inv_d) {
      DenseMatrix W = MakeDSTMatrix(n);
      inv_d.SetSize(n);
      for (int j = 0; j < n; ++j) {
        const real_t d = real_t(1.0) + real_t(0.15) * static_cast<real_t>(j);
        inv_d(j) = real_t(1.0) / (d * d);
        for (int i = 0; i < n; ++i) {
          W(i, j) *= d;
        }
      }
      return W;
    };

    Vector inv0, inv1, inv2;
    const DenseMatrix W0 = make_scaled_dst(12, inv0);
    const DenseMatrix W1 = make_scaled_dst(13, inv1);
    const DenseMatrix W2 = make_scaled_dst(14, inv2);
    const Vector inv_diag = KronVecArray<3>({inv0, inv1, inv2});

    KroneckerProductOperator op(W0, W1, W2);
    for (size_t k = 0; k < 3; ++k) {
      REQUIRE(op.GetSymmetryType(k, false) == Sym3::ColEvenOdd);
      REQUIRE(op.GetSymmetryType(k, true) == Sym3::RowEvenOdd);
    }

    op.SetBackend(backend);
    op.SetUseEvenOdd(use_eo);

    // W = V D => W D^{-2} W^T = V V^T = I.
    const Vector x = MakeRandomVector(op.Width(), 1050);
    Vector tmp(op.Height()), y(op.Width());
    op.MultTranspose(x, tmp);
    for (int i = 0; i < tmp.Size(); ++i) {
      tmp(i) *= inv_diag(i);
    }
    op.Mult(tmp, y);
    REQUIRE(RelativeError(y, x) <= ContractionTol(op.Width(), op.Height()));
  }
}

TEST_CASE("KroneckerProductOperator even-odd symmetry detection and kernels",
          "[Linalg][Kronecker]") {
  using namespace kronecker_test;

  const int b_idx = GENERATE(0, 1, 2);
  const bool use_eo = GENERATE(false, true);

  SECTION("RowEvenOdd and ColEvenOdd mixed, all-even, and all-odd") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;
    const B2 backend2 =
        b_idx == 0 ? B2::Auto : (b_idx == 1 ? B2::BLAS : B2::BlockedCPP);

    const int mode = GENERATE(0, 1, 2);
    CAPTURE(mode);

    const DenseMatrix D0 = MakeRowEvenOddMatrix(10, 11, 1100 + mode, mode);
    const DenseMatrix D1 = MakeColEvenOddMatrix(11, 9, 1110 + mode, mode);
    KroneckerProductOperator op(D0, D1);

    REQUIRE(op.GetSymmetryType(0, false) == Sym2::RowEvenOdd);
    REQUIRE(op.GetSymmetryType(0, true) == Sym2::ColEvenOdd);
    REQUIRE(op.GetSymmetryType(1, false) == Sym2::ColEvenOdd);
    REQUIRE(op.GetSymmetryType(1, true) == Sym2::RowEvenOdd);

    VerifyExplicitConfig(op, ExplicitKronecker(D0, D1), backend2, use_eo,
                         1120 + mode);
  }

  SECTION("3D RowEvenOdd and ColEvenOdd crossing batch chunking (R > 128)") {
    using B3 = KroneckerProductOperator<3>::ContractionBackend;
    using Sym3 = KroneckerProductOperator<3>::SymmetryType;
    const B3 backend3 =
        b_idx == 0 ? B3::Auto : (b_idx == 1 ? B3::BLAS : B3::BlockedCPP);

    const int mode = GENERATE(0, 1, 2);
    CAPTURE(mode);

    const DenseMatrix D0 = MakeRowEvenOddMatrix(11, 13, 1200 + mode, mode);
    const DenseMatrix D1 = MakeColEvenOddMatrix(12, 14, 1210 + mode, mode);
    const DenseMatrix D2 = MakeRowEvenOddMatrix(13, 12, 1220 + mode, mode);
    KroneckerProductOperator op(D0, D1, D2);

    REQUIRE(op.GetSymmetryType(0, false) == Sym3::RowEvenOdd);
    REQUIRE(op.GetSymmetryType(1, false) == Sym3::ColEvenOdd);
    REQUIRE(op.GetSymmetryType(2, false) == Sym3::RowEvenOdd);
    REQUIRE(op.GetSymmetryType(0, true) == Sym3::ColEvenOdd);
    REQUIRE(op.GetSymmetryType(1, true) == Sym3::RowEvenOdd);
    REQUIRE(op.GetSymmetryType(2, true) == Sym3::ColEvenOdd);

    VerifyLargeConfig<3>(op, {D0, D1, D2}, backend3, use_eo, 1230 + mode);
  }

  SECTION("Centrosymmetric square, rectangular, and R=127/128/129 chunks") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;
    const B2 backend2 =
        b_idx == 0 ? B2::Auto : (b_idx == 1 ? B2::BLAS : B2::BlockedCPP);

    {
      const DenseMatrix D0 = MakeCentrosymmetricMatrix(10, 13, 1300);
      const DenseMatrix D1 = MakeCentrosymmetricMatrix(13, 10, 1301);
      KroneckerProductOperator op(D0, D1);
      REQUIRE(op.GetSymmetryType(0, false) == Sym2::Centrosymmetric);
      REQUIRE(op.GetSymmetryType(0, true) == Sym2::Centrosymmetric);
      REQUIRE(op.GetSymmetryType(1, false) == Sym2::Centrosymmetric);
      REQUIRE(op.GetSymmetryType(1, true) == Sym2::Centrosymmetric);
      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1), backend2, use_eo,
                           1302);
    }

    const int r_dim = GENERATE(127, 128, 129);
    CAPTURE(r_dim);
    {
      const DenseMatrix D0 = MakeCentrosymmetricMatrix(10, 11, 1310);
      const DenseMatrix D1 = MakeRandomMatrix(r_dim, r_dim, 1311 + r_dim);
      KroneckerProductOperator op(D0, D1);
      REQUIRE(op.GetSymmetryType(0, false) == Sym2::Centrosymmetric);
      REQUIRE(op.GetSymmetryType(1, false) == Sym2::None);
      VerifyLargeConfig<2>(op, {D0, D1}, backend2, use_eo, 1320 + r_dim);
    }
  }

  SECTION("Minimum dimension threshold (7 vs 8) and rectangular reordering") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;
    const B2 backend2 =
        b_idx == 0 ? B2::Auto : (b_idx == 1 ? B2::BLAS : B2::BlockedCPP);

    const std::vector<std::pair<int, int>> dims = {{7, 7}, {7, 8}, {8, 7},
                                                   {8, 8}, {8, 9}, {9, 8}};
    int seed = 1400;
    for (const auto &sh : dims) {
      CAPTURE(sh.first, sh.second);
      const DenseMatrix D0 =
          MakeCentrosymmetricMatrix(sh.first, sh.second, seed++);
      const DenseMatrix D1 =
          MakeRowEvenOddMatrix(sh.second, sh.first, seed++, 0);
      KroneckerProductOperator op(D0, D1);

      const bool active = (sh.first >= 8 && sh.second >= 8);
      REQUIRE(op.GetSymmetryType(0, false) ==
              (active ? Sym2::Centrosymmetric : Sym2::None));
      REQUIRE(op.GetSymmetryType(1, false) ==
              (active ? Sym2::RowEvenOdd : Sym2::None));

      VerifyExplicitConfig(op, ExplicitKronecker(D0, D1), backend2, use_eo,
                           seed++);
    }

    // Rectangular symmetric factors forcing non-trivial mode reordering.
    {
      using B3 = KroneckerProductOperator<3>::ContractionBackend;
      using Sym3 = KroneckerProductOperator<3>::SymmetryType;
      const B3 backend3 =
          b_idx == 0 ? B3::Auto : (b_idx == 1 ? B3::BLAS : B3::BlockedCPP);

      const DenseMatrix D0 = MakeRowEvenOddMatrix(16, 8, 1450, 0);
      const DenseMatrix D1 = MakeCentrosymmetricMatrix(8, 17, 1451);
      const DenseMatrix D2 = MakeColEvenOddMatrix(12, 12, 1452, 0);
      KroneckerProductOperator op(D0, D1, D2);

      REQUIRE(op.GetSymmetryType(0, false) == Sym3::RowEvenOdd);
      REQUIRE(op.GetSymmetryType(1, false) == Sym3::Centrosymmetric);
      REQUIRE(op.GetSymmetryType(2, false) == Sym3::ColEvenOdd);

      VerifyLargeConfig<3>(op, {D0, D1, D2}, backend3, use_eo, 1453);
    }
  }

  SECTION("False-positive rejection: perturbed and anti-centrosymmetric") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;
    const B2 backend2 =
        b_idx == 0 ? B2::Auto : (b_idx == 1 ? B2::BLAS : B2::BlockedCPP);

    const real_t perturb = real_t(1e5) * std::numeric_limits<real_t>::epsilon();

    DenseMatrix D_row = MakeRowEvenOddMatrix(10, 10, 1500, 0);
    D_row(3, 1) += perturb;

    DenseMatrix D_col = MakeColEvenOddMatrix(10, 10, 1501, 0);
    D_col(2, 4) -= perturb;

    DenseMatrix D_cen = MakeCentrosymmetricMatrix(11, 11, 1502);
    D_cen(0, 0) += perturb;

    const DenseMatrix D_anti = MakeAntiCentrosymmetricMatrix(10, 11, 1503);

    {
      KroneckerProductOperator op(D_row, D_col);
      REQUIRE(op.GetSymmetryType(0, false) == Sym2::None);
      REQUIRE(op.GetSymmetryType(0, true) == Sym2::None);
      REQUIRE(op.GetSymmetryType(1, false) == Sym2::None);
      REQUIRE(op.GetSymmetryType(1, true) == Sym2::None);
      VerifyExplicitConfig(op, ExplicitKronecker(D_row, D_col), backend2,
                           use_eo, 1504);
    }
    {
      KroneckerProductOperator op(D_cen, D_anti);
      REQUIRE(op.GetSymmetryType(0, false) == Sym2::None);
      REQUIRE(op.GetSymmetryType(0, true) == Sym2::None);
      REQUIRE(op.GetSymmetryType(1, false) == Sym2::None);
      REQUIRE(op.GetSymmetryType(1, true) == Sym2::None);
      VerifyExplicitConfig(op, ExplicitKronecker(D_cen, D_anti), backend2,
                           use_eo, 1505);
    }
  }

  SECTION("Scale invariance for very small and very large symmetric factors") {
    using B2 = KroneckerProductOperator<2>::ContractionBackend;
    using Sym2 = KroneckerProductOperator<2>::SymmetryType;
    const B2 backend2 =
        b_idx == 0 ? B2::Auto : (b_idx == 1 ? B2::BLAS : B2::BlockedCPP);

    const real_t small_scale =
        std::is_same_v<real_t, float> ? real_t(1e-6) : real_t(1e-30);
    const real_t large_scale =
        std::is_same_v<real_t, float> ? real_t(1e+6) : real_t(1e+30);
    const real_t s = GENERATE_COPY(small_scale, large_scale);
    CAPTURE(s, BackendName<2>(backend2), use_eo);

    DenseMatrix D0 = MakeCentrosymmetricMatrix(10, 10, 1600);
    DenseMatrix D1 = MakeRowEvenOddMatrix(11, 10, 1601, 0);
    D0 *= s;
    D1 *= s;

    KroneckerProductOperator op(D0, D1);
    REQUIRE(op.GetSymmetryType(0, false) == Sym2::Centrosymmetric);
    REQUIRE(op.GetSymmetryType(1, false) == Sym2::RowEvenOdd);

    op.SetBackend(backend2);
    op.SetUseEvenOdd(use_eo);

    const DenseMatrix K = ExplicitKronecker(D0, D1);
    const Vector x = MakeRandomVector(op.Width(), 1602);
    const Vector z = MakeRandomVector(op.Height(), 1603);

    Vector y_ref(op.Height()), w_ref(op.Width());
    K.Mult(x, y_ref);
    K.MultTranspose(z, w_ref);

    Vector y(op.Height()), w(op.Width());
    y = real_t(999.0);
    w = real_t(-999.0);
    op.Mult(x, y);
    op.MultTranspose(z, w);

    const real_t tol = ContractionTol(op.Width(), op.Height());
    REQUIRE(StrictRelativeError(y, y_ref) <= tol);
    REQUIRE(StrictRelativeError(w, w_ref) <= tol);
  }
}

} // namespace mfem
