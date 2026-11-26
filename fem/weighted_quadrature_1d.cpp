// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#include "../general/error.hpp"
#include "intrules.hpp"
#include "weighted_quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef MFEM_USE_LAPACK
extern "C" void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,
                        double *b, int *ldb, double *s, double *rcond,
                        int *rank, double *work, int *lwork, int *info);

extern "C" void sgelss_(int *m, int *n, int *nrhs, float *a, int *lda, float *b,
                        int *ldb, float *s, float *rcond, int *rank,
                        float *work, int *lwork, int *info);
#endif

namespace mfem {

namespace {

#ifdef MFEM_USE_LAPACK
void LapackGelss(int *m, int *n, int *nrhs, double *a, int *lda, double *b,
                 int *ldb, double *s, double *rcond, int *rank, double *work,
                 int *lwork, int *info) {
  dgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}

void LapackGelss(int *m, int *n, int *nrhs, float *a, int *lda, float *b,
                 int *ldb, float *s, float *rcond, int *rank, float *work,
                 int *lwork, int *info) {
  sgelss_(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info);
}

/// Minimum-norm least-squares solve via LAPACK SVD (dgelss_/sgelss_).
void SolveMinNormLS(int m, int n, int nrhs, Vector &A, Vector &B, Vector &S,
                    Vector &work) {
  int lda = std::max(1, m);
  int ldb = std::max(1, std::max(m, n));
  int rank = 0;
  int info = 0;
  real_t rcond = static_cast<real_t>(1e-12);

  const int min_mn = std::max(1, std::min(m, n));
  if (S.Size() < min_mn) {
    S.SetSize(min_mn);
  }

  int lwork = -1;
  real_t work_query = 0.0;
  LapackGelss(&m, &n, &nrhs, A.GetData(), &lda, B.GetData(), &ldb, S.GetData(),
              &rcond, &rank, &work_query, &lwork, &info);
  MFEM_VERIFY(info == 0,
              "LAPACK gelss workspace query failed in BuildWQWeights1D.");

  lwork = std::max(1, static_cast<int>(std::ceil(work_query)));
  if (work.Size() < lwork) {
    work.SetSize(lwork);
  }
  lwork = work.Size();

  LapackGelss(&m, &n, &nrhs, A.GetData(), &lda, B.GetData(), &ldb, S.GetData(),
              &rcond, &rank, work.GetData(), &lwork, &info);
  MFEM_VERIFY(info == 0, "LAPACK gelss solve failed in BuildWQWeights1D.");
}
#endif

/// Validate open KnotVector and extract non-empty span indices and lengths.
void ExtractNonEmptySpans(const KnotVector &kv, Array<int> &span_knot_idx,
                          Vector &span_lengths) {
  const int p = kv.GetOrder();
  const int n = kv.GetNCP();
  const int m = kv.Size();

  MFEM_VERIFY(p >= 0, "KnotVector order must be non-negative.");
  MFEM_VERIFY(n > 0, "KnotVector must have at least one control point.");
  MFEM_VERIFY(m == n + p + 1,
              "Expected open KnotVector with Size() == GetNCP() + "
              "GetOrder() + 1.");

  for (int k = 0; k < m - 1; k++) {
    MFEM_VERIFY(kv[k + 1] >= kv[k], "KnotVector knots must be non-decreasing.");
  }

  const real_t domain_len = kv[m - 1] - kv[0];
  MFEM_VERIFY(domain_len > 0.0,
              "KnotVector parametric domain must have positive length.");
  const real_t knot_tol = static_cast<real_t>(1e-14) * domain_len;

  for (int k = 1; k <= p; k++) {
    MFEM_VERIFY(std::abs(kv[k] - kv[0]) <= knot_tol,
                "KnotVector must be open (clamped) at the left endpoint.");
    MFEM_VERIFY(std::abs(kv[m - 1 - k] - kv[m - 1]) <= knot_tol,
                "KnotVector must be open (clamped) at the right endpoint.");
  }

  span_knot_idx.SetSize(0);
  for (int k = p; k < n; k++) {
    if (kv[k + 1] - kv[k] > knot_tol) {
      span_knot_idx.Append(k);
    }
  }

  const int num_spans = span_knot_idx.Size();
  MFEM_VERIFY(num_spans > 0, "KnotVector has no non-empty knot spans.");

  span_lengths.SetSize(num_spans);
  for (int s = 0; s < num_spans; s++) {
    const int k = span_knot_idx[s];
    span_lengths[s] = kv[k + 1] - kv[k];
  }
}

/// Evaluate the p + 1 non-vanishing B-spline values B_{k-p+r}(u) and exact
/// global parametric derivatives dB_{k-p+r}/du(u) on non-empty span k.
void EvalBSplineAndDeriv(const KnotVector &kv, int p, int k, real_t u,
                         Vector &N, Vector &dN, Vector &Nm1, Vector &left,
                         Vector &right) {
  if (p == 0) {
    N[0] = 1.0;
    dN[0] = 0.0;
    return;
  }

  N[0] = 1.0;
  for (int d = 1; d <= p; d++) {
    if (d == p) {
      for (int r = 0; r < p; r++) {
        Nm1[r] = N[r];
      }
    }

    left[d] = u - kv[k + 1 - d];
    right[d] = kv[k + d] - u;

    real_t saved = 0.0;
    for (int r = 0; r < d; r++) {
      const real_t denom = right[r + 1] + left[d - r];
      const real_t temp = N[r] / denom;
      N[r] = saved + right[r + 1] * temp;
      saved = left[d - r] * temp;
    }
    N[d] = saved;
  }

  const real_t p_real = static_cast<real_t>(p);
  for (int r = 0; r <= p; r++) {
    real_t deriv = 0.0;
    if (r > 0) {
      const real_t denom = kv[k + r] - kv[k - p + r];
      deriv += (p_real / denom) * Nm1[r - 1];
    }
    if (r < p) {
      const real_t denom = kv[k + r + 1] - kv[k - p + r + 1];
      deriv -= (p_real / denom) * Nm1[r];
    }
    dN[r] = deriv;
  }
}

} // namespace

void BuildWQRule1D(const KnotVector &kv, const WQOptions &opts,
                   WQRule1D &rule) {
  Array<int> span_knot_idx;
  Vector span_lengths;
  ExtractNonEmptySpans(kv, span_knot_idx, span_lengths);

  const int p = kv.GetOrder();
  const int num_dofs = kv.GetNCP();
  const int num_spans = span_knot_idx.Size();

  rule.order = p;
  rule.num_dofs = num_dofs;
  rule.num_spans = num_spans;
  rule.span_lengths = span_lengths;
  rule.span_q_begin.SetSize(num_spans);
  rule.span_q_end.SetSize(num_spans);

  const int bdr_pts =
      (opts.boundary_pts_per_span > 0) ? opts.boundary_pts_per_span : (p + 1);
  MFEM_VERIFY(bdr_pts >= 1, "WQ boundary points per span must be at least 1.");

  Array<int> span_npts(num_spans);
  int total_qpts = 0;
  for (int s = 0; s < num_spans; s++) {
    int n_s = 0;
    if (num_spans == 1) {
      n_s = p + 1;
    } else if (s == 0 || s == num_spans - 1) {
      n_s = bdr_pts;
    } else {
      const int left_mult = span_knot_idx[s] - span_knot_idx[s - 1];
      const int right_mult = span_knot_idx[s + 1] - span_knot_idx[s];
      const int max_adj_mult = std::max(left_mult, right_mult);
      n_s = std::min(p + 1,
                     std::max(opts.interior_pts_per_span, max_adj_mult + 1));
    }

    const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 2 * n_s - 1);
    const int n_ir = ir.GetNPoints();
    MFEM_VERIFY(n_ir >= n_s,
                "Unexpected Gauss-Legendre point count in BuildWQRule1D.");

    span_npts[s] = n_ir;
    rule.span_q_begin[s] = total_qpts;
    total_qpts += n_ir;
    rule.span_q_end[s] = total_qpts;
  }

  rule.num_qpts = total_qpts;
  rule.xi.SetSize(total_qpts);
  rule.span_of_qpt.SetSize(total_qpts);
  rule.local_xi.SetSize(total_qpts);

  std::vector<real_t> local_pts;
  for (int s = 0; s < num_spans; s++) {
    const int n_ir = span_npts[s];
    const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 2 * n_ir - 1);
    MFEM_VERIFY(ir.GetNPoints() == n_ir,
                "Inconsistent Gauss-Legendre rule size in BuildWQRule1D.");

    local_pts.resize(static_cast<std::size_t>(n_ir));
    for (int q = 0; q < n_ir; q++) {
      local_pts[static_cast<std::size_t>(q)] = ir.IntPoint(q).x;
    }
    std::sort(local_pts.begin(), local_pts.end());

    const int k = span_knot_idx[s];
    const real_t u_a = kv[k];
    const real_t h = span_lengths[s];
    const int q0 = rule.span_q_begin[s];

    for (int q = 0; q < n_ir; q++) {
      const real_t x_loc = local_pts[static_cast<std::size_t>(q)];
      const int gq = q0 + q;
      rule.span_of_qpt[gq] = s;
      rule.local_xi[gq] = x_loc;
      rule.xi[gq] = u_a + h * x_loc;
    }
  }

  rule.q_begin.SetSize(num_dofs);
  rule.q_end.SetSize(num_dofs);
  rule.j_begin.SetSize(num_dofs);
  rule.j_end.SetSize(num_dofs);
  rule.trial_offsets.SetSize(num_dofs + 1);

  const int *sk_begin = span_knot_idx.GetData();
  const int *sk_end = sk_begin + num_spans;

  rule.trial_offsets[0] = 0;
  for (int i = 0; i < num_dofs; i++) {
    const int sb =
        static_cast<int>(std::lower_bound(sk_begin, sk_end, i) - sk_begin);
    const int se =
        static_cast<int>(std::upper_bound(sk_begin, sk_end, i + p) - sk_begin);
    MFEM_VERIFY(sb < se, "Basis function has empty support in BuildWQRule1D.");

    rule.q_begin[i] = rule.span_q_begin[sb];
    rule.q_end[i] = rule.span_q_end[se - 1];
    rule.j_begin[i] = span_knot_idx[sb] - p;
    rule.j_end[i] = span_knot_idx[se - 1] + 1;
    rule.trial_offsets[i + 1] =
        rule.trial_offsets[i] + (rule.q_end[i] - rule.q_begin[i]);
  }

  const int total_entries = rule.trial_offsets[num_dofs];
  rule.B_vals.SetSize(total_entries);
  rule.D_vals.SetSize(total_entries);
  rule.B_vals = 0.0;
  rule.D_vals = 0.0;

  Vector N(p + 1), dN(p + 1), Nm1(std::max(1, p)), left(p + 1), right(p + 1);
  for (int s = 0; s < num_spans; s++) {
    const int k = span_knot_idx[s];
    const int qb = rule.span_q_begin[s];
    const int qe = rule.span_q_end[s];

    for (int q = qb; q < qe; q++) {
      EvalBSplineAndDeriv(kv, p, k, rule.xi[q], N, dN, Nm1, left, right);
      for (int r = 0; r <= p; r++) {
        const int j = k - p + r;
        const int pos = rule.trial_offsets[j] + (q - rule.q_begin[j]);
        rule.B_vals[pos] = N[r];
        rule.D_vals[pos] = dN[r];
      }
    }
  }
}

void BuildWQWeights1D(const KnotVector &kv, const WQRule1D &rule,
                      WQWeights1D &weights) {
#ifndef MFEM_USE_LAPACK
  static_cast<void>(kv);
  static_cast<void>(rule);
  static_cast<void>(weights);
  MFEM_ABORT("BuildWQWeights1D requires MFEM to be built with LAPACK "
             "(MFEM_USE_LAPACK).");
#else
  Array<int> span_knot_idx;
  Vector span_lengths;
  ExtractNonEmptySpans(kv, span_knot_idx, span_lengths);

  const int p = kv.GetOrder();
  const int num_dofs = kv.GetNCP();
  const int num_spans = span_knot_idx.Size();

  MFEM_VERIFY(rule.order == p && rule.num_dofs == num_dofs &&
                  rule.num_spans == num_spans,
              "WQRule1D metadata does not match KnotVector in "
              "BuildWQWeights1D.");
  MFEM_VERIFY(
      rule.trial_offsets.Size() == num_dofs + 1 &&
          rule.q_begin.Size() == num_dofs && rule.q_end.Size() == num_dofs &&
          rule.j_begin.Size() == num_dofs && rule.j_end.Size() == num_dofs,
      "Incomplete WQRule1D arrays in BuildWQWeights1D.");

  const int total_entries = rule.trial_offsets[num_dofs];
  MFEM_VERIFY(rule.B_vals.Size() == total_entries &&
                  rule.D_vals.Size() == total_entries,
              "WQRule1D trial table sizes mismatch in BuildWQWeights1D.");

  weights.offsets = rule.trial_offsets;
  weights.W00.SetSize(total_entries);
  weights.W11.SetSize(total_entries);
  weights.W01.SetSize(total_entries);
  weights.W10.SetSize(total_entries);
  weights.W00 = 0.0;
  weights.W11 = 0.0;
  weights.W01 = 0.0;
  weights.W10 = 0.0;

  const IntegrationRule &gl = IntRules.Get(Geometry::SEGMENT, 2 * p + 1);
  const int n_gl = gl.GetNPoints();

  const int *sk_begin = span_knot_idx.GetData();
  const int *sk_end = sk_begin + num_spans;

  Vector N(p + 1), dN(p + 1), Nm1(std::max(1, p)), left(p + 1), right(p + 1);
  Vector r00, r11, r01, r10;
  Vector M_B, M_D, rhs, S, work;

  for (int i = 0; i < num_dofs; i++) {
    const int j0 = rule.j_begin[i];
    const int j1 = rule.j_end[i];
    const int q0 = rule.q_begin[i];
    const int q1 = rule.q_end[i];
    const int n_j = j1 - j0;
    const int n_q = q1 - q0;

    MFEM_VERIFY(n_j > 0 && n_q > 0,
                "Invalid trial/quadrature range for basis function in "
                "BuildWQWeights1D.");

    r00.SetSize(n_j);
    r11.SetSize(n_j);
    r01.SetSize(n_j);
    r10.SetSize(n_j);
    r00 = 0.0;
    r11 = 0.0;
    r01 = 0.0;
    r10 = 0.0;

    const int sb =
        static_cast<int>(std::lower_bound(sk_begin, sk_end, i) - sk_begin);
    const int se =
        static_cast<int>(std::upper_bound(sk_begin, sk_end, i + p) - sk_begin);
    MFEM_VERIFY(sb < se,
                "Basis function has empty support in BuildWQWeights1D.");

    // Exact 1D moments on supp(B_i) via (p + 1)-point Gauss-Legendre.
    for (int s = sb; s < se; s++) {
      const int k = span_knot_idx[s];
      const real_t u_a = kv[k];
      const real_t h = span_lengths[s];
      const int local_i = i - (k - p);

      for (int g = 0; g < n_gl; g++) {
        const IntegrationPoint &ip = gl.IntPoint(g);
        const real_t u = u_a + h * ip.x;
        const real_t w = h * ip.weight;

        EvalBSplineAndDeriv(kv, p, k, u, N, dN, Nm1, left, right);
        const real_t bi = N[local_i];
        const real_t di = dN[local_i];

        for (int r = 0; r <= p; r++) {
          const int j_loc = (k - p + r) - j0;
          r00[j_loc] += w * bi * N[r];
          r11[j_loc] += w * di * dN[r];
          r01[j_loc] += w * bi * dN[r];
          r10[j_loc] += w * di * N[r];
        }
      }
    }

    // Assemble column-major n_j x n_q system matrices M_B and M_D.
    M_B.SetSize(n_j * n_q);
    M_D.SetSize(n_j * n_q);
    M_B = 0.0;
    M_D = 0.0;

    for (int j_loc = 0; j_loc < n_j; j_loc++) {
      const int j = j0 + j_loc;
      const int qa = std::max(q0, rule.q_begin[j]);
      const int qb = std::min(q1, rule.q_end[j]);
      const int j_off = rule.trial_offsets[j];
      const int j_qb = rule.q_begin[j];

      for (int q = qa; q < qb; q++) {
        const int q_loc = q - q0;
        const int src = j_off + (q - j_qb);
        const int idx = j_loc + q_loc * n_j;
        M_B[idx] = rule.B_vals[src];
        M_D[idx] = rule.D_vals[src];
      }
    }

    const int ldb = std::max(n_j, n_q);
    const int row_off = weights.offsets[i];
    rhs.SetSize(ldb * 2);

    // System B: solve M_B * [W00, W10] = [r00, r10].
    rhs = 0.0;
    for (int j_loc = 0; j_loc < n_j; j_loc++) {
      rhs[j_loc] = r00[j_loc];
      rhs[ldb + j_loc] = r10[j_loc];
    }
    SolveMinNormLS(n_j, n_q, 2, M_B, rhs, S, work);
    for (int q_loc = 0; q_loc < n_q; q_loc++) {
      weights.W00[row_off + q_loc] = rhs[q_loc];
      weights.W10[row_off + q_loc] = rhs[ldb + q_loc];
    }

    // System D: solve M_D * [W11, W01] = [r11, r01] (rank-deficient by 1).
    if (p > 0) {
      rhs = 0.0;
      for (int j_loc = 0; j_loc < n_j; j_loc++) {
        rhs[j_loc] = r11[j_loc];
        rhs[ldb + j_loc] = r01[j_loc];
      }
      SolveMinNormLS(n_j, n_q, 2, M_D, rhs, S, work);
      for (int q_loc = 0; q_loc < n_q; q_loc++) {
        weights.W11[row_off + q_loc] = rhs[q_loc];
        weights.W01[row_off + q_loc] = rhs[ldb + q_loc];
      }
    }
  }
#endif
}

} // namespace mfem
