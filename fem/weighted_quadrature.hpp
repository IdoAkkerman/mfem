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

#ifndef MFEM_WEIGHTED_QUADRATURE
#define MFEM_WEIGHTED_QUADRATURE

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "../linalg/sparsemat.hpp"
#include "../linalg/vector.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/nurbs.hpp"
#include "coefficient.hpp"
#include "fespace.hpp"

#include <vector>

namespace mfem {

/// Configuration options for Weighted Quadrature (WQ) diffusion assembly.
struct WQOptions {
  /// Symmetrize assembled matrix via 0.5 * (A + A^T) for CG/MINRES.
  bool symmetrize = true;
  /// Points per non-empty interior span for maximal C^{p-1} continuity.
  int interior_pts_per_span = 2;
  /// Points on first/last non-empty spans; <= 0 selects p + 1 points.
  int boundary_pts_per_span = 0;
  /// Optional drop tolerance applied during assembly.
  real_t drop_tol = 0.0;
};

/// Wall-clock timings (seconds) for stages of WQ diffusion assembly.
struct WQTimings {
  double rules_1d = 0.0;
  double weights_1d = 0.0;
  double metric_eval = 0.0;
  double patch_assemble = 0.0;
  double symmetrize = 0.0;
  double bc_elim = 0.0;
};

/// 1D Weighted Quadrature points, basis support ranges, and trial tables.
struct WQRule1D {
  int order = 0;
  int num_dofs = 0;
  /// Number of non-empty knot spans.
  int num_spans = 0;
  int num_qpts = 0;

  /// Global 1D parametric quadrature points (size num_qpts).
  Vector xi;
  /// Non-empty span index in 0..num_spans-1 owning each qpt (size num_qpts).
  Array<int> span_of_qpt;
  /// Reference coordinate in [0, 1] inside span_of_qpt[q] (size num_qpts).
  Vector local_xi;
  /// Span lengths u_{s+1} - u_s for each non-empty span s (size num_spans).
  Vector span_lengths;
  /// Quadrature point range [span_q_begin[s], span_q_end[s]) per span.
  Array<int> span_q_begin;
  Array<int> span_q_end;

  /// Contiguous active qpt range [q_begin[i], q_end[i]) for basis i.
  Array<int> q_begin;
  Array<int> q_end;
  /// Contiguous overlapping trial basis range [j_begin[i], j_end[i]) for i.
  Array<int> j_begin;
  Array<int> j_end;

  /// Row offsets for compact trial tables (size num_dofs + 1).
  Array<int> trial_offsets;
  /// Trial values B_j(xi_q) at trial_offsets[j] + (q - q_begin[j]).
  Vector B_vals;
  /// Global parametric derivatives dB_j/du(xi_q) at the same offsets.
  Vector D_vals;
};

/// 1D Weighted Quadrature weight tables sharing WQRule1D::trial_offsets layout.
struct WQWeights1D {
  /// Row offsets (size num_dofs + 1), matching WQRule1D::trial_offsets.
  Array<int> offsets;
  /// Weights exact for \int B_i(u) B_j(u) du pairing with trial B_j.
  Vector W00;
  /// Weights exact for \int B_i'(u) B_j'(u) du pairing with trial B_j'.
  Vector W11;
  /// Weights exact for \int B_i(u) B_j'(u) du pairing with trial B_j'.
  Vector W01;
  /// Weights exact for \int B_i'(u) B_j(u) du pairing with trial B_j.
  Vector W10;
};

/// Packed symmetric diffusion pullback metric evaluated on patch WQ points.
struct WQPatchMetric {
  /// Spatial dimension (2 or 3).
  int dim = 0;
  /// Number of 1D WQ points per parametric direction (size dim).
  Array<int> q_sizes;
  /// True if all off-diagonal G_{ab} (a != b) vanish within tolerance.
  bool is_diagonal = false;
  /// False when all patch weights are equal within 1e-12 * max_w.
  bool is_rational = false;
  /// Packed symmetric metric entries (direction 0 fastest; stride 3 in 2D,
  /// stride 6 in 3D). Scaled by 1/W^2 when is_rational == true.
  Vector G;
  /// Convection vector v_a per WQ point when is_rational == true (stride dim,
  /// size dim * Q_total; empty when !is_rational).
  Vector V;
  /// Reaction scalar c per WQ point when is_rational == true (size Q_total;
  /// empty when !is_rational).
  Vector C;
  /// Patch-local control point weights ordered in patch tensor order
  /// (size \prod_c rules[c].num_dofs; populated when is_rational == true).
  Vector patch_weights;
};

/// Build continuity-aware 1D WQ quadrature points and banded trial tables.
void BuildWQRule1D(const KnotVector &kv, const WQOptions &opts, WQRule1D &rule);

/// Solve 1D WQ weight systems for test-trial pairings (W00, W11, W01, W10).
void BuildWQWeights1D(const KnotVector &kv, const WQRule1D &rule,
                      WQWeights1D &weights);

/// Evaluate the symmetric diffusion metric tensor on a NURBS patch's WQ grid.
void EvaluatePatchMetric(const Mesh &mesh, int patch,
                         const std::vector<WQRule1D> &rules, Coefficient &coeff,
                         WQPatchMetric &metric);

/// Assemble a single NURBS patch's diffusion contribution using 1D WQ rules.
void AssemblePatchDiffusionWQ(const FiniteElementSpace &fes, int patch,
                              const std::vector<WQRule1D> &rules,
                              const std::vector<WQWeights1D> &weights,
                              const WQPatchMetric &metric,
                              const WQOptions &opts, SparseMatrix &A);

/// Assemble global WQ diffusion linear system with essential BC elimination.
void AssembleDiffusionWQ(const FiniteElementSpace &fes, Coefficient &coeff,
                         const Array<int> &ess_tdof_list, Vector &x, Vector &b,
                         SparseMatrix &A, Vector &X, Vector &B,
                         const WQOptions &opts = WQOptions(),
                         WQTimings *timings = nullptr);

} // namespace mfem

#endif
