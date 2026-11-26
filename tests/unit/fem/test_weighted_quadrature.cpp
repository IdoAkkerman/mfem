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

#include "fem/weighted_quadrature.hpp"
#include "mfem.hpp"
#include "unit_tests.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;

namespace {

std::string FindMeshFile(const std::string &rel_from_root,
                         const std::string &rel_from_data = "") {
  std::vector<std::string> candidates;
  if (!rel_from_data.empty()) {
    candidates.push_back(::mfem_data_dir + "/" + rel_from_data);
  }
  candidates.push_back(::mfem_data_dir + "/../" + rel_from_root);
  candidates.push_back("../../" + rel_from_root);
  candidates.push_back("../../../" + rel_from_root);
  candidates.push_back(rel_from_root);

  for (const std::string &path : candidates) {
    std::ifstream ifs(path);
    if (ifs.good()) {
      return path;
    }
  }
  return "";
}

std::unique_ptr<KnotVector> MakeOpenKnotVector(int order, int num_spans,
                                               int interior_mult) {
  MFEM_VERIFY(order >= 1 && num_spans >= 1 && interior_mult >= 1 &&
                  interior_mult <= order,
              "Invalid parameters in MakeOpenKnotVector.");

  const int n_internal = (num_spans - 1) * interior_mult;
  const int ncp = order + 1 + n_internal;
  auto kv = std::make_unique<KnotVector>(order, ncp);

  int idx = 0;
  for (int k = 0; k <= order; k++) {
    (*kv)[idx++] = real_t(0.0);
  }
  for (int s = 1; s < num_spans; s++) {
    const real_t u = static_cast<real_t>(s) / static_cast<real_t>(num_spans);
    for (int m = 0; m < interior_mult; m++) {
      (*kv)[idx++] = u;
    }
  }
  for (int k = 0; k <= order; k++) {
    (*kv)[idx++] = real_t(1.0);
  }
  kv->GetElements();
  return kv;
}

void RefineNURBSTestMesh(Mesh &mesh, int order, int ref_levels, bool p_refine) {
  if (p_refine) {
    for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
    }
    mesh.DegreeElevate(order, order);
  } else {
    mesh.DegreeElevate(order, order);
    for (int l = 0; l < ref_levels; l++) {
      mesh.UniformRefinement();
    }
  }
}

#ifdef MFEM_USE_LAPACK

real_t SparseMaxAbs(const SparseMatrix &A) {
  REQUIRE(A.Finalized());
  const int nnz = A.NumNonZeroElems();
  const real_t *data = A.GetData();
  real_t max_val = 0.0;
  for (int k = 0; k < nnz; k++) {
    max_val = std::max(max_val, std::abs(data[k]));
  }
  return max_val;
}

real_t SparseFrobNorm(const SparseMatrix &A) {
  REQUIRE(A.Finalized());
  const int nnz = A.NumNonZeroElems();
  const real_t *data = A.GetData();
  real_t sum_sq = 0.0;
  for (int k = 0; k < nnz; k++) {
    sum_sq += data[k] * data[k];
  }
  return std::sqrt(sum_sq);
}

void ComputeSparseDiff(const SparseMatrix &A, const SparseMatrix &B,
                       real_t &max_diff, real_t &frob_diff) {
  REQUIRE(A.Height() == B.Height());
  REQUIRE(A.Width() == B.Width());
  REQUIRE(A.Finalized());
  REQUIRE(B.Finalized());

  std::unique_ptr<SparseMatrix> diff(Add(real_t(1.0), A, real_t(-1.0), B));
  max_diff = SparseMaxAbs(*diff);
  frob_diff = SparseFrobNorm(*diff);
}

real_t SparseSymmetryError(const SparseMatrix &A) {
  REQUIRE(A.Height() == A.Width());
  REQUIRE(A.Finalized());

  std::unique_ptr<SparseMatrix> AT(Transpose(A));
  std::unique_ptr<SparseMatrix> skew(Add(real_t(1.0), A, real_t(-1.0), *AT));
  return SparseMaxAbs(*skew);
}

void SolveWithCG(const SparseMatrix &A, const Vector &B, Vector &X,
                 real_t rtol = real_t(1e-12)) {
  DSmoother prec(A);
  CGSolver cg;
  cg.SetRelTol(rtol);
  cg.SetAbsTol(real_t(0.0));
  cg.SetMaxIter(20000);
  cg.SetPrintLevel(0);
  cg.SetPreconditioner(prec);
  cg.SetOperator(A);

  X.SetSize(A.Width());
  X = real_t(0.0);
  cg.Mult(B, X);
  REQUIRE(cg.GetConverged());
}

void BuildPatchRulesAndWeights(FiniteElementSpace &fespace, int patch,
                               const WQOptions &opts,
                               std::vector<WQRule1D> &rules,
                               std::vector<WQWeights1D> &weights) {
  NURBSExtension *next = fespace.GetNURBSext();
  REQUIRE(next != nullptr);

  Array<KnotVector *> kvs;
  next->GetPatchKnotVectors(patch, kvs);
  const int dim = kvs.Size();
  REQUIRE(dim == fespace.GetMesh()->Dimension());

  rules.resize(static_cast<std::size_t>(dim));
  weights.resize(static_cast<std::size_t>(dim));
  for (int d = 0; d < dim; d++) {
    BuildWQRule1D(*kvs[d], opts, rules[static_cast<std::size_t>(d)]);
    BuildWQWeights1D(*kvs[d], rules[static_cast<std::size_t>(d)],
                     weights[static_cast<std::size_t>(d)]);
  }
}

void ComputeExact1DMatrices(const KnotVector &kv, DenseMatrix &M00,
                            DenseMatrix &M11, DenseMatrix &M01,
                            DenseMatrix &M10) {
  const int p = kv.GetOrder();
  const int n = kv.GetNCP();
  M00.SetSize(n, n);
  M11.SetSize(n, n);
  M01.SetSize(n, n);
  M10.SetSize(n, n);
  M00 = real_t(0.0);
  M11 = real_t(0.0);
  M01 = real_t(0.0);
  M10 = real_t(0.0);

  const IntegrationRule &gl = IntRules.Get(Geometry::SEGMENT, 2 * p + 3);
  Vector shape(p + 1), dshape(p + 1);

  for (int k = p; k < n; k++) {
    const real_t h = kv[k + 1] - kv[k];
    if (h <= real_t(1e-14)) {
      continue;
    }

    for (int g = 0; g < gl.GetNPoints(); g++) {
      const IntegrationPoint &ip = gl.IntPoint(g);
      kv.CalcShape(shape, k, ip.x);
      kv.CalcDShape(dshape, k, ip.x);

      // Convert dB/dxi on [0, 1] to global parametric derivative dB/du.
      dshape /= h;
      const real_t w = h * ip.weight;

      for (int r = 0; r <= p; r++) {
        const int i = k - p + r;
        for (int c = 0; c <= p; c++) {
          const int j = k - p + c;
          M00(i, j) += w * shape[r] * shape[c];
          M11(i, j) += w * dshape[r] * dshape[c];
          M01(i, j) += w * shape[r] * dshape[c];
          M10(i, j) += w * dshape[r] * shape[c];
        }
      }
    }
  }
}

void CheckWQ1DExactness(const KnotVector &kv, const WQOptions &opts) {
  WQRule1D rule;
  WQWeights1D weights;
  BuildWQRule1D(kv, opts, rule);
  BuildWQWeights1D(kv, rule, weights);

  const int n = kv.GetNCP();
  REQUIRE(rule.order == kv.GetOrder());
  REQUIRE(rule.num_dofs == n);
  REQUIRE(weights.offsets.Size() == n + 1);

  DenseMatrix M00, M11, M01, M10;
  ComputeExact1DMatrices(kv, M00, M11, M01, M10);

  const real_t norm00 = M00.FNorm();
  const real_t norm11 = M11.FNorm();
  const real_t norm01 = M01.FNorm();
  const real_t norm10 = M10.FNorm();
  REQUIRE(norm00 > real_t(0.0));
  REQUIRE(norm11 > real_t(0.0));
  REQUIRE(norm01 > real_t(0.0));
  REQUIRE(norm10 > real_t(0.0));

  real_t max_err00 = 0.0;
  real_t max_err11 = 0.0;
  real_t max_err01 = 0.0;
  real_t max_err10 = 0.0;

  for (int i = 0; i < n; i++) {
    const int j0 = rule.j_begin[i];
    const int j1 = rule.j_end[i];
    const int qi0 = rule.q_begin[i];
    const int qi1 = rule.q_end[i];
    const int off_i = weights.offsets[i];

    for (int j = j0; j < j1; j++) {
      const int qj0 = rule.q_begin[j];
      const int qj1 = rule.q_end[j];
      const int off_j = rule.trial_offsets[j];

      const int qa = std::max(qi0, qj0);
      const int qb = std::min(qi1, qj1);

      real_t q00 = 0.0;
      real_t q11 = 0.0;
      real_t q01 = 0.0;
      real_t q10 = 0.0;
      for (int q = qa; q < qb; q++) {
        const int wi = off_i + (q - qi0);
        const int tj = off_j + (q - qj0);
        const real_t bj = rule.B_vals[tj];
        const real_t dj = rule.D_vals[tj];
        q00 += weights.W00[wi] * bj;
        q11 += weights.W11[wi] * dj;
        q01 += weights.W01[wi] * dj;
        q10 += weights.W10[wi] * bj;
      }

      max_err00 = std::max(max_err00, std::abs(q00 - M00(i, j)));
      max_err11 = std::max(max_err11, std::abs(q11 - M11(i, j)));
      max_err01 = std::max(max_err01, std::abs(q01 - M01(i, j)));
      max_err10 = std::max(max_err10, std::abs(q10 - M10(i, j)));
    }
  }

#ifdef MFEM_USE_SINGLE
  const real_t tol = real_t(1e-4);
#else
  const real_t tol = real_t(1e-11);
#endif
  CAPTURE(max_err00, norm00, max_err11, norm11, max_err01, norm01, max_err10,
          norm10);
  REQUIRE(max_err00 < tol * norm00);
  REQUIRE(max_err11 < tol * norm11);
  REQUIRE(max_err01 < tol * norm01);
  REQUIRE(max_err10 < tol * norm10);
}

real_t SolveAndCompareL2Error(Mesh &mesh, int order, Coefficient &coeff,
                              const WQOptions &opts, real_t &sym_ratio,
                              bool &is_diag_out) {
  NURBSFECollection fec(order, mesh.Dimension());
  FiniteElementSpace fespace(&mesh, &fec);

  Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  std::vector<WQRule1D> rules;
  std::vector<WQWeights1D> weights;
  BuildPatchRulesAndWeights(fespace, 0, opts, rules, weights);

  WQPatchMetric metric;
  EvaluatePatchMetric(mesh, 0, rules, coeff, metric);
  is_diag_out = metric.is_diagonal;

  ConstantCoefficient one(real_t(1.0));
  LinearForm b(&fespace);
  b.AddDomainIntegrator(new DomainLFIntegrator(one));
  b.Assemble();

  BilinearForm a_leg(&fespace);
  a_leg.SetAssemblyLevel(AssemblyLevel::LEGACY);
  a_leg.AddDomainIntegrator(new DiffusionIntegrator(coeff));
  a_leg.Assemble();
  a_leg.Finalize();

  GridFunction u_leg(&fespace);
  u_leg = real_t(0.0);
  Vector b_leg(b);
  OperatorHandle A_leg_op;
  Vector X_leg, B_leg;
  a_leg.FormLinearSystem(ess_tdof_list, u_leg, b_leg, A_leg_op, X_leg, B_leg);
  const SparseMatrix *A_leg = A_leg_op.As<SparseMatrix>();
  REQUIRE(A_leg != nullptr);

  GridFunction u_wq(&fespace);
  u_wq = real_t(0.0);
  Vector b_wq(b);
  SparseMatrix A_wq;
  Vector X_wq, B_wq;
  WQTimings timings;
  AssembleDiffusionWQ(fespace, coeff, ess_tdof_list, u_wq, b_wq, A_wq, X_wq,
                      B_wq, opts, &timings);

  const real_t max_A = SparseMaxAbs(A_wq);
  REQUIRE(max_A > real_t(0.0));
  sym_ratio = SparseSymmetryError(A_wq) / max_A;

  SolveWithCG(*A_leg, B_leg, X_leg, real_t(1e-12));
  SolveWithCG(A_wq, B_wq, X_wq, real_t(1e-12));

  a_leg.RecoverFEMSolution(X_leg, b_leg, u_leg);
  u_wq.SetFromTrueDofs(X_wq);

  GridFunction diff_u(&fespace);
  subtract(u_wq, u_leg, diff_u);

  ConstantCoefficient zero(real_t(0.0));
  const real_t l2_err = diff_u.ComputeL2Error(zero);
  const real_t l2_ref = u_leg.ComputeL2Error(zero);
  REQUIRE(l2_ref > real_t(0.0));
  return l2_err / l2_ref;
}

real_t SpatiallyVaryingCoeff(const Vector &x) {
  const real_t pi = static_cast<real_t>(M_PI);
  return real_t(1.0) + real_t(0.5) * std::sin(pi * x[0]) * std::cos(pi * x[1]) +
         real_t(0.25) * x[0] * x[0] * x[1];
}

#endif // MFEM_USE_LAPACK

} // namespace

TEST_CASE("WQ 1D rule structure and partition of unity",
          "[FEM][WeightedQuadrature]") {
  const int order = GENERATE(1, 2, 4);
  const int ref_levels = GENERATE(0, 2);
  const bool p_refine = GENERATE(false, true);
  CAPTURE(order, ref_levels, p_refine);

  const int num_spans = 1 << ref_levels;
  const int mult = p_refine ? order : 1;
  auto kv = MakeOpenKnotVector(order, num_spans, mult);

  WQOptions opts;
  WQRule1D rule;
  BuildWQRule1D(*kv, opts, rule);

  REQUIRE(rule.order == order);
  REQUIRE(rule.num_dofs == kv->GetNCP());
  REQUIRE(rule.num_spans == num_spans);
  REQUIRE(rule.num_qpts > 0);

  Vector sum_B(rule.num_qpts);
  Vector sum_D(rule.num_qpts);
  sum_B = real_t(0.0);
  sum_D = real_t(0.0);

  for (int j = 0; j < rule.num_dofs; j++) {
    const int q0 = rule.q_begin[j];
    const int q1 = rule.q_end[j];
    const int off = rule.trial_offsets[j];
    for (int q = q0; q < q1; q++) {
      sum_B[q] += rule.B_vals[off + (q - q0)];
      sum_D[q] += rule.D_vals[off + (q - q0)];
    }
  }

  const real_t tol = real_t(100.0) * std::numeric_limits<real_t>::epsilon();
  for (int q = 0; q < rule.num_qpts; q++) {
    CAPTURE(q, sum_B[q], sum_D[q]);
    REQUIRE(sum_B[q] == MFEM_Approx(real_t(1.0), tol, tol));
    REQUIRE(sum_D[q] == MFEM_Approx(real_t(0.0), tol, tol));
  }
}

#ifdef MFEM_USE_LAPACK

TEST_CASE("WQ 1D exactness of weight tables", "[FEM][WeightedQuadrature]") {
  SECTION("k-refinement and p-refinement sweep") {
    const int order = GENERATE(1, 2, 3, 5, 7);
    const int ref_levels = GENERATE(0, 2, 3);
    const bool p_refine = GENERATE(false, true);
    CAPTURE(order, ref_levels, p_refine);

    const int num_spans = 1 << ref_levels;
    const int mult = p_refine ? order : 1;
    auto kv = MakeOpenKnotVector(order, num_spans, mult);

    WQOptions opts;
    CheckWQ1DExactness(*kv, opts);
  }

  SECTION("Intermediate interior knot multiplicity") {
    const int order = GENERATE(3, 5, 7);
    const int mult = GENERATE(2, 3);
    if (mult < order) {
      CAPTURE(order, mult);
      auto kv = MakeOpenKnotVector(order, 4, mult);
      WQOptions opts;
      CheckWQ1DExactness(*kv, opts);
    }
  }
}

TEST_CASE("WQ 2D affine NURBS matches LEGACY", "[FEM][WeightedQuadrature]") {
  const std::string mesh_path =
      FindMeshFile("data/square-nurbs.mesh", "square-nurbs.mesh");
  REQUIRE(!mesh_path.empty());

  const int order = GENERATE(2, 4, 6);
  const int ref_levels = GENERATE(1, 2);
  const bool p_refine = GENERATE(false, true);
  CAPTURE(order, ref_levels, p_refine);

  Mesh mesh(mesh_path.c_str(), 1, 1);
  RefineNURBSTestMesh(mesh, order, ref_levels, p_refine);

  NURBSFECollection fec(order, mesh.Dimension());
  FiniteElementSpace fespace(&mesh, &fec);

  ConstantCoefficient one(real_t(1.0));
  WQOptions opts;
  opts.symmetrize = true;

  // Verify diagonal metric detection and uneliminated patch assembly.
  std::vector<WQRule1D> rules;
  std::vector<WQWeights1D> weights;
  BuildPatchRulesAndWeights(fespace, 0, opts, rules, weights);

  WQPatchMetric metric;
  EvaluatePatchMetric(mesh, 0, rules, one, metric);
  REQUIRE(metric.dim == 2);
  REQUIRE(metric.is_diagonal);

  SparseMatrix A_wq_raw;
  AssemblePatchDiffusionWQ(fespace, 0, rules, weights, metric, opts, A_wq_raw);
  if (!A_wq_raw.Finalized()) {
    A_wq_raw.Finalize();
  }

  BilinearForm a_leg(&fespace);
  a_leg.SetAssemblyLevel(AssemblyLevel::LEGACY);
  a_leg.AddDomainIntegrator(new DiffusionIntegrator(one));
  a_leg.Assemble();
  a_leg.Finalize();

  const SparseMatrix &A_leg_raw = a_leg.SpMat();
  const real_t max_raw = SparseMaxAbs(A_leg_raw);
  const real_t frob_raw = SparseFrobNorm(A_leg_raw);
  REQUIRE(max_raw > real_t(0.0));

  real_t max_diff_raw = 0.0;
  real_t frob_diff_raw = 0.0;
  ComputeSparseDiff(A_wq_raw, A_leg_raw, max_diff_raw, frob_diff_raw);

#ifdef MFEM_USE_SINGLE
  const real_t tol = real_t(1e-4);
#else
  const real_t tol = real_t(1e-11);
#endif
  CAPTURE(max_diff_raw, max_raw, frob_diff_raw, frob_raw);
  REQUIRE(max_diff_raw < tol * max_raw);
  REQUIRE(frob_diff_raw < tol * frob_raw);

  // Verify full system assembly after Dirichlet BC elimination.
  Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  LinearForm b(&fespace);
  b.AddDomainIntegrator(new DomainLFIntegrator(one));
  b.Assemble();

  GridFunction x_leg(&fespace);
  x_leg = real_t(0.0);
  Vector b_leg(b);
  OperatorHandle A_leg_op;
  Vector X_leg, B_leg;
  a_leg.FormLinearSystem(ess_tdof_list, x_leg, b_leg, A_leg_op, X_leg, B_leg);
  const SparseMatrix *A_leg = A_leg_op.As<SparseMatrix>();
  REQUIRE(A_leg != nullptr);

  GridFunction x_wq(&fespace);
  x_wq = real_t(0.0);
  Vector b_wq(b);
  SparseMatrix A_wq;
  Vector X_wq, B_wq;
  WQTimings timings;
  AssembleDiffusionWQ(fespace, one, ess_tdof_list, x_wq, b_wq, A_wq, X_wq, B_wq,
                      opts, &timings);

  const real_t max_elim = SparseMaxAbs(*A_leg);
  const real_t frob_elim = SparseFrobNorm(*A_leg);
  real_t max_diff_elim = 0.0;
  real_t frob_diff_elim = 0.0;
  ComputeSparseDiff(A_wq, *A_leg, max_diff_elim, frob_diff_elim);

  CAPTURE(max_diff_elim, max_elim, frob_diff_elim, frob_elim);
  REQUIRE(max_diff_elim < tol * max_elim);
  REQUIRE(frob_diff_elim < tol * frob_elim);

  Vector B_diff(B_wq);
  B_diff -= B_leg;
  REQUIRE(B_diff.Norml2() < tol * std::max(real_t(1.0), B_leg.Norml2()));

  REQUIRE(timings.rules_1d >= 0.0);
  REQUIRE(timings.weights_1d >= 0.0);
  REQUIRE(timings.metric_eval >= 0.0);
  REQUIRE(timings.patch_assemble >= 0.0);
  REQUIRE(timings.symmetrize >= 0.0);
  REQUIRE(timings.bc_elim >= 0.0);
}

TEST_CASE("WQ 2D distorted non-affine NURBS consistency",
          "[FEM][WeightedQuadrature]") {
  const std::string mesh_path =
      FindMeshFile("miniapps/nurbs/meshes/square-nurbs-distorted.mesh");
  REQUIRE(!mesh_path.empty());

  ConstantCoefficient one(real_t(1.0));
  WQOptions opts;
  opts.symmetrize = true;

#ifdef MFEM_USE_SINGLE
  const real_t sym_tol = real_t(1e-5);
#else
  const real_t sym_tol = real_t(1e-13);
#endif

  const int orders[3] = {2, 4, 5};
  const int refs[2] = {2, 3};
  real_t rel_l2[3][2] = {{0.0}};

  for (int oi = 0; oi < 3; oi++) {
    const int order = orders[oi];
    for (int ri = 0; ri < 2; ri++) {
      const int ref_levels = refs[ri];
      CAPTURE(order, ref_levels);

      Mesh mesh(mesh_path.c_str(), 1, 1);
      RefineNURBSTestMesh(mesh, order, ref_levels, false);

      real_t sym_ratio = 0.0;
      bool is_diag = true;
      rel_l2[oi][ri] =
          SolveAndCompareL2Error(mesh, order, one, opts, sym_ratio, is_diag);

      CAPTURE(is_diag, sym_ratio, rel_l2[oi][ri]);
      REQUIRE_FALSE(is_diag);
      REQUIRE(sym_ratio <= sym_tol);
    }

    // Error must shrink under h-refinement for each order.
    CAPTURE(order, rel_l2[oi][0], rel_l2[oi][1]);
    REQUIRE(rel_l2[oi][1] < rel_l2[oi][0]);
  }

  // Error must shrink rapidly with p at ref = 3 and meet quantitative bounds.
  CAPTURE(rel_l2[0][1], rel_l2[1][1], rel_l2[2][1]);
  REQUIRE(rel_l2[1][1] < rel_l2[0][1]);
  REQUIRE(rel_l2[2][1] < rel_l2[1][1]);
  REQUIRE(rel_l2[0][1] < real_t(1e-4));
  REQUIRE(rel_l2[2][1] < real_t(1e-7));
}

TEST_CASE("WQ non-zero Dirichlet BC elimination", "[FEM][WeightedQuadrature]") {
  const std::string mesh_path =
      FindMeshFile("data/square-nurbs.mesh", "square-nurbs.mesh");
  REQUIRE(!mesh_path.empty());

  const int order = GENERATE(2, 4);
  const bool p_refine = GENERATE(false, true);
  const int ref_levels = 2;
  CAPTURE(order, ref_levels, p_refine);

  Mesh mesh(mesh_path.c_str(), 1, 1);
  RefineNURBSTestMesh(mesh, order, ref_levels, p_refine);

  NURBSFECollection fec(order, mesh.Dimension());
  FiniteElementSpace fespace(&mesh, &fec);

  Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
  REQUIRE(ess_tdof_list.Size() > 0);

  ConstantCoefficient one(real_t(1.0));
  BilinearForm a_leg(&fespace);
  a_leg.SetAssemblyLevel(AssemblyLevel::LEGACY);
  a_leg.AddDomainIntegrator(new DiffusionIntegrator(one));
  a_leg.Assemble();
  a_leg.Finalize();

  GridFunction x_init(&fespace);
  x_init.Randomize(1);

  Vector b_init(fespace.GetVSize());
  b_init.Randomize(2);

  GridFunction x_leg(x_init);
  Vector b_leg(b_init);
  OperatorHandle A_leg_op;
  Vector X_leg, B_leg;
  a_leg.FormLinearSystem(ess_tdof_list, x_leg, b_leg, A_leg_op, X_leg, B_leg);
  const SparseMatrix *A_leg = A_leg_op.As<SparseMatrix>();
  REQUIRE(A_leg != nullptr);

  GridFunction x_wq(x_init);
  Vector b_wq(b_init);
  SparseMatrix A_wq;
  Vector X_wq, B_wq;
  WQOptions opts;
  opts.symmetrize = true;
  AssembleDiffusionWQ(fespace, one, ess_tdof_list, x_wq, b_wq, A_wq, X_wq, B_wq,
                      opts, nullptr);

#ifdef MFEM_USE_SINGLE
  const real_t tol = real_t(1e-4);
#else
  const real_t tol = real_t(1e-11);
#endif

  real_t max_diff = 0.0;
  real_t frob_diff = 0.0;
  ComputeSparseDiff(A_wq, *A_leg, max_diff, frob_diff);
  const real_t max_ref = SparseMaxAbs(*A_leg);
  REQUIRE(max_diff < tol * max_ref);

  Vector B_diff(B_wq);
  B_diff -= B_leg;
  REQUIRE(B_leg.Norml2() > real_t(0.0));
  REQUIRE(B_diff.Norml2() < tol * B_leg.Norml2());

  Vector X_diff(X_wq);
  X_diff -= X_leg;
  REQUIRE(X_leg.Norml2() > real_t(0.0));
  REQUIRE(X_diff.Norml2() < tol * X_leg.Norml2());

  // Verify identical linear system solutions with non-zero Dirichlet data.
  SolveWithCG(*A_leg, B_leg, X_leg, real_t(1e-12));
  SolveWithCG(A_wq, B_wq, X_wq, real_t(1e-12));

  Vector sol_diff(X_wq);
  sol_diff -= X_leg;
  REQUIRE(sol_diff.Norml2() < real_t(1e-10) * X_leg.Norml2());
}

TEST_CASE("WQ 3D affine NURBS matches LEGACY", "[FEM][WeightedQuadrature]") {
  const std::string mesh_path =
      FindMeshFile("data/cube-nurbs.mesh", "cube-nurbs.mesh");
  REQUIRE(!mesh_path.empty());

  const int order = GENERATE(2, 3, 4);
  const int ref_levels = GENERATE(1, 2);
  CAPTURE(order, ref_levels);

  Mesh mesh(mesh_path.c_str(), 1, 1);
  RefineNURBSTestMesh(mesh, order, ref_levels, false);

  NURBSFECollection fec(order, mesh.Dimension());
  FiniteElementSpace fespace(&mesh, &fec);

  ConstantCoefficient one(real_t(1.0));
  WQOptions opts;
  opts.symmetrize = true;

  std::vector<WQRule1D> rules;
  std::vector<WQWeights1D> weights;
  BuildPatchRulesAndWeights(fespace, 0, opts, rules, weights);

  WQPatchMetric metric;
  EvaluatePatchMetric(mesh, 0, rules, one, metric);
  REQUIRE(metric.dim == 3);
  REQUIRE(metric.is_diagonal);

  Array<int> ess_bdr(mesh.bdr_attributes.Max());
  ess_bdr = 1;
  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  LinearForm b(&fespace);
  b.AddDomainIntegrator(new DomainLFIntegrator(one));
  b.Assemble();

  BilinearForm a_leg(&fespace);
  a_leg.SetAssemblyLevel(AssemblyLevel::LEGACY);
  a_leg.AddDomainIntegrator(new DiffusionIntegrator(one));
  a_leg.Assemble();
  a_leg.Finalize();

  GridFunction x_leg(&fespace);
  x_leg = real_t(0.0);
  Vector b_leg(b);
  OperatorHandle A_leg_op;
  Vector X_leg, B_leg;
  a_leg.FormLinearSystem(ess_tdof_list, x_leg, b_leg, A_leg_op, X_leg, B_leg);
  const SparseMatrix *A_leg = A_leg_op.As<SparseMatrix>();
  REQUIRE(A_leg != nullptr);

  GridFunction x_wq(&fespace);
  x_wq = real_t(0.0);
  Vector b_wq(b);
  SparseMatrix A_wq;
  Vector X_wq, B_wq;
  AssembleDiffusionWQ(fespace, one, ess_tdof_list, x_wq, b_wq, A_wq, X_wq, B_wq,
                      opts, nullptr);

#ifdef MFEM_USE_SINGLE
  const real_t mat_tol = real_t(1e-4);
  const real_t sol_tol = real_t(1e-4);
#else
  const real_t mat_tol = real_t(1e-11);
  const real_t sol_tol = real_t(1e-10);
#endif

  real_t max_diff = 0.0;
  real_t frob_diff = 0.0;
  ComputeSparseDiff(A_wq, *A_leg, max_diff, frob_diff);
  const real_t max_ref = SparseMaxAbs(*A_leg);
  const real_t frob_ref = SparseFrobNorm(*A_leg);

  CAPTURE(max_diff, max_ref, frob_diff, frob_ref);
  REQUIRE(max_diff < mat_tol * max_ref);
  REQUIRE(frob_diff < mat_tol * frob_ref);

  SolveWithCG(*A_leg, B_leg, X_leg, real_t(1e-12));
  SolveWithCG(A_wq, B_wq, X_wq, real_t(1e-12));

  Vector sol_diff(X_wq);
  sol_diff -= X_leg;
  REQUIRE(X_leg.Norml2() > real_t(0.0));
  REQUIRE(sol_diff.Norml2() < sol_tol * X_leg.Norml2());
}

TEST_CASE("WQ non-unit coefficient and error handling",
          "[FEM][WeightedQuadrature]") {
  const std::string mesh_path =
      FindMeshFile("data/square-nurbs.mesh", "square-nurbs.mesh");
  REQUIRE(!mesh_path.empty());

  SECTION("Spatially varying FunctionCoefficient converges to LEGACY") {
    FunctionCoefficient var_coeff(SpatiallyVaryingCoeff);
    const bool sym = GENERATE(true, false);
    WQOptions opts;
    opts.symmetrize = sym;

    const int orders[2] = {2, 4};
    const int refs[2] = {2, 3};
    real_t rel_l2[2][2] = {{0.0}};

    for (int oi = 0; oi < 2; oi++) {
      const int order = orders[oi];
      for (int ri = 0; ri < 2; ri++) {
        const int ref_levels = refs[ri];
        CAPTURE(sym, order, ref_levels);

        Mesh mesh(mesh_path.c_str(), 1, 1);
        RefineNURBSTestMesh(mesh, order, ref_levels, false);

        real_t sym_ratio = 0.0;
        bool is_diag = false;
        rel_l2[oi][ri] = SolveAndCompareL2Error(mesh, order, var_coeff, opts,
                                                sym_ratio, is_diag);

        CAPTURE(is_diag, sym_ratio, rel_l2[oi][ri]);
        REQUIRE(is_diag);
        if (sym) {
          REQUIRE(sym_ratio <= real_t(1e-13));
        }
      }
      REQUIRE(rel_l2[oi][1] < rel_l2[oi][0]);
    }

    REQUIRE(rel_l2[1][1] < rel_l2[0][1]);
    REQUIRE(rel_l2[0][1] < real_t(1e-4));
    REQUIRE(rel_l2[1][1] < real_t(1e-7));
  }

  SECTION("Non-unit constant scaling and drop_tol") {
    Mesh mesh(mesh_path.c_str(), 1, 1);
    RefineNURBSTestMesh(mesh, 3, 2, false);

    NURBSFECollection fec(3, mesh.Dimension());
    FiniteElementSpace fespace(&mesh, &fec);

    ConstantCoefficient alpha(real_t(3.5));
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    Array<int> ess_tdof_list;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    BilinearForm a_leg(&fespace);
    a_leg.SetAssemblyLevel(AssemblyLevel::LEGACY);
    a_leg.AddDomainIntegrator(new DiffusionIntegrator(alpha));
    a_leg.Assemble();
    a_leg.Finalize();

    Vector b_init(fespace.GetVSize());
    b_init = real_t(1.0);
    GridFunction x_leg(&fespace), x_wq(&fespace);
    x_leg = real_t(0.0);
    x_wq = real_t(0.0);
    Vector b_leg(b_init), b_wq(b_init);

    OperatorHandle A_leg_op;
    Vector X_leg, B_leg;
    a_leg.FormLinearSystem(ess_tdof_list, x_leg, b_leg, A_leg_op, X_leg, B_leg);
    const SparseMatrix *A_leg = A_leg_op.As<SparseMatrix>();
    REQUIRE(A_leg != nullptr);

    WQOptions opts;
    opts.symmetrize = true;
    opts.drop_tol = real_t(1e-14);

    SparseMatrix A_wq;
    Vector X_wq, B_wq;
    AssembleDiffusionWQ(fespace, alpha, ess_tdof_list, x_wq, b_wq, A_wq, X_wq,
                        B_wq, opts, nullptr);

    real_t max_diff = 0.0;
    real_t frob_diff = 0.0;
    ComputeSparseDiff(A_wq, *A_leg, max_diff, frob_diff);
    REQUIRE(max_diff < real_t(1e-11) * SparseMaxAbs(*A_leg));
  }

#ifdef MFEM_USE_EXCEPTIONS
  SECTION("Invalid inputs throw ErrorException") {
    // Non-NURBS mesh rejected by AssembleDiffusionWQ.
    Mesh cart_mesh = Mesh::MakeCartesian2D(2, 2, Element::QUADRILATERAL);
    H1_FECollection h1_fec(2, 2);
    FiniteElementSpace h1_fes(&cart_mesh, &h1_fec);
    ConstantCoefficient one(real_t(1.0));
    Array<int> ess_tdofs;
    Vector x(h1_fes.GetVSize()), b(h1_fes.GetVSize()), X, B;
    x = real_t(0.0);
    b = real_t(1.0);
    SparseMatrix A;
    WQOptions opts;

    REQUIRE_THROWS_AS(AssembleDiffusionWQ(h1_fes, one, ess_tdofs, x, b, A, X, B,
                                          opts, nullptr),
                      ErrorException);

    // Unclamped (non-open) KnotVector rejected by BuildWQRule1D.
    KnotVector bad_kv(2, 4);
    bad_kv[0] = real_t(0.0);
    bad_kv[1] = real_t(0.1);
    bad_kv[2] = real_t(0.2);
    bad_kv[3] = real_t(0.5);
    bad_kv[4] = real_t(0.8);
    bad_kv[5] = real_t(0.9);
    bad_kv[6] = real_t(1.0);
    WQRule1D rule;
    REQUIRE_THROWS_AS(BuildWQRule1D(bad_kv, opts, rule), ErrorException);

    // Mismatched WQRule1D rejected by BuildWQWeights1D.
    auto kv1 = MakeOpenKnotVector(2, 2, 1);
    auto kv2 = MakeOpenKnotVector(3, 2, 1);
    BuildWQRule1D(*kv1, opts, rule);
    WQWeights1D weights;
    REQUIRE_THROWS_AS(BuildWQWeights1D(*kv2, rule, weights), ErrorException);
  }
#endif // MFEM_USE_EXCEPTIONS
}

#endif // MFEM_USE_LAPACK
