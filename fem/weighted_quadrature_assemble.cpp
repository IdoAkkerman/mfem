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
#include "eltrans.hpp"
#include "fe/fe_nurbs.hpp"
#include "gridfunc.hpp"
#include "weighted_quadrature.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

namespace mfem {

namespace {

/// Extract knot indices of non-empty knot spans from an open KnotVector.
void ExtractNonEmptySpanIndices(const KnotVector &kv,
                                Array<int> &span_knot_idx) {
  const int p = kv.GetOrder();
  const int n = kv.GetNCP();
  const int m = kv.Size();
  MFEM_VERIFY(m == n + p + 1 && m > 0,
              "Invalid open KnotVector in Weighted Quadrature assembly.");

  const real_t domain_len = kv[m - 1] - kv[0];
  MFEM_VERIFY(domain_len > 0.0, "KnotVector domain length must be positive.");
  const real_t knot_tol = static_cast<real_t>(1e-14) * domain_len;

  span_knot_idx.SetSize(0);
  for (int k = p; k < n; k++) {
    if (kv[k + 1] - kv[k] > knot_tol) {
      span_knot_idx.Append(k);
    }
  }
}

/// Retrieve element range and KnotVectors for a patch via public
/// NURBSExtension::LoadFE.
void GetPatchInfoFromNURBSExt(const NURBSExtension *nurbs_ext,
                              int total_mesh_ne, int dim, int patch,
                              int &patch_el_offset, int &patch_ne,
                              Array<const KnotVector *> &patch_kvs) {
  MFEM_VERIFY(nurbs_ext != nullptr, "Null NURBSExtension pointer.");
  MFEM_VERIFY(nurbs_ext->GetNE() == nurbs_ext->GetGNE(),
              "Weighted Quadrature requires a serial NURBS mesh with all "
              "elements active.");

  const int np = nurbs_ext->GetNP();
  MFEM_VERIFY(patch >= 0 && patch < np, "NURBS patch index out of range.");

  std::unique_ptr<NURBS2DFiniteElement> nfe2d;
  std::unique_ptr<NURBS3DFiniteElement> nfe3d;
  const FiniteElement *fe = nullptr;
  const NURBSFiniteElement *nfe = nullptr;
  if (dim == 2) {
    nfe2d.reset(new NURBS2DFiniteElement(1));
    fe = nfe2d.get();
    nfe = nfe2d.get();
  } else {
    nfe3d.reset(new NURBS3DFiniteElement(1));
    fe = nfe3d.get();
    nfe = nfe3d.get();
  }

  patch_el_offset = 0;
  patch_ne = 0;
  int el = 0;

  for (int q = 0; q < np; q++) {
    MFEM_VERIFY(el < total_mesh_ne,
                "Element index out of range while scanning NURBS patches.");
    nurbs_ext->LoadFE(el, fe);
    MFEM_VERIFY(nfe->GetPatch() == q,
                "Unexpected non-sequential patch element ordering in "
                "NURBSExtension.");

    const Array<const KnotVector *> &q_kv = nfe->KnotVectors();
    MFEM_VERIFY(q_kv.Size() == dim,
                "Unexpected patch KnotVector count in NURBSExtension.");

    int ne_q = 1;
    for (int d = 0; d < dim; d++) {
      ne_q *= q_kv[d]->GetNE();
    }
    if (q == patch) {
      patch_el_offset = el;
      patch_ne = ne_q;
      patch_kvs = q_kv;
    }
    el += ne_q;
  }

  MFEM_VERIFY(el == total_mesh_ne,
              "Sum of NURBS patch elements does not match mesh.GetNE().");
}

/// Retrieve element range and space KnotVectors for a patch via fes.GetFE().
void GetPatchInfoFromFESpace(const FiniteElementSpace &fes, int dim, int patch,
                             int &patch_el_offset, int &patch_ne,
                             Array<const KnotVector *> &patch_kvs) {
  const Mesh *mesh = fes.GetMesh();
  MFEM_VERIFY(mesh != nullptr && mesh->NURBSExt() != nullptr,
              "Weighted Quadrature requires a NURBS Mesh.");

  const NURBSExtension *mesh_ext = mesh->NURBSExt();
  MFEM_VERIFY(mesh_ext->GetNE() == mesh_ext->GetGNE(),
              "Weighted Quadrature requires a serial NURBS mesh with all "
              "elements active.");

  const int np = mesh_ext->GetNP();
  const int total_mesh_ne = mesh->GetNE();
  MFEM_VERIFY(patch >= 0 && patch < np, "NURBS patch index out of range.");

  patch_el_offset = 0;
  patch_ne = 0;
  int el = 0;

  for (int q = 0; q < np; q++) {
    MFEM_VERIFY(el < total_mesh_ne,
                "Element index out of range while scanning NURBS patches.");
    const FiniteElement *fe = fes.GetFE(el);
    const auto *nfe = dynamic_cast<const NURBSFiniteElement *>(fe);
    MFEM_VERIFY(nfe != nullptr,
                "Expected NURBSFiniteElement from FiniteElementSpace.");
    MFEM_VERIFY(nfe->GetPatch() == q,
                "Unexpected non-sequential patch element ordering in "
                "FiniteElementSpace.");

    const Array<const KnotVector *> &q_kv = nfe->KnotVectors();
    MFEM_VERIFY(q_kv.Size() == dim,
                "Unexpected patch KnotVector count in FiniteElementSpace.");

    int ne_q = 1;
    for (int d = 0; d < dim; d++) {
      ne_q *= q_kv[d]->GetNE();
    }
    if (q == patch) {
      patch_el_offset = el;
      patch_ne = ne_q;
      patch_kvs = q_kv;
    }
    el += ne_q;
  }

  MFEM_VERIFY(el == total_mesh_ne,
              "Sum of NURBS patch elements does not match mesh.GetNE().");
}

/// Extract patch control-point weights in patch tensor order (dir 0 fastest)
/// from mesh.NURBSExt()->GetWeights(), verify w > 0, and detect rationality.
void ExtractPatchWeights(const Mesh &mesh,
                         const Array<const KnotVector *> &patch_kv,
                         int patch_el_offset, bool &is_rational,
                         Vector &patch_weights) {
  is_rational = false;
  patch_weights.SetSize(0);

  const NURBSExtension *nurbs_ext = mesh.NURBSExt();
  const Vector &weights = nurbs_ext->GetWeights();
  if (weights.Size() == 0) {
    return;
  }

  const int dim = patch_kv.Size();
  std::vector<Array<int>> span_knot_idx(static_cast<std::size_t>(dim));
  int num_patch_cp = 1;
  for (int d = 0; d < dim; d++) {
    ExtractNonEmptySpanIndices(*patch_kv[d],
                               span_knot_idx[static_cast<std::size_t>(d)]);
    MFEM_VERIFY(span_knot_idx[static_cast<std::size_t>(d)].Size() ==
                    patch_kv[d]->GetNE(),
                "Non-empty span count mismatch in ExtractPatchWeights.");
    num_patch_cp *= patch_kv[d]->GetNCP();
  }

  Vector w_patch(num_patch_cp);
  real_t *wp = w_patch.HostWrite();
  std::vector<char> seen(static_cast<std::size_t>(num_patch_cp), 0);

  const real_t *w_data = weights.HostRead();
  const GridFunction *nodes = mesh.GetNodes();
  const FiniteElementSpace *geo_fes =
      (nodes != nullptr) ? nodes->FESpace() : nullptr;
  const bool use_geo_fes =
      (geo_fes != nullptr && geo_fes->GetNDofs() == weights.Size());

  std::unique_ptr<NURBS2DFiniteElement> nfe2d;
  std::unique_ptr<NURBS3DFiniteElement> nfe3d;
  const FiniteElement *fe = nullptr;
  const NURBSFiniteElement *nfe = nullptr;
  if (!use_geo_fes) {
    if (dim == 2) {
      nfe2d.reset(new NURBS2DFiniteElement(1));
      fe = nfe2d.get();
      nfe = nfe2d.get();
    } else {
      nfe3d.reset(new NURBS3DFiniteElement(1));
      fe = nfe3d.get();
      nfe = nfe3d.get();
    }
  }

  auto record_weight = [&](int p_idx, real_t w) {
    MFEM_VERIFY(w > static_cast<real_t>(0.0),
                "Non-positive NURBS weight detected on patch.");
    const std::size_t idx = static_cast<std::size_t>(p_idx);
    if (seen[idx]) {
      const real_t scale = std::max(wp[p_idx], w);
      MFEM_VERIFY(std::abs(wp[p_idx] - w) <= static_cast<real_t>(1e-14) * scale,
                  "Inconsistent NURBS control-point weight on patch.");
    } else {
      wp[p_idx] = w;
      seen[idx] = 1;
    }
  };

  Array<int> el_dofs;
  if (dim == 2) {
    const int p0 = patch_kv[0]->GetOrder();
    const int p1 = patch_kv[1]->GetOrder();
    const int n0 = patch_kv[0]->GetNCP();
    const int ns0 = patch_kv[0]->GetNE();
    const int ns1 = patch_kv[1]->GetNE();
    const int expected_nd = (p0 + 1) * (p1 + 1);

    for (int s1 = 0; s1 < ns1; s1++) {
      const int k1 = span_knot_idx[1][s1];
      for (int s0 = 0; s0 < ns0; s0++) {
        const int k0 = span_knot_idx[0][s0];
        const int el = patch_el_offset + s0 + ns0 * s1;
        const real_t *el_w = nullptr;

        if (use_geo_fes) {
          geo_fes->GetElementDofs(el, el_dofs);
          MFEM_VERIFY(el_dofs.Size() == expected_nd,
                      "Unexpected 2D geometry element DoF count.");
        } else {
          nurbs_ext->LoadFE(el, fe);
          MFEM_VERIFY(nfe->Weights().Size() == expected_nd,
                      "Unexpected 2D NURBS element weight count.");
          el_w = nfe->Weights().HostRead();
        }

        for (int r1 = 0; r1 <= p1; r1++) {
          const int i1 = k1 - p1 + r1;
          for (int r0 = 0; r0 <= p0; r0++) {
            const int i0 = k0 - p0 + r0;
            const int p_idx = i0 + n0 * i1;
            const int loc = r0 + (p0 + 1) * r1;
            real_t w = static_cast<real_t>(0.0);
            if (use_geo_fes) {
              const int dof = el_dofs[loc];
              MFEM_VERIFY(dof >= 0 && dof < weights.Size(),
                          "Geometry control-point index out of range.");
              w = w_data[dof];
            } else {
              w = el_w[loc];
            }
            record_weight(p_idx, w);
          }
        }
      }
    }
  } else {
    const int p0 = patch_kv[0]->GetOrder();
    const int p1 = patch_kv[1]->GetOrder();
    const int p2 = patch_kv[2]->GetOrder();
    const int n0 = patch_kv[0]->GetNCP();
    const int n1 = patch_kv[1]->GetNCP();
    const int ns0 = patch_kv[0]->GetNE();
    const int ns1 = patch_kv[1]->GetNE();
    const int ns2 = patch_kv[2]->GetNE();
    const int expected_nd = (p0 + 1) * (p1 + 1) * (p2 + 1);

    for (int s2 = 0; s2 < ns2; s2++) {
      const int k2 = span_knot_idx[2][s2];
      for (int s1 = 0; s1 < ns1; s1++) {
        const int k1 = span_knot_idx[1][s1];
        for (int s0 = 0; s0 < ns0; s0++) {
          const int k0 = span_knot_idx[0][s0];
          const int el = patch_el_offset + s0 + ns0 * (s1 + ns1 * s2);
          const real_t *el_w = nullptr;

          if (use_geo_fes) {
            geo_fes->GetElementDofs(el, el_dofs);
            MFEM_VERIFY(el_dofs.Size() == expected_nd,
                        "Unexpected 3D geometry element DoF count.");
          } else {
            nurbs_ext->LoadFE(el, fe);
            MFEM_VERIFY(nfe->Weights().Size() == expected_nd,
                        "Unexpected 3D NURBS element weight count.");
            el_w = nfe->Weights().HostRead();
          }

          for (int r2 = 0; r2 <= p2; r2++) {
            const int i2 = k2 - p2 + r2;
            for (int r1 = 0; r1 <= p1; r1++) {
              const int i1 = k1 - p1 + r1;
              for (int r0 = 0; r0 <= p0; r0++) {
                const int i0 = k0 - p0 + r0;
                const int p_idx = i0 + n0 * (i1 + n1 * i2);
                const int loc = r0 + (p0 + 1) * (r1 + (p1 + 1) * r2);
                real_t w = static_cast<real_t>(0.0);
                if (use_geo_fes) {
                  const int dof = el_dofs[loc];
                  MFEM_VERIFY(dof >= 0 && dof < weights.Size(),
                              "Geometry control-point index out of "
                              "range.");
                  w = w_data[dof];
                } else {
                  w = el_w[loc];
                }
                record_weight(p_idx, w);
              }
            }
          }
        }
      }
    }
  }

  real_t min_w = std::numeric_limits<real_t>::max();
  real_t max_w = static_cast<real_t>(0.0);
  for (int k = 0; k < num_patch_cp; k++) {
    MFEM_VERIFY(seen[static_cast<std::size_t>(k)],
                "Unmapped control-point weight on NURBS patch.");
    min_w = std::min(min_w, wp[k]);
    max_w = std::max(max_w, wp[k]);
  }

  is_rational = (max_w - min_w > static_cast<real_t>(1e-12) * max_w);
  if (is_rational) {
    patch_weights.Swap(w_patch);
  }
}

/// Evaluate NURBS denominator W(q) and global parametric derivatives dW/du_a(q)
/// at all tensor-product WQ points via directional banded sum factorization.
void EvaluateNURBSDenominator(const std::vector<WQRule1D> &rules,
                              const Vector &patch_weights,
                              std::vector<real_t> &W,
                              std::vector<std::vector<real_t>> &dW) {
  const int dim = static_cast<int>(rules.size());
  const real_t *pw = patch_weights.HostRead();

  const int n0 = rules[0].num_dofs;
  const int Q0 = rules[0].num_qpts;
  const int *qb0 = rules[0].q_begin.HostRead();
  const int *qe0 = rules[0].q_end.HostRead();
  const int *off0 = rules[0].trial_offsets.HostRead();
  const real_t *B0 = rules[0].B_vals.HostRead();
  const real_t *D0 = rules[0].D_vals.HostRead();

  const int n1 = rules[1].num_dofs;
  const int Q1 = rules[1].num_qpts;
  const int *qb1 = rules[1].q_begin.HostRead();
  const int *qe1 = rules[1].q_end.HostRead();
  const int *off1 = rules[1].trial_offsets.HostRead();
  const real_t *B1 = rules[1].B_vals.HostRead();
  const real_t *D1 = rules[1].D_vals.HostRead();

  dW.resize(static_cast<std::size_t>(dim));

  if (dim == 2) {
    const std::size_t Q0u = static_cast<std::size_t>(Q0);
    const std::size_t Q01 = Q0u * static_cast<std::size_t>(Q1);
    const std::size_t T_sz = Q0u * static_cast<std::size_t>(n1);

    std::vector<real_t> T_B0(T_sz, static_cast<real_t>(0.0));
    std::vector<real_t> T_D0(T_sz, static_cast<real_t>(0.0));

    for (int k1 = 0; k1 < n1; k1++) {
      const std::size_t t_off = Q0u * static_cast<std::size_t>(k1);
      const int w_off = n0 * k1;
      for (int k0 = 0; k0 < n0; k0++) {
        const real_t w = pw[w_off + k0];
        const int qb = qb0[k0];
        const int qe = qe0[k0];
        const int off = off0[k0];
        for (int q0 = qb; q0 < qe; q0++) {
          const int idx = off + (q0 - qb);
          const std::size_t dst = t_off + static_cast<std::size_t>(q0);
          T_B0[dst] += w * B0[idx];
          T_D0[dst] += w * D0[idx];
        }
      }
    }

    W.assign(Q01, static_cast<real_t>(0.0));
    dW[0].assign(Q01, static_cast<real_t>(0.0));
    dW[1].assign(Q01, static_cast<real_t>(0.0));

    for (int k1 = 0; k1 < n1; k1++) {
      const int qb = qb1[k1];
      const int qe = qe1[k1];
      const int off = off1[k1];
      const std::size_t t_off = Q0u * static_cast<std::size_t>(k1);
      const real_t *tb = T_B0.data() + t_off;
      const real_t *td = T_D0.data() + t_off;

      for (int q1 = qb; q1 < qe; q1++) {
        const int idx = off + (q1 - qb);
        const real_t b1 = B1[idx];
        const real_t d1 = D1[idx];
        const std::size_t out_off = Q0u * static_cast<std::size_t>(q1);
        for (std::size_t q0 = 0; q0 < Q0u; q0++) {
          W[out_off + q0] += tb[q0] * b1;
          dW[0][out_off + q0] += td[q0] * b1;
          dW[1][out_off + q0] += tb[q0] * d1;
        }
      }
    }
  } else {
    const int n2 = rules[2].num_dofs;
    const int Q2 = rules[2].num_qpts;
    const int *qb2 = rules[2].q_begin.HostRead();
    const int *qe2 = rules[2].q_end.HostRead();
    const int *off2 = rules[2].trial_offsets.HostRead();
    const real_t *B2 = rules[2].B_vals.HostRead();
    const real_t *D2 = rules[2].D_vals.HostRead();

    const std::size_t Q0u = static_cast<std::size_t>(Q0);
    const std::size_t Q1u = static_cast<std::size_t>(Q1);
    const std::size_t Q01 = Q0u * Q1u;
    const std::size_t Q_total = Q01 * static_cast<std::size_t>(Q2);
    const std::size_t T1_sz =
        Q0u * static_cast<std::size_t>(n1) * static_cast<std::size_t>(n2);
    const std::size_t T2_sz = Q01 * static_cast<std::size_t>(n2);

    // Stage 1: contract direction 0.
    std::vector<real_t> T1_B(T1_sz, static_cast<real_t>(0.0));
    std::vector<real_t> T1_D0(T1_sz, static_cast<real_t>(0.0));

    for (int k2 = 0; k2 < n2; k2++) {
      for (int k1 = 0; k1 < n1; k1++) {
        const int k12 = k1 + n1 * k2;
        const std::size_t t_off = Q0u * static_cast<std::size_t>(k12);
        const int w_off = n0 * k12;
        for (int k0 = 0; k0 < n0; k0++) {
          const real_t w = pw[w_off + k0];
          const int qb = qb0[k0];
          const int qe = qe0[k0];
          const int off = off0[k0];
          for (int q0 = qb; q0 < qe; q0++) {
            const int idx = off + (q0 - qb);
            const std::size_t dst = t_off + static_cast<std::size_t>(q0);
            T1_B[dst] += w * B0[idx];
            T1_D0[dst] += w * D0[idx];
          }
        }
      }
    }

    // Stage 2: contract direction 1.
    std::vector<real_t> T2_B(T2_sz, static_cast<real_t>(0.0));
    std::vector<real_t> T2_D0(T2_sz, static_cast<real_t>(0.0));
    std::vector<real_t> T2_D1(T2_sz, static_cast<real_t>(0.0));

    for (int k2 = 0; k2 < n2; k2++) {
      const std::size_t k2_q01 = Q01 * static_cast<std::size_t>(k2);
      for (int k1 = 0; k1 < n1; k1++) {
        const int qb = qb1[k1];
        const int qe = qe1[k1];
        const int off = off1[k1];
        const std::size_t t1_off = Q0u * static_cast<std::size_t>(k1 + n1 * k2);
        const real_t *t1_b = T1_B.data() + t1_off;
        const real_t *t1_d0 = T1_D0.data() + t1_off;

        for (int q1 = qb; q1 < qe; q1++) {
          const int idx = off + (q1 - qb);
          const real_t b1 = B1[idx];
          const real_t d1 = D1[idx];
          const std::size_t t2_off =
              k2_q01 + Q0u * static_cast<std::size_t>(q1);
          for (std::size_t q0 = 0; q0 < Q0u; q0++) {
            T2_B[t2_off + q0] += t1_b[q0] * b1;
            T2_D0[t2_off + q0] += t1_d0[q0] * b1;
            T2_D1[t2_off + q0] += t1_b[q0] * d1;
          }
        }
      }
    }

    T1_B.clear();
    T1_D0.clear();

    // Stage 3: contract direction 2.
    W.assign(Q_total, static_cast<real_t>(0.0));
    dW[0].assign(Q_total, static_cast<real_t>(0.0));
    dW[1].assign(Q_total, static_cast<real_t>(0.0));
    dW[2].assign(Q_total, static_cast<real_t>(0.0));

    for (int k2 = 0; k2 < n2; k2++) {
      const int qb = qb2[k2];
      const int qe = qe2[k2];
      const int off = off2[k2];
      const std::size_t t2_off = Q01 * static_cast<std::size_t>(k2);
      const real_t *t2_b = T2_B.data() + t2_off;
      const real_t *t2_d0 = T2_D0.data() + t2_off;
      const real_t *t2_d1 = T2_D1.data() + t2_off;

      for (int q2 = qb; q2 < qe; q2++) {
        const int idx = off + (q2 - qb);
        const real_t b2 = B2[idx];
        const real_t d2 = D2[idx];
        const std::size_t out_off = Q01 * static_cast<std::size_t>(q2);
        for (std::size_t idx01 = 0; idx01 < Q01; idx01++) {
          W[out_off + idx01] += t2_b[idx01] * b2;
          dW[0][out_off + idx01] += t2_d0[idx01] * b2;
          dW[1][out_off + idx01] += t2_d1[idx01] * b2;
          dW[2][out_off + idx01] += t2_b[idx01] * d2;
        }
      }
    }
  }
}

/// Build patch-local tensor-product DoF (direction 0 fastest) to global DoF
/// map using element DoF tables on the patch.
void BuildPatchDofMap(const FiniteElementSpace &fes,
                      const Array<const KnotVector *> &kvs,
                      const std::vector<WQRule1D> &rules, int patch_el_offset,
                      Array<int> &patch_dofs) {
  const int dim = static_cast<int>(rules.size());
  MFEM_VERIFY(kvs.Size() == dim,
              "Patch KnotVector dimension mismatch in BuildPatchDofMap.");

  std::vector<Array<int>> span_knot_idx(static_cast<std::size_t>(dim));
  int num_patch_dofs = 1;
  for (int d = 0; d < dim; d++) {
    ExtractNonEmptySpanIndices(*kvs[d],
                               span_knot_idx[static_cast<std::size_t>(d)]);
    MFEM_VERIFY(kvs[d]->GetOrder() == rules[d].order &&
                    kvs[d]->GetNCP() == rules[d].num_dofs &&
                    span_knot_idx[static_cast<std::size_t>(d)].Size() ==
                        rules[d].num_spans,
                "WQRule1D does not match patch KnotVector in "
                "AssemblePatchDiffusionWQ.");
    num_patch_dofs *= rules[d].num_dofs;
  }

  patch_dofs.SetSize(num_patch_dofs);
  patch_dofs = -1;

  Array<int> el_dofs;
  if (dim == 2) {
    const int p0 = rules[0].order;
    const int p1 = rules[1].order;
    const int n0 = rules[0].num_dofs;
    const int ns0 = rules[0].num_spans;
    const int ns1 = rules[1].num_spans;
    const int expected_nd = (p0 + 1) * (p1 + 1);

    for (int s1 = 0; s1 < ns1; s1++) {
      const int k1 = span_knot_idx[1][s1];
      for (int s0 = 0; s0 < ns0; s0++) {
        const int k0 = span_knot_idx[0][s0];
        const int el = patch_el_offset + s0 + ns0 * s1;
        fes.GetElementDofs(el, el_dofs);
        MFEM_VERIFY(el_dofs.Size() == expected_nd,
                    "Unexpected 2D NURBS element DoF count.");

        for (int r1 = 0; r1 <= p1; r1++) {
          const int i1 = k1 - p1 + r1;
          for (int r0 = 0; r0 <= p0; r0++) {
            const int i0 = k0 - p0 + r0;
            const int p_idx = i0 + n0 * i1;
            const int g_dof = el_dofs[r0 + (p0 + 1) * r1];
            MFEM_VERIFY(patch_dofs[p_idx] == -1 || patch_dofs[p_idx] == g_dof,
                        "Inconsistent 2D patch-to-global DoF mapping.");
            patch_dofs[p_idx] = g_dof;
          }
        }
      }
    }
  } else {
    const int p0 = rules[0].order;
    const int p1 = rules[1].order;
    const int p2 = rules[2].order;
    const int n0 = rules[0].num_dofs;
    const int n1 = rules[1].num_dofs;
    const int ns0 = rules[0].num_spans;
    const int ns1 = rules[1].num_spans;
    const int ns2 = rules[2].num_spans;
    const int expected_nd = (p0 + 1) * (p1 + 1) * (p2 + 1);

    for (int s2 = 0; s2 < ns2; s2++) {
      const int k2 = span_knot_idx[2][s2];
      for (int s1 = 0; s1 < ns1; s1++) {
        const int k1 = span_knot_idx[1][s1];
        for (int s0 = 0; s0 < ns0; s0++) {
          const int k0 = span_knot_idx[0][s0];
          const int el = patch_el_offset + s0 + ns0 * (s1 + ns1 * s2);
          fes.GetElementDofs(el, el_dofs);
          MFEM_VERIFY(el_dofs.Size() == expected_nd,
                      "Unexpected 3D NURBS element DoF count.");

          for (int r2 = 0; r2 <= p2; r2++) {
            const int i2 = k2 - p2 + r2;
            for (int r1 = 0; r1 <= p1; r1++) {
              const int i1 = k1 - p1 + r1;
              for (int r0 = 0; r0 <= p0; r0++) {
                const int i0 = k0 - p0 + r0;
                const int p_idx = i0 + n0 * (i1 + n1 * i2);
                const int loc = r0 + (p0 + 1) * (r1 + (p1 + 1) * r2);
                const int g_dof = el_dofs[loc];
                MFEM_VERIFY(patch_dofs[p_idx] == -1 ||
                                patch_dofs[p_idx] == g_dof,
                            "Inconsistent 3D patch-to-global DoF "
                            "mapping.");
                patch_dofs[p_idx] = g_dof;
              }
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < num_patch_dofs; i++) {
    MFEM_VERIFY(patch_dofs[i] >= 0,
                "Unmapped patch DoF detected in BuildPatchDofMap.");
  }
}

/// Check if patch_dofs is a bijection onto 0..num_tdofs-1 (and whether it is
/// the identity map).
bool CheckBijectiveDofMap(const Array<int> &patch_dofs, int num_tdofs,
                          bool &is_identity) {
  is_identity = false;
  if (patch_dofs.Size() != num_tdofs) {
    return false;
  }

  std::vector<char> seen(static_cast<std::size_t>(num_tdofs), 0);
  bool identity = true;
  for (int i = 0; i < num_tdofs; i++) {
    const int d = patch_dofs[i];
    if (d < 0 || d >= num_tdofs || seen[static_cast<std::size_t>(d)]) {
      return false;
    }
    seen[static_cast<std::size_t>(d)] = 1;
    if (d != i) {
      identity = false;
    }
  }

  is_identity = identity;
  return true;
}

/// 2D sum-factorized WQ diffusion assembly (O(N_dof * p^3)).
void AssemblePatch2D(const std::vector<WQRule1D> &rules,
                     const std::vector<WQWeights1D> &weights,
                     const WQPatchMetric &metric, const WQOptions &opts,
                     const Array<int> &patch_dofs, bool use_fast_csr,
                     const int *I_csr, int *J_csr, real_t *Data_csr,
                     SparseMatrix &A) {
  const int n0 = rules[0].num_dofs;
  const int n1 = rules[1].num_dofs;
  const int Q0 = rules[0].num_qpts;
  const bool rational = metric.is_rational;
  const bool use_general = (!metric.is_diagonal || rational);
  const real_t drop_tol = opts.drop_tol;

  const int *qb0 = rules[0].q_begin.HostRead();
  const int *qe0 = rules[0].q_end.HostRead();
  const int *jb0 = rules[0].j_begin.HostRead();
  const int *je0 = rules[0].j_end.HostRead();
  const int *off0 = rules[0].trial_offsets.HostRead();
  const real_t *B0 = rules[0].B_vals.HostRead();
  const real_t *D0 = rules[0].D_vals.HostRead();
  const real_t *W00_0 = weights[0].W00.HostRead();
  const real_t *W11_0 = weights[0].W11.HostRead();
  const real_t *W01_0 = weights[0].W01.HostRead();
  const real_t *W10_0 = weights[0].W10.HostRead();

  const int *qb1 = rules[1].q_begin.HostRead();
  const int *qe1 = rules[1].q_end.HostRead();
  const int *jb1 = rules[1].j_begin.HostRead();
  const int *je1 = rules[1].j_end.HostRead();
  const int *off1 = rules[1].trial_offsets.HostRead();
  const real_t *B1 = rules[1].B_vals.HostRead();
  const real_t *D1 = rules[1].D_vals.HostRead();
  const real_t *W00_1 = weights[1].W00.HostRead();
  const real_t *W11_1 = weights[1].W11.HostRead();
  const real_t *W01_1 = weights[1].W01.HostRead();
  const real_t *W10_1 = weights[1].W10.HostRead();

  const real_t *G = metric.G.HostRead();
  const real_t *V = rational ? metric.V.HostRead() : nullptr;
  const real_t *C = rational ? metric.C.HostRead() : nullptr;
  const real_t *pw = rational ? metric.patch_weights.HostRead() : nullptr;
  const int *pdofs = patch_dofs.HostRead();

  std::vector<real_t> C11(static_cast<std::size_t>(Q0));
  std::vector<real_t> C00(static_cast<std::size_t>(Q0));
  std::vector<real_t> C10, C01;
  if (use_general) {
    C10.resize(static_cast<std::size_t>(Q0));
    C01.resize(static_cast<std::size_t>(Q0));
  }

  for (int i1 = 0; i1 < n1; i1++) {
    const int qb_i1 = qb1[i1];
    const int qe_i1 = qe1[i1];
    const int off_i1 = off1[i1];
    const int jb_i1 = jb1[i1];
    const int je_i1 = je1[i1];

    for (int j1 = jb_i1; j1 < je_i1; j1++) {
      const int qb_j1 = qb1[j1];
      const int off_j1 = off1[j1];
      const int q1_lo = std::max(qb_i1, qb_j1);
      const int q1_hi = std::min(qe_i1, qe1[j1]);

      std::fill(C11.begin(), C11.end(), static_cast<real_t>(0.0));
      std::fill(C00.begin(), C00.end(), static_cast<real_t>(0.0));
      if (use_general) {
        std::fill(C10.begin(), C10.end(), static_cast<real_t>(0.0));
        std::fill(C01.begin(), C01.end(), static_cast<real_t>(0.0));
      }

      // Contract direction 1 into 1D buffers along q0.
      for (int q1 = q1_lo; q1 < q1_hi; q1++) {
        const int ii1 = off_i1 + (q1 - qb_i1);
        const int jj1 = off_j1 + (q1 - qb_j1);
        const real_t wb1 = W00_1[ii1] * B1[jj1];
        const real_t wd1 = W11_1[ii1] * D1[jj1];
        const std::size_t q_row =
            static_cast<std::size_t>(Q0) * static_cast<std::size_t>(q1);
        const real_t *g_row = G + static_cast<std::size_t>(3) * q_row;

        if (!use_general) {
          for (int q0 = 0; q0 < Q0; q0++) {
            const std::size_t q0u = static_cast<std::size_t>(q0);
            const std::size_t g3 = static_cast<std::size_t>(3) * q0u;
            C11[q0u] += wb1 * g_row[g3];
            C00[q0u] += wd1 * g_row[g3 + 2];
          }
        } else if (!rational) {
          const real_t w01d1 = W01_1[ii1] * D1[jj1];
          const real_t w10b1 = W10_1[ii1] * B1[jj1];
          for (int q0 = 0; q0 < Q0; q0++) {
            const std::size_t q0u = static_cast<std::size_t>(q0);
            const std::size_t g3 = static_cast<std::size_t>(3) * q0u;
            const real_t g01 = g_row[g3 + 1];
            C11[q0u] += wb1 * g_row[g3];
            C00[q0u] += wd1 * g_row[g3 + 2];
            C10[q0u] += w01d1 * g01;
            C01[q0u] += w10b1 * g01;
          }
        } else {
          const real_t w01d1 = W01_1[ii1] * D1[jj1];
          const real_t w10b1 = W10_1[ii1] * B1[jj1];
          const real_t wmix1 = w10b1 + w01d1;
          const real_t *v_row = V + static_cast<std::size_t>(2) * q_row;
          const real_t *c_row = C + q_row;
          for (int q0 = 0; q0 < Q0; q0++) {
            const std::size_t q0u = static_cast<std::size_t>(q0);
            const std::size_t g3 = static_cast<std::size_t>(3) * q0u;
            const std::size_t v2 = static_cast<std::size_t>(2) * q0u;
            const real_t g01 = g_row[g3 + 1];
            const real_t v0_wb = v_row[v2] * wb1;
            const real_t v1 = v_row[v2 + 1];
            C11[q0u] += wb1 * g_row[g3];
            C00[q0u] += wd1 * g_row[g3 + 2] + v1 * wmix1 + c_row[q0u] * wb1;
            C10[q0u] += w01d1 * g01 + v0_wb;
            C01[q0u] += w10b1 * g01 + v0_wb;
          }
        }
      }

      // Contract direction 0 into matrix entries.
      const int j1_loc = j1 - jb_i1;
      for (int i0 = 0; i0 < n0; i0++) {
        const int i_local = i0 + n0 * i1;
        const int row = pdofs[i_local];
        const real_t wi = rational ? pw[i_local] : static_cast<real_t>(1.0);
        const int qb_i0 = qb0[i0];
        const int qe_i0 = qe0[i0];
        const int off_i0 = off0[i0];
        const int jb_i0 = jb0[i0];
        const int je_i0 = je0[i0];
        const int nj0 = je_i0 - jb_i0;
        const int row_base = use_fast_csr ? (I_csr[row] + j1_loc * nj0) : 0;

        for (int j0 = jb_i0; j0 < je_i0; j0++) {
          const int j_local = j0 + n0 * j1;
          const int col = pdofs[j_local];
          const int qb_j0 = qb0[j0];
          const int off_j0 = off0[j0];
          const int q0_lo = std::max(qb_i0, qb_j0);
          const int q0_hi = std::min(qe_i0, qe0[j0]);

          real_t val = static_cast<real_t>(0.0);
          if (!use_general) {
            for (int q0 = q0_lo; q0 < q0_hi; q0++) {
              const int ii0 = off_i0 + (q0 - qb_i0);
              const int jj0 = off_j0 + (q0 - qb_j0);
              const std::size_t q0u = static_cast<std::size_t>(q0);
              val += C11[q0u] * W11_0[ii0] * D0[jj0] +
                     C00[q0u] * W00_0[ii0] * B0[jj0];
            }
          } else {
            for (int q0 = q0_lo; q0 < q0_hi; q0++) {
              const int ii0 = off_i0 + (q0 - qb_i0);
              const int jj0 = off_j0 + (q0 - qb_j0);
              const std::size_t q0u = static_cast<std::size_t>(q0);
              val += C11[q0u] * W11_0[ii0] * D0[jj0] +
                     C00[q0u] * W00_0[ii0] * B0[jj0] +
                     C10[q0u] * W10_0[ii0] * B0[jj0] +
                     C01[q0u] * W01_0[ii0] * D0[jj0];
            }
          }

          if (rational) {
            val *= wi * pw[j_local];
          }

          const bool dropped = (drop_tol > static_cast<real_t>(0.0) &&
                                row != col && std::abs(val) <= drop_tol);
          if (use_fast_csr) {
            const int pos = row_base + (j0 - jb_i0);
            J_csr[pos] = col;
            Data_csr[pos] = dropped ? static_cast<real_t>(0.0) : val;
          } else if (!dropped) {
            A.Add(row, col, val);
          }
        }
      }
    }
  }
}

/// 3D sum-factorized WQ diffusion assembly (O(N_dof * p^4)).
void AssemblePatch3D(const std::vector<WQRule1D> &rules,
                     const std::vector<WQWeights1D> &weights,
                     const WQPatchMetric &metric, const WQOptions &opts,
                     const Array<int> &patch_dofs, bool use_fast_csr,
                     const int *I_csr, int *J_csr, real_t *Data_csr,
                     SparseMatrix &A) {
  const int n0 = rules[0].num_dofs;
  const int n1 = rules[1].num_dofs;
  const int n2 = rules[2].num_dofs;
  const int Q0 = rules[0].num_qpts;
  const int Q1 = rules[1].num_qpts;
  const std::size_t Q01 =
      static_cast<std::size_t>(Q0) * static_cast<std::size_t>(Q1);
  const bool rational = metric.is_rational;
  const bool use_general = (!metric.is_diagonal || rational);
  const real_t drop_tol = opts.drop_tol;

  const int *qb0 = rules[0].q_begin.HostRead();
  const int *qe0 = rules[0].q_end.HostRead();
  const int *jb0 = rules[0].j_begin.HostRead();
  const int *je0 = rules[0].j_end.HostRead();
  const int *off0 = rules[0].trial_offsets.HostRead();
  const real_t *B0 = rules[0].B_vals.HostRead();
  const real_t *D0 = rules[0].D_vals.HostRead();
  const real_t *W00_0 = weights[0].W00.HostRead();
  const real_t *W11_0 = weights[0].W11.HostRead();
  const real_t *W01_0 = weights[0].W01.HostRead();
  const real_t *W10_0 = weights[0].W10.HostRead();

  const int *qb1 = rules[1].q_begin.HostRead();
  const int *qe1 = rules[1].q_end.HostRead();
  const int *jb1 = rules[1].j_begin.HostRead();
  const int *je1 = rules[1].j_end.HostRead();
  const int *off1 = rules[1].trial_offsets.HostRead();
  const real_t *B1 = rules[1].B_vals.HostRead();
  const real_t *D1 = rules[1].D_vals.HostRead();
  const real_t *W00_1 = weights[1].W00.HostRead();
  const real_t *W11_1 = weights[1].W11.HostRead();
  const real_t *W01_1 = weights[1].W01.HostRead();
  const real_t *W10_1 = weights[1].W10.HostRead();

  const int *qb2 = rules[2].q_begin.HostRead();
  const int *qe2 = rules[2].q_end.HostRead();
  const int *jb2 = rules[2].j_begin.HostRead();
  const int *je2 = rules[2].j_end.HostRead();
  const int *off2 = rules[2].trial_offsets.HostRead();
  const real_t *B2 = rules[2].B_vals.HostRead();
  const real_t *D2 = rules[2].D_vals.HostRead();
  const real_t *W00_2 = weights[2].W00.HostRead();
  const real_t *W11_2 = weights[2].W11.HostRead();
  const real_t *W01_2 = weights[2].W01.HostRead();
  const real_t *W10_2 = weights[2].W10.HostRead();

  const real_t *G = metric.G.HostRead();
  const real_t *V = rational ? metric.V.HostRead() : nullptr;
  const real_t *C = rational ? metric.C.HostRead() : nullptr;
  const real_t *pw = rational ? metric.patch_weights.HostRead() : nullptr;
  const int *pdofs = patch_dofs.HostRead();

  // Stage-1 2D slices (size Q0 * Q1).
  std::vector<real_t> C2_00(Q01), C2_11(Q01), C2_22(Q01);
  std::vector<real_t> C2_01, C2_02, C2_20, C2_12, C2_21;
  if (use_general) {
    C2_01.resize(Q01);
    C2_02.resize(Q01);
    C2_20.resize(Q01);
    C2_12.resize(Q01);
    C2_21.resize(Q01);
  }

  // Stage-2 1D buffers (size Q0).
  std::vector<real_t> C1_11(static_cast<std::size_t>(Q0));
  std::vector<real_t> C1_00(static_cast<std::size_t>(Q0));
  std::vector<real_t> C1_10, C1_01;
  if (use_general) {
    C1_10.resize(static_cast<std::size_t>(Q0));
    C1_01.resize(static_cast<std::size_t>(Q0));
  }

  for (int i2 = 0; i2 < n2; i2++) {
    const int qb_i2 = qb2[i2];
    const int qe_i2 = qe2[i2];
    const int off_i2 = off2[i2];
    const int jb_i2 = jb2[i2];
    const int je_i2 = je2[i2];

    for (int j2 = jb_i2; j2 < je_i2; j2++) {
      const int qb_j2 = qb2[j2];
      const int off_j2 = off2[j2];
      const int q2_lo = std::max(qb_i2, qb_j2);
      const int q2_hi = std::min(qe_i2, qe2[j2]);

      std::fill(C2_00.begin(), C2_00.end(), static_cast<real_t>(0.0));
      std::fill(C2_11.begin(), C2_11.end(), static_cast<real_t>(0.0));
      std::fill(C2_22.begin(), C2_22.end(), static_cast<real_t>(0.0));
      if (use_general) {
        std::fill(C2_01.begin(), C2_01.end(), static_cast<real_t>(0.0));
        std::fill(C2_02.begin(), C2_02.end(), static_cast<real_t>(0.0));
        std::fill(C2_20.begin(), C2_20.end(), static_cast<real_t>(0.0));
        std::fill(C2_12.begin(), C2_12.end(), static_cast<real_t>(0.0));
        std::fill(C2_21.begin(), C2_21.end(), static_cast<real_t>(0.0));
      }

      // Stage 1: contract direction 2 over q2 into 2D slices.
      for (int q2 = q2_lo; q2 < q2_hi; q2++) {
        const int ii2 = off_i2 + (q2 - qb_i2);
        const int jj2 = off_j2 + (q2 - qb_j2);
        const real_t wb2 = W00_2[ii2] * B2[jj2];
        const real_t wd2 = W11_2[ii2] * D2[jj2];
        const std::size_t q_slice = Q01 * static_cast<std::size_t>(q2);
        const real_t *g_slice = G + static_cast<std::size_t>(6) * q_slice;

        if (!use_general) {
          for (std::size_t idx = 0; idx < Q01; idx++) {
            const std::size_t g6 = static_cast<std::size_t>(6) * idx;
            C2_00[idx] += wb2 * g_slice[g6];
            C2_11[idx] += wb2 * g_slice[g6 + 3];
            C2_22[idx] += wd2 * g_slice[g6 + 5];
          }
        } else if (!rational) {
          const real_t w01d2 = W01_2[ii2] * D2[jj2];
          const real_t w10b2 = W10_2[ii2] * B2[jj2];
          for (std::size_t idx = 0; idx < Q01; idx++) {
            const std::size_t g6 = static_cast<std::size_t>(6) * idx;
            const real_t g02 = g_slice[g6 + 2];
            const real_t g12 = g_slice[g6 + 4];
            C2_00[idx] += wb2 * g_slice[g6];
            C2_11[idx] += wb2 * g_slice[g6 + 3];
            C2_22[idx] += wd2 * g_slice[g6 + 5];
            C2_01[idx] += wb2 * g_slice[g6 + 1];
            C2_02[idx] += w01d2 * g02;
            C2_20[idx] += w10b2 * g02;
            C2_12[idx] += w01d2 * g12;
            C2_21[idx] += w10b2 * g12;
          }
        } else {
          const real_t w01d2 = W01_2[ii2] * D2[jj2];
          const real_t w10b2 = W10_2[ii2] * B2[jj2];
          const real_t wmix2 = w10b2 + w01d2;
          const real_t *v_slice = V + static_cast<std::size_t>(3) * q_slice;
          const real_t *c_slice = C + q_slice;
          for (std::size_t idx = 0; idx < Q01; idx++) {
            const std::size_t g6 = static_cast<std::size_t>(6) * idx;
            const std::size_t v3 = static_cast<std::size_t>(3) * idx;
            const real_t g02 = g_slice[g6 + 2];
            const real_t g12 = g_slice[g6 + 4];
            const real_t v0_wb = v_slice[v3] * wb2;
            const real_t v1_wb = v_slice[v3 + 1] * wb2;
            const real_t v2 = v_slice[v3 + 2];
            C2_00[idx] += wb2 * g_slice[g6];
            C2_11[idx] += wb2 * g_slice[g6 + 3];
            C2_22[idx] +=
                wd2 * g_slice[g6 + 5] + v2 * wmix2 + c_slice[idx] * wb2;
            C2_01[idx] += wb2 * g_slice[g6 + 1];
            C2_02[idx] += w01d2 * g02 + v0_wb;
            C2_20[idx] += w10b2 * g02 + v0_wb;
            C2_12[idx] += w01d2 * g12 + v1_wb;
            C2_21[idx] += w10b2 * g12 + v1_wb;
          }
        }
      }

      const int j2_loc = j2 - jb_i2;

      // Stage 2: contract direction 1 over q1 into 1D buffers.
      for (int i1 = 0; i1 < n1; i1++) {
        const int qb_i1 = qb1[i1];
        const int qe_i1 = qe1[i1];
        const int off_i1 = off1[i1];
        const int jb_i1 = jb1[i1];
        const int je_i1 = je1[i1];
        const int nj1 = je_i1 - jb_i1;

        for (int j1 = jb_i1; j1 < je_i1; j1++) {
          const int qb_j1 = qb1[j1];
          const int off_j1 = off1[j1];
          const int q1_lo = std::max(qb_i1, qb_j1);
          const int q1_hi = std::min(qe_i1, qe1[j1]);

          std::fill(C1_11.begin(), C1_11.end(), static_cast<real_t>(0.0));
          std::fill(C1_00.begin(), C1_00.end(), static_cast<real_t>(0.0));
          if (use_general) {
            std::fill(C1_10.begin(), C1_10.end(), static_cast<real_t>(0.0));
            std::fill(C1_01.begin(), C1_01.end(), static_cast<real_t>(0.0));
          }

          for (int q1 = q1_lo; q1 < q1_hi; q1++) {
            const int ii1 = off_i1 + (q1 - qb_i1);
            const int jj1 = off_j1 + (q1 - qb_j1);
            const real_t wb1 = W00_1[ii1] * B1[jj1];
            const real_t wd1 = W11_1[ii1] * D1[jj1];
            const std::size_t base =
                static_cast<std::size_t>(Q0) * static_cast<std::size_t>(q1);

            if (!use_general) {
              for (int q0 = 0; q0 < Q0; q0++) {
                const std::size_t q0u = static_cast<std::size_t>(q0);
                const std::size_t idx = base + q0u;
                C1_11[q0u] += C2_00[idx] * wb1;
                C1_00[q0u] += C2_11[idx] * wd1 + C2_22[idx] * wb1;
              }
            } else {
              const real_t w01d1 = W01_1[ii1] * D1[jj1];
              const real_t w10b1 = W10_1[ii1] * B1[jj1];
              for (int q0 = 0; q0 < Q0; q0++) {
                const std::size_t q0u = static_cast<std::size_t>(q0);
                const std::size_t idx = base + q0u;
                const real_t c01 = C2_01[idx];
                C1_11[q0u] += C2_00[idx] * wb1;
                C1_00[q0u] += C2_11[idx] * wd1 + C2_22[idx] * wb1 +
                              C2_12[idx] * w10b1 + C2_21[idx] * w01d1;
                C1_10[q0u] += c01 * w01d1 + C2_02[idx] * wb1;
                C1_01[q0u] += c01 * w10b1 + C2_20[idx] * wb1;
              }
            }
          }

          // Stage 3: contract direction 0 into matrix entries.
          const int j1_loc = j1 - jb_i1;
          const int j21_loc = j2_loc * nj1 + j1_loc;
          const int row_12 = n0 * (i1 + n1 * i2);
          const int col_12 = n0 * (j1 + n1 * j2);

          for (int i0 = 0; i0 < n0; i0++) {
            const int i_local = i0 + row_12;
            const int row = pdofs[i_local];
            const real_t wi = rational ? pw[i_local] : static_cast<real_t>(1.0);
            const int qb_i0 = qb0[i0];
            const int qe_i0 = qe0[i0];
            const int off_i0 = off0[i0];
            const int jb_i0 = jb0[i0];
            const int je_i0 = je0[i0];
            const int nj0 = je_i0 - jb_i0;
            const int row_base =
                use_fast_csr ? (I_csr[row] + j21_loc * nj0) : 0;

            for (int j0 = jb_i0; j0 < je_i0; j0++) {
              const int j_local = j0 + col_12;
              const int col = pdofs[j_local];
              const int qb_j0 = qb0[j0];
              const int off_j0 = off0[j0];
              const int q0_lo = std::max(qb_i0, qb_j0);
              const int q0_hi = std::min(qe_i0, qe0[j0]);

              real_t val = static_cast<real_t>(0.0);
              if (!use_general) {
                for (int q0 = q0_lo; q0 < q0_hi; q0++) {
                  const int ii0 = off_i0 + (q0 - qb_i0);
                  const int jj0 = off_j0 + (q0 - qb_j0);
                  const std::size_t q0u = static_cast<std::size_t>(q0);
                  val += C1_11[q0u] * W11_0[ii0] * D0[jj0] +
                         C1_00[q0u] * W00_0[ii0] * B0[jj0];
                }
              } else {
                for (int q0 = q0_lo; q0 < q0_hi; q0++) {
                  const int ii0 = off_i0 + (q0 - qb_i0);
                  const int jj0 = off_j0 + (q0 - qb_j0);
                  const std::size_t q0u = static_cast<std::size_t>(q0);
                  val += C1_11[q0u] * W11_0[ii0] * D0[jj0] +
                         C1_00[q0u] * W00_0[ii0] * B0[jj0] +
                         C1_10[q0u] * W10_0[ii0] * B0[jj0] +
                         C1_01[q0u] * W01_0[ii0] * D0[jj0];
                }
              }

              if (rational) {
                val *= wi * pw[j_local];
              }

              const bool dropped = (drop_tol > static_cast<real_t>(0.0) &&
                                    row != col && std::abs(val) <= drop_tol);
              if (use_fast_csr) {
                const int pos = row_base + (j0 - jb_i0);
                J_csr[pos] = col;
                Data_csr[pos] = dropped ? static_cast<real_t>(0.0) : val;
              } else if (!dropped) {
                A.Add(row, col, val);
              }
            }
          }
        }
      }
    }
  }
}

} // namespace

void EvaluatePatchMetric(const Mesh &mesh, int patch,
                         const std::vector<WQRule1D> &rules, Coefficient &coeff,
                         WQPatchMetric &metric) {
  const int dim = mesh.Dimension();
  MFEM_VERIFY(dim == 2 || dim == 3,
              "EvaluatePatchMetric only supports 2D and 3D meshes.");
  MFEM_VERIFY(mesh.NURBSExt() != nullptr,
              "EvaluatePatchMetric requires a NURBS mesh.");
  MFEM_VERIFY(static_cast<int>(rules.size()) == dim,
              "Number of 1D WQ rules must match mesh dimension.");

  int patch_el_offset = 0;
  int patch_ne = 0;
  Array<const KnotVector *> patch_kv;
  GetPatchInfoFromNURBSExt(mesh.NURBSExt(), mesh.GetNE(), dim, patch,
                           patch_el_offset, patch_ne, patch_kv);

  for (int d = 0; d < dim; d++) {
    MFEM_VERIFY(patch_kv[d]->GetNE() == rules[d].num_spans,
                "WQRule1D span count does not match mesh patch.");
  }

  ExtractPatchWeights(mesh, patch_kv, patch_el_offset, metric.is_rational,
                      metric.patch_weights);

  if (metric.is_rational) {
    for (int d = 0; d < dim; d++) {
      MFEM_VERIFY(patch_kv[d]->GetOrder() == rules[d].order &&
                      patch_kv[d]->GetNCP() == rules[d].num_dofs,
                  "Rational NURBS Weighted Quadrature requires WQRule1D "
                  "to match the mesh patch KnotVectors.");
    }
  }

  metric.dim = dim;
  metric.q_sizes.SetSize(dim);
  std::size_t total_qpts = 1;
  for (int d = 0; d < dim; d++) {
    metric.q_sizes[d] = rules[d].num_qpts;
    total_qpts *= static_cast<std::size_t>(rules[d].num_qpts);
  }

  const int stride = (dim == 2) ? 3 : 6;
  metric.G.SetSize(static_cast<int>(total_qpts) * stride);
  real_t *G_data = metric.G.HostWrite();

  real_t max_diag = static_cast<real_t>(0.0);
  real_t max_offdiag = static_cast<real_t>(0.0);

  IsoparametricTransformation el_tr;
  IntegrationPoint ip;

  if (dim == 2) {
    const int Q0 = rules[0].num_qpts;
    const int ns0 = rules[0].num_spans;
    const int ns1 = rules[1].num_spans;
    const real_t *lxi0 = rules[0].local_xi.HostRead();
    const real_t *lxi1 = rules[1].local_xi.HostRead();
    const real_t *h0_arr = rules[0].span_lengths.HostRead();
    const real_t *h1_arr = rules[1].span_lengths.HostRead();

    for (int s1 = 0; s1 < ns1; s1++) {
      const real_t du1 = h1_arr[s1];
      const int q1_b = rules[1].span_q_begin[s1];
      const int q1_e = rules[1].span_q_end[s1];

      for (int s0 = 0; s0 < ns0; s0++) {
        const real_t du0 = h0_arr[s0];
        const int q0_b = rules[0].span_q_begin[s0];
        const int q0_e = rules[0].span_q_end[s0];
        const int el = patch_el_offset + s0 + ns0 * s1;

        mesh.GetElementTransformation(el, &el_tr);
        const real_t inv_vol = static_cast<real_t>(1.0) / (du0 * du1);
        const real_t s00 = du0 * du0 * inv_vol;
        const real_t s01 = du0 * du1 * inv_vol;
        const real_t s11 = du1 * du1 * inv_vol;

        for (int q1 = q1_b; q1 < q1_e; q1++) {
          ip.y = lxi1[q1];
          for (int q0 = q0_b; q0 < q0_e; q0++) {
            ip.x = lxi0[q0];
            el_tr.SetIntPoint(&ip);

            const real_t alpha = coeff.Eval(el_tr, ip);
            const real_t detJ = std::abs(el_tr.Weight());
            const DenseMatrix &invJ = el_tr.InverseJacobian();
            const int sdim = invJ.Width();

            real_t m00 = static_cast<real_t>(0.0);
            real_t m01 = static_cast<real_t>(0.0);
            real_t m11 = static_cast<real_t>(0.0);
            for (int k = 0; k < sdim; k++) {
              const real_t j0k = invJ(0, k);
              const real_t j1k = invJ(1, k);
              m00 += j0k * j0k;
              m01 += j0k * j1k;
              m11 += j1k * j1k;
            }

            const real_t aw = alpha * detJ;
            const real_t g00 = aw * s00 * m00;
            const real_t g01 = aw * s01 * m01;
            const real_t g11 = aw * s11 * m11;

            const std::size_t q_idx =
                static_cast<std::size_t>(q0) +
                static_cast<std::size_t>(Q0) * static_cast<std::size_t>(q1);
            real_t *dst = G_data + static_cast<std::size_t>(3) * q_idx;
            dst[0] = g00;
            dst[1] = g01;
            dst[2] = g11;

            max_diag = std::max(max_diag, std::abs(g00));
            max_diag = std::max(max_diag, std::abs(g11));
            max_offdiag = std::max(max_offdiag, std::abs(g01));
          }
        }
      }
    }
  } else {
    const int Q0 = rules[0].num_qpts;
    const int Q1 = rules[1].num_qpts;
    const int ns0 = rules[0].num_spans;
    const int ns1 = rules[1].num_spans;
    const int ns2 = rules[2].num_spans;
    const real_t *lxi0 = rules[0].local_xi.HostRead();
    const real_t *lxi1 = rules[1].local_xi.HostRead();
    const real_t *lxi2 = rules[2].local_xi.HostRead();
    const real_t *h0_arr = rules[0].span_lengths.HostRead();
    const real_t *h1_arr = rules[1].span_lengths.HostRead();
    const real_t *h2_arr = rules[2].span_lengths.HostRead();

    for (int s2 = 0; s2 < ns2; s2++) {
      const real_t du2 = h2_arr[s2];
      const int q2_b = rules[2].span_q_begin[s2];
      const int q2_e = rules[2].span_q_end[s2];

      for (int s1 = 0; s1 < ns1; s1++) {
        const real_t du1 = h1_arr[s1];
        const int q1_b = rules[1].span_q_begin[s1];
        const int q1_e = rules[1].span_q_end[s1];

        for (int s0 = 0; s0 < ns0; s0++) {
          const real_t du0 = h0_arr[s0];
          const int q0_b = rules[0].span_q_begin[s0];
          const int q0_e = rules[0].span_q_end[s0];
          const int el = patch_el_offset + s0 + ns0 * (s1 + ns1 * s2);

          mesh.GetElementTransformation(el, &el_tr);
          const real_t inv_vol = static_cast<real_t>(1.0) / (du0 * du1 * du2);
          const real_t s00 = du0 * du0 * inv_vol;
          const real_t s01 = du0 * du1 * inv_vol;
          const real_t s02 = du0 * du2 * inv_vol;
          const real_t s11 = du1 * du1 * inv_vol;
          const real_t s12 = du1 * du2 * inv_vol;
          const real_t s22 = du2 * du2 * inv_vol;

          for (int q2 = q2_b; q2 < q2_e; q2++) {
            ip.z = lxi2[q2];
            for (int q1 = q1_b; q1 < q1_e; q1++) {
              ip.y = lxi1[q1];
              for (int q0 = q0_b; q0 < q0_e; q0++) {
                ip.x = lxi0[q0];
                el_tr.SetIntPoint(&ip);

                const real_t alpha = coeff.Eval(el_tr, ip);
                const real_t detJ = std::abs(el_tr.Weight());
                const DenseMatrix &invJ = el_tr.InverseJacobian();
                const int sdim = invJ.Width();

                real_t m00 = static_cast<real_t>(0.0);
                real_t m01 = static_cast<real_t>(0.0);
                real_t m02 = static_cast<real_t>(0.0);
                real_t m11 = static_cast<real_t>(0.0);
                real_t m12 = static_cast<real_t>(0.0);
                real_t m22 = static_cast<real_t>(0.0);
                for (int k = 0; k < sdim; k++) {
                  const real_t j0k = invJ(0, k);
                  const real_t j1k = invJ(1, k);
                  const real_t j2k = invJ(2, k);
                  m00 += j0k * j0k;
                  m01 += j0k * j1k;
                  m02 += j0k * j2k;
                  m11 += j1k * j1k;
                  m12 += j1k * j2k;
                  m22 += j2k * j2k;
                }

                const real_t aw = alpha * detJ;
                const real_t g00 = aw * s00 * m00;
                const real_t g01 = aw * s01 * m01;
                const real_t g02 = aw * s02 * m02;
                const real_t g11 = aw * s11 * m11;
                const real_t g12 = aw * s12 * m12;
                const real_t g22 = aw * s22 * m22;

                const std::size_t q_idx =
                    static_cast<std::size_t>(q0) +
                    static_cast<std::size_t>(Q0) *
                        (static_cast<std::size_t>(q1) +
                         static_cast<std::size_t>(Q1) *
                             static_cast<std::size_t>(q2));
                real_t *dst = G_data + static_cast<std::size_t>(6) * q_idx;
                dst[0] = g00;
                dst[1] = g01;
                dst[2] = g02;
                dst[3] = g11;
                dst[4] = g12;
                dst[5] = g22;

                max_diag = std::max(max_diag, std::abs(g00));
                max_diag = std::max(max_diag, std::abs(g11));
                max_diag = std::max(max_diag, std::abs(g22));
                max_offdiag = std::max(max_offdiag, std::abs(g01));
                max_offdiag = std::max(max_offdiag, std::abs(g02));
                max_offdiag = std::max(max_offdiag, std::abs(g12));
              }
            }
          }
        }
      }
    }
  }

  metric.is_diagonal = (max_offdiag <= static_cast<real_t>(1e-14) * max_diag);

  if (metric.is_rational) {
    std::vector<real_t> W;
    std::vector<std::vector<real_t>> dW;
    EvaluateNURBSDenominator(rules, metric.patch_weights, W, dW);

    metric.V.SetSize(static_cast<int>(total_qpts) * dim);
    metric.C.SetSize(static_cast<int>(total_qpts));
    real_t *V_data = metric.V.HostWrite();
    real_t *C_data = metric.C.HostWrite();

    if (dim == 2) {
      const real_t *W0 = dW[0].data();
      const real_t *W1 = dW[1].data();
      for (std::size_t q = 0; q < total_qpts; q++) {
        const real_t wq = W[q];
        MFEM_VERIFY(wq > static_cast<real_t>(0.0),
                    "Non-positive NURBS denominator W(q) at WQ point.");
        const real_t inv_W = static_cast<real_t>(1.0) / wq;
        const real_t inv_W2 = inv_W * inv_W;
        const real_t gam0 = W0[q] * inv_W;
        const real_t gam1 = W1[q] * inv_W;

        real_t *g_ptr = G_data + static_cast<std::size_t>(3) * q;
        g_ptr[0] *= inv_W2;
        g_ptr[1] *= inv_W2;
        g_ptr[2] *= inv_W2;

        const real_t v0 = -(g_ptr[0] * gam0 + g_ptr[1] * gam1);
        const real_t v1 = -(g_ptr[1] * gam0 + g_ptr[2] * gam1);
        V_data[static_cast<std::size_t>(2) * q] = v0;
        V_data[static_cast<std::size_t>(2) * q + 1] = v1;
        C_data[q] = -(v0 * gam0 + v1 * gam1);
      }
    } else {
      const real_t *W0 = dW[0].data();
      const real_t *W1 = dW[1].data();
      const real_t *W2 = dW[2].data();
      for (std::size_t q = 0; q < total_qpts; q++) {
        const real_t wq = W[q];
        MFEM_VERIFY(wq > static_cast<real_t>(0.0),
                    "Non-positive NURBS denominator W(q) at WQ point.");
        const real_t inv_W = static_cast<real_t>(1.0) / wq;
        const real_t inv_W2 = inv_W * inv_W;
        const real_t gam0 = W0[q] * inv_W;
        const real_t gam1 = W1[q] * inv_W;
        const real_t gam2 = W2[q] * inv_W;

        real_t *g_ptr = G_data + static_cast<std::size_t>(6) * q;
        g_ptr[0] *= inv_W2;
        g_ptr[1] *= inv_W2;
        g_ptr[2] *= inv_W2;
        g_ptr[3] *= inv_W2;
        g_ptr[4] *= inv_W2;
        g_ptr[5] *= inv_W2;

        const real_t v0 =
            -(g_ptr[0] * gam0 + g_ptr[1] * gam1 + g_ptr[2] * gam2);
        const real_t v1 =
            -(g_ptr[1] * gam0 + g_ptr[3] * gam1 + g_ptr[4] * gam2);
        const real_t v2 =
            -(g_ptr[2] * gam0 + g_ptr[4] * gam1 + g_ptr[5] * gam2);
        V_data[static_cast<std::size_t>(3) * q] = v0;
        V_data[static_cast<std::size_t>(3) * q + 1] = v1;
        V_data[static_cast<std::size_t>(3) * q + 2] = v2;
        C_data[q] = -(v0 * gam0 + v1 * gam1 + v2 * gam2);
      }
    }
  } else {
    metric.V.SetSize(0);
    metric.C.SetSize(0);
    metric.patch_weights.SetSize(0);
  }
}

void AssemblePatchDiffusionWQ(const FiniteElementSpace &fes, int patch,
                              const std::vector<WQRule1D> &rules,
                              const std::vector<WQWeights1D> &weights,
                              const WQPatchMetric &metric,
                              const WQOptions &opts, SparseMatrix &A) {
  const Mesh *mesh = fes.GetMesh();
  MFEM_VERIFY(mesh != nullptr && mesh->NURBSExt() != nullptr,
              "AssemblePatchDiffusionWQ requires a NURBS Mesh.");
  MFEM_VERIFY(fes.GetVDim() == 1,
              "AssemblePatchDiffusionWQ only supports scalar spaces.");
  MFEM_VERIFY(fes.GetProlongationMatrix() == nullptr,
              "AssemblePatchDiffusionWQ requires an unconstrained serial "
              "NURBS space (GetProlongationMatrix() == nullptr).");

  const int dim = mesh->Dimension();
  MFEM_VERIFY(dim == 2 || dim == 3,
              "AssemblePatchDiffusionWQ only supports 2D and 3D.");
  MFEM_VERIFY(static_cast<int>(rules.size()) == dim &&
                  static_cast<int>(weights.size()) == dim,
              "1D rule/weight vector sizes must match spatial dimension.");
  MFEM_VERIFY(metric.dim == dim && metric.q_sizes.Size() == dim,
              "WQPatchMetric dimension mismatch.");

  int total_qpts = 1;
  int num_patch_dofs = 1;
  for (int d = 0; d < dim; d++) {
    MFEM_VERIFY(metric.q_sizes[d] == rules[d].num_qpts,
                "WQPatchMetric q_sizes mismatch with WQRule1D.");
    MFEM_VERIFY(weights[d].offsets.Size() == rules[d].num_dofs + 1 &&
                    weights[d].W00.Size() == rules[d].B_vals.Size() &&
                    weights[d].W11.Size() == rules[d].B_vals.Size() &&
                    weights[d].W01.Size() == rules[d].B_vals.Size() &&
                    weights[d].W10.Size() == rules[d].B_vals.Size(),
                "WQWeights1D table sizes mismatch with WQRule1D.");
    total_qpts *= rules[d].num_qpts;
    num_patch_dofs *= rules[d].num_dofs;
  }
  const int stride = (dim == 2) ? 3 : 6;
  MFEM_VERIFY(metric.G.Size() == total_qpts * stride,
              "WQPatchMetric::G size mismatch.");

  if (metric.is_rational) {
    const NURBSExtension *fes_ext = fes.GetNURBSext();
    bool shares_weights = (fes_ext == nullptr || fes_ext == mesh->NURBSExt());
    if (!shares_weights) {
      const Vector &mw = mesh->NURBSExt()->GetWeights();
      const Vector &fw = fes_ext->GetWeights();
      if (mw.Size() == fw.Size()) {
        const real_t *mw_d = mw.HostRead();
        const real_t *fw_d = fw.HostRead();
        real_t max_w = static_cast<real_t>(0.0);
        for (int k = 0; k < mw.Size(); k++) {
          max_w = std::max(max_w, std::abs(mw_d[k]));
        }
        const real_t tol = static_cast<real_t>(1e-12) * max_w;
        shares_weights = true;
        for (int k = 0; k < mw.Size(); k++) {
          if (std::abs(mw_d[k] - fw_d[k]) > tol) {
            shares_weights = false;
            break;
          }
        }
      }
    }
    MFEM_VERIFY(shares_weights,
                "Rational NURBS Weighted Quadrature requires the "
                "FiniteElementSpace to share the mesh NURBS weights.");
    MFEM_VERIFY(metric.V.Size() == total_qpts * dim,
                "WQPatchMetric::V size mismatch.");
    MFEM_VERIFY(metric.C.Size() == total_qpts,
                "WQPatchMetric::C size mismatch.");
    MFEM_VERIFY(metric.patch_weights.Size() == num_patch_dofs,
                "WQPatchMetric::patch_weights size mismatch.");
  }

  int patch_el_offset = 0;
  int patch_ne = 0;
  Array<const KnotVector *> patch_kvs;
  GetPatchInfoFromFESpace(fes, dim, patch, patch_el_offset, patch_ne,
                          patch_kvs);

  Array<int> patch_dofs;
  BuildPatchDofMap(fes, patch_kvs, rules, patch_el_offset, patch_dofs);

  const int num_patches = mesh->NURBSExt()->GetNP();
  const int num_tdofs = fes.GetTrueVSize();
  const bool a_no_entries =
      A.Empty() || (!A.Finalized() && A.NumNonZeroElems() == 0);

  bool is_identity_map = false;
  const bool use_fast_csr =
      (num_patches == 1) && a_no_entries &&
      CheckBijectiveDofMap(patch_dofs, num_tdofs, is_identity_map);

  int *I_csr = nullptr;
  int *J_csr = nullptr;
  real_t *Data_csr = nullptr;

  if (use_fast_csr) {
    I_csr = new int[num_tdofs + 1];
    I_csr[0] = 0;
    for (int r = 1; r <= num_tdofs; r++) {
      I_csr[r] = 0;
    }

    if (dim == 2) {
      const int n0 = rules[0].num_dofs;
      const int n1 = rules[1].num_dofs;
      for (int i1 = 0; i1 < n1; i1++) {
        const int nj1 = rules[1].j_end[i1] - rules[1].j_begin[i1];
        for (int i0 = 0; i0 < n0; i0++) {
          const int nj0 = rules[0].j_end[i0] - rules[0].j_begin[i0];
          const int row = patch_dofs[i0 + n0 * i1];
          I_csr[row + 1] = nj0 * nj1;
        }
      }
    } else {
      const int n0 = rules[0].num_dofs;
      const int n1 = rules[1].num_dofs;
      const int n2 = rules[2].num_dofs;
      for (int i2 = 0; i2 < n2; i2++) {
        const int nj2 = rules[2].j_end[i2] - rules[2].j_begin[i2];
        for (int i1 = 0; i1 < n1; i1++) {
          const int nj1 = rules[1].j_end[i1] - rules[1].j_begin[i1];
          const int nj12 = nj1 * nj2;
          for (int i0 = 0; i0 < n0; i0++) {
            const int nj0 = rules[0].j_end[i0] - rules[0].j_begin[i0];
            const int row = patch_dofs[i0 + n0 * (i1 + n1 * i2)];
            I_csr[row + 1] = nj0 * nj12;
          }
        }
      }
    }

    for (int r = 0; r < num_tdofs; r++) {
      I_csr[r + 1] += I_csr[r];
    }

    const int nnz = I_csr[num_tdofs];
    J_csr = new int[nnz];
    Data_csr = new real_t[nnz];
  } else {
    if (A.Empty()) {
      SparseMatrix tmp(num_tdofs);
      A.Swap(tmp);
    }
    MFEM_VERIFY(A.Height() == num_tdofs && A.Width() == num_tdofs,
                "Target SparseMatrix dimensions do not match FESpace.");
    MFEM_VERIFY(!A.Finalized(),
                "Multi-patch WQ assembly requires an unfinalized "
                "SparseMatrix.");
  }

  if (dim == 2) {
    AssemblePatch2D(rules, weights, metric, opts, patch_dofs, use_fast_csr,
                    I_csr, J_csr, Data_csr, A);
  } else {
    AssemblePatch3D(rules, weights, metric, opts, patch_dofs, use_fast_csr,
                    I_csr, J_csr, Data_csr, A);
  }

  if (use_fast_csr) {
    SparseMatrix csr(I_csr, J_csr, Data_csr, num_tdofs, num_tdofs, true, true,
                     is_identity_map);
    if (!is_identity_map) {
      csr.SortColumnIndices();
    }
    A.Swap(csr);
  }
}

void AssembleDiffusionWQ(const FiniteElementSpace &fes, Coefficient &coeff,
                         const Array<int> &ess_tdof_list, Vector &x, Vector &b,
                         SparseMatrix &A, Vector &X, Vector &B,
                         const WQOptions &opts, WQTimings *timings) {
  using Clock = std::chrono::steady_clock;
  using Seconds = std::chrono::duration<double>;

  const Mesh *mesh = fes.GetMesh();
  MFEM_VERIFY(mesh != nullptr, "FiniteElementSpace has a null Mesh.");
  MFEM_VERIFY(mesh->NURBSExt() != nullptr,
              "AssembleDiffusionWQ requires a NURBS Mesh.");
  MFEM_VERIFY(fes.GetVDim() == 1,
              "AssembleDiffusionWQ only supports scalar (vdim == 1) spaces.");
  MFEM_VERIFY(fes.GetProlongationMatrix() == nullptr,
              "AssembleDiffusionWQ requires an unconstrained serial NURBS "
              "space.");

  const int dim = mesh->Dimension();
  MFEM_VERIFY(dim == 2 || dim == 3,
              "AssembleDiffusionWQ only supports 2D and 3D meshes.");

  const int num_patches = mesh->NURBSExt()->GetNP();
  const int num_tdofs = fes.GetTrueVSize();

  MFEM_VERIFY(x.Size() == num_tdofs && b.Size() == num_tdofs,
              "Input vectors x and b must match fes.GetTrueVSize().");

  if (timings != nullptr) {
    *timings = WQTimings();
  }

  if (num_patches == 1) {
    A.Clear();
  } else {
    SparseMatrix tmp(num_tdofs);
    A.Swap(tmp);
  }

  Array<const KnotVector *> kvs;
  std::vector<WQRule1D> rules(static_cast<std::size_t>(dim));
  std::vector<WQWeights1D> weights(static_cast<std::size_t>(dim));
  WQPatchMetric metric;

  for (int patch = 0; patch < num_patches; patch++) {
    int patch_el_offset = 0;
    int patch_ne = 0;
    GetPatchInfoFromFESpace(fes, dim, patch, patch_el_offset, patch_ne, kvs);

    const auto t0 = Clock::now();
    for (int d = 0; d < dim; d++) {
      BuildWQRule1D(*kvs[d], opts, rules[static_cast<std::size_t>(d)]);
    }
    const auto t1 = Clock::now();
    for (int d = 0; d < dim; d++) {
      BuildWQWeights1D(*kvs[d], rules[static_cast<std::size_t>(d)],
                       weights[static_cast<std::size_t>(d)]);
    }
    const auto t2 = Clock::now();
    EvaluatePatchMetric(*mesh, patch, rules, coeff, metric);
    const auto t3 = Clock::now();
    AssemblePatchDiffusionWQ(fes, patch, rules, weights, metric, opts, A);
    const auto t4 = Clock::now();

    if (timings != nullptr) {
      timings->rules_1d += Seconds(t1 - t0).count();
      timings->weights_1d += Seconds(t2 - t1).count();
      timings->metric_eval += Seconds(t3 - t2).count();
      timings->patch_assemble += Seconds(t4 - t3).count();
    }
  }

  if (!A.Finalized()) {
    const auto t_fin0 = Clock::now();
    A.Finalize();
    const auto t_fin1 = Clock::now();
    if (timings != nullptr) {
      timings->patch_assemble += Seconds(t_fin1 - t_fin0).count();
    }
  }

  if (opts.symmetrize) {
    const auto t_sym0 = Clock::now();
    std::unique_ptr<SparseMatrix> At(Transpose(A));
    std::unique_ptr<SparseMatrix> Asym(
        Add(static_cast<real_t>(0.5), A, static_cast<real_t>(0.5), *At));
    A.Swap(*Asym);
    const auto t_sym1 = Clock::now();
    if (timings != nullptr) {
      timings->symmetrize += Seconds(t_sym1 - t_sym0).count();
    }
  }

  // Single O(nnz) pass for essential Dirichlet boundary condition elimination.
  const auto t_bc0 = Clock::now();

  if (&X != &x) {
    X = x;
  }
  if (&B != &b) {
    B = b;
  }

  if (ess_tdof_list.Size() > 0) {
    Array<char> is_ess(num_tdofs);
    is_ess = 0;
    for (int k = 0; k < ess_tdof_list.Size(); k++) {
      const int tdof = ess_tdof_list[k];
      MFEM_VERIFY(tdof >= 0 && tdof < num_tdofs,
                  "Essential tdof index out of range in "
                  "AssembleDiffusionWQ.");
      is_ess[tdof] = 1;
    }

    MFEM_VERIFY(A.Finalized(),
                "Assembled SparseMatrix must be finalized before BC "
                "elimination.");

    const int *I = A.HostReadI();
    const int *J = A.HostReadJ();
    real_t *Data = A.HostReadWriteData();
    const real_t *X_data = X.HostRead();
    real_t *B_data = B.HostReadWrite();

    for (int i = 0; i < num_tdofs; i++) {
      const int row_begin = I[i];
      const int row_end = I[i + 1];

      if (is_ess[i]) {
        bool found_diag = false;
        for (int k = row_begin; k < row_end; k++) {
          if (J[k] == i) {
            Data[k] = static_cast<real_t>(1.0);
            found_diag = true;
          } else {
            Data[k] = static_cast<real_t>(0.0);
          }
        }
        MFEM_VERIFY(found_diag,
                    "Missing diagonal entry on essential boundary row in "
                    "AssembleDiffusionWQ.");
        B_data[i] = X_data[i];
      } else {
        real_t rhs_shift = static_cast<real_t>(0.0);
        for (int k = row_begin; k < row_end; k++) {
          const int j = J[k];
          if (is_ess[j]) {
            rhs_shift += Data[k] * X_data[j];
            Data[k] = static_cast<real_t>(0.0);
          }
        }
        B_data[i] -= rhs_shift;
      }
    }
  }

  const auto t_bc1 = Clock::now();
  if (timings != nullptr) {
    timings->bc_elim += Seconds(t_bc1 - t_bc0).count();
  }
}

} // namespace mfem
