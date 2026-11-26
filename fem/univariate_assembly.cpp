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

#include "univariate_assembly.hpp"

#include "../general/error.hpp"
#include "../linalg/sparsemat.hpp"
#include "../linalg/vector.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/nurbs.hpp"
#include "bilinearform.hpp"
#include "bilininteg.hpp"
#include "fe/fe_nurbs.hpp"
#include "fe_coll.hpp"
#include "gridfunc.hpp"

#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace mfem {

namespace {

void BuildTensorToDofMap(const FiniteElementSpace &fespace,
                         const Array<const KnotVector *> &kvs,
                         const Array<int> &full_sizes_1d,
                         Array<int> &tensor_to_dof) {
  const int dim = full_sizes_1d.Size();
  Array<int> orders(dim);
  Array<int> ne_1d(dim);
  Array<int> strides(dim);
  std::vector<Array<int>> span_idx(dim);

  int ne_total = 1;
  int n_loc = 1;
  int n_total = 1;

  for (int d = 0; d < dim; d++) {
    const KnotVector *kv = kvs[d];
    const int p_d = kv->GetOrder();
    const int n_d = full_sizes_1d[d];

    orders[d] = p_d;
    strides[d] = n_total;
    n_total *= n_d;
    n_loc *= (p_d + 1);

    span_idx[d].Reserve(kv->GetNE());
    for (int s = p_d; s < n_d; s++) {
      if ((*kv)[s] != (*kv)[s + 1]) {
        span_idx[d].Append(s);
      }
    }
    MFEM_VERIFY(span_idx[d].Size() == kv->GetNE(),
                "KnotVector non-zero span count does not match GetNE().");
    ne_1d[d] = span_idx[d].Size();
    ne_total *= ne_1d[d];
  }

  MFEM_VERIFY(ne_total == fespace.GetNE(),
              "Single-patch knot-span element count does not match mesh.");

  tensor_to_dof.SetSize(n_total);
  tensor_to_dof = -1;

  Array<int> e_idx(dim);
  Array<int> elem_dofs;

  for (int e = 0; e < ne_total; e++) {
    int rem = e;
    for (int d = 0; d < dim; d++) {
      e_idx[d] = rem % ne_1d[d];
      rem /= ne_1d[d];
    }

    fespace.GetElementDofs(e, elem_dofs);
    MFEM_VERIFY(elem_dofs.Size() == n_loc,
                "Unexpected element DoF count on NURBS element.");

    for (int loc = 0; loc < n_loc; loc++) {
      int lrem = loc;
      int tidx = 0;
      for (int d = 0; d < dim; d++) {
        const int l_d = lrem % (orders[d] + 1);
        lrem /= (orders[d] + 1);
        const int cp_d = (span_idx[d][e_idx[d]] - orders[d]) + l_d;
        MFEM_VERIFY(cp_d >= 0 && cp_d < full_sizes_1d[d],
                    "Control-point index out of range.");
        tidx += cp_d * strides[d];
      }

      const int dof = elem_dofs[loc];
      MFEM_VERIFY(dof >= 0 && dof < n_total,
                  "Element DoF index out of expected range.");
      if (tensor_to_dof[tidx] == -1) {
        tensor_to_dof[tidx] = dof;
      } else {
        MFEM_VERIFY(tensor_to_dof[tidx] == dof,
                    "Inconsistent tensor-product DoF mapping.");
      }
    }
  }

  Array<bool> seen(n_total);
  seen = false;
  for (int i = 0; i < n_total; i++) {
    const int dof = tensor_to_dof[i];
    MFEM_VERIFY(dof >= 0 && dof < n_total && !seen[dof],
                "Tensor-product DoF mapping is not a bijection.");
    seen[dof] = true;
  }
}

void ExtractPrincipalSubmatrix(const SparseMatrix &A,
                               const Array<int> &cp_to_dof,
                               const Array<int> &free_dofs, DenseMatrix &sub) {
  const int n_free = free_dofs.Size();
  sub.SetSize(n_free, n_free);
  sub = 0.0;

  if (n_free == 0) {
    return;
  }

  Array<int> dof_to_free(A.Height());
  dof_to_free = -1;
  for (int a = 0; a < n_free; a++) {
    const int dof = cp_to_dof[free_dofs[a]];
    MFEM_VERIFY(dof >= 0 && dof < A.Height(),
                "1D DoF index out of range during submatrix extraction.");
    dof_to_free[dof] = a;
  }

  const int *I = A.HostReadI();
  const int *J = A.HostReadJ();
  const real_t *data = A.HostReadData();

  for (int a = 0; a < n_free; a++) {
    const int row = cp_to_dof[free_dofs[a]];
    for (int k = I[row]; k < I[row + 1]; k++) {
      const int b = dof_to_free[J[k]];
      if (b >= 0) {
        sub(a, b) += data[k];
      }
    }
  }
}

void AssembleUnivariateDirection(const KnotVector *kv,
                                 const Array<int> &free_dofs,
                                 DenseMatrix &K_out, DenseMatrix &M_out) {
  const int p = kv->GetOrder();
  const int n = kv->GetNCP();
  MFEM_VERIFY(p >= 1 && n >= p + 1 && kv->Size() == n + p + 1,
              "Invalid 1D KnotVector for univariate assembly.");

  // Greville abscissae give exact identity map x(u) = u for any B-spline.
  Vector greville(n);
  for (int i = 0; i < n; i++) {
    real_t sum = 0.0;
    for (int j = 1; j <= p; j++) {
      sum += (*kv)[i + j];
    }
    greville(i) = sum / static_cast<real_t>(p);
  }

  std::stringstream ss;
  ss << std::setprecision(std::numeric_limits<real_t>::max_digits10);
  ss << "MFEM NURBS mesh v1.0\n\n"
     << "dimension\n1\n\n"
     << "elements\n1\n1 1 0 1\n\n"
     << "boundary\n2\n1 0 0\n2 0 1\n\n"
     << "edges\n1\n0 0 1\n\n"
     << "vertices\n2\n\n"
     << "knotvectors\n1\n";
  kv->Print(ss);
  ss << "\nweights\n";
  for (int i = 0; i < n; i++) {
    ss << "1\n";
  }
  ss << "\nFiniteElementSpace\n"
     << "FiniteElementCollection: NURBS\n"
     << "VDim: 1\n"
     << "Ordering: 1\n\n"
     << greville(0) << '\n'
     << greville(n - 1) << '\n';
  for (int i = 1; i < n - 1; i++) {
    ss << greville(i) << '\n';
  }

  Mesh mesh_1d(ss, 1, 1);
  NURBSFECollection fec_1d;
  FiniteElementSpace fespace_1d(&mesh_1d, &fec_1d);

  Array<const KnotVector *> kvs_1d(1);
  kvs_1d[0] = kv;
  Array<int> sizes_1d(1);
  sizes_1d[0] = n;
  Array<int> cp_to_dof_1d;
  BuildTensorToDofMap(fespace_1d, kvs_1d, sizes_1d, cp_to_dof_1d);

  GridFunction *nodes_1d = mesh_1d.GetNodes();
  MFEM_VERIFY(nodes_1d != nullptr && nodes_1d->Size() == n,
              "Failed to initialize 1D NURBS mesh nodes.");
  for (int i = 0; i < n; i++) {
    (*nodes_1d)(cp_to_dof_1d[i]) = greville(i);
  }

  BilinearForm k_form(&fespace_1d);
  k_form.AddDomainIntegrator(new DiffusionIntegrator);
  k_form.Assemble();
  k_form.Finalize();

  BilinearForm m_form(&fespace_1d);
  m_form.AddDomainIntegrator(new MassIntegrator);
  m_form.Assemble();
  m_form.Finalize();

  ExtractPrincipalSubmatrix(k_form.SpMat(), cp_to_dof_1d, free_dofs, K_out);
  ExtractPrincipalSubmatrix(m_form.SpMat(), cp_to_dof_1d, free_dofs, M_out);
}

} // namespace

void AssembleUnivariateLaplace(const FiniteElementSpace &fespace,
                               const Array<int> &ess_bdr,
                               UnivariateFactors &factors) {
  const Mesh *mesh = fespace.GetMesh();
  MFEM_VERIFY(mesh != nullptr && mesh->NURBSExt() != nullptr,
              "AssembleUnivariateLaplace requires a NURBS mesh.");
  MFEM_VERIFY(mesh->NURBSExt()->GetNP() == 1,
              "AssembleUnivariateLaplace requires a single-patch NURBS mesh.");

  const NURBSExtension *nurbs_ext = fespace.GetNURBSext();
  MFEM_VERIFY(nurbs_ext != nullptr,
              "FiniteElementSpace must use a NURBSFECollection.");
  MFEM_VERIFY(nurbs_ext->GetNP() == 1,
              "NURBSExtension must represent a single tensor-product patch.");
  MFEM_VERIFY(fespace.GetVDim() == 1,
              "AssembleUnivariateLaplace supports scalar spaces (vdim == 1).");
  MFEM_VERIFY(fespace.GetVSize() == fespace.GetTrueVSize(),
              "AssembleUnivariateLaplace requires conforming scalar DoFs.");

  const int dim = mesh->Dimension();
  MFEM_VERIFY(dim >= 1 && dim <= 3, "Mesh dimension must be 1, 2, or 3.");
  MFEM_VERIFY(mesh->GetNE() > 0,
              "NURBS mesh must contain at least one element.");

  const FiniteElement *fe = fespace.GetFE(0);
  MFEM_VERIFY(fe != nullptr,
              "FiniteElementSpace returned a null finite element.");
  const NURBSFiniteElement *nfe = dynamic_cast<const NURBSFiniteElement *>(fe);
  MFEM_VERIFY(nfe != nullptr,
              "FiniteElementSpace must provide NURBS finite elements.");
  nfe->Reset();
  nurbs_ext->LoadFE(0, fe);

  const Array<const KnotVector *> &fe_kvs = nfe->KnotVectors();
  MFEM_VERIFY(fe_kvs.Size() == dim,
              "NURBSFiniteElement knot vector count does not match dim.");

  Array<const KnotVector *> kvs(dim);
  factors.full_sizes_1d.SetSize(dim);
  Array<int> strides(dim);
  std::int64_t prod_n = 1;

  for (int d = 0; d < dim; d++) {
    kvs[d] = fe_kvs[d];
    MFEM_VERIFY(kvs[d] != nullptr, "Null KnotVector in NURBSExtension.");
    const int n_d = kvs[d]->GetNCP();
    const int p_d = kvs[d]->GetOrder();
    MFEM_VERIFY(p_d >= 1 && n_d >= p_d + 1,
                "Invalid 1D KnotVector order or control-point count.");
    factors.full_sizes_1d[d] = n_d;
    strides[d] = static_cast<int>(prod_n);
    prod_n *= static_cast<std::int64_t>(n_d);
  }

  const int n_total = fespace.GetTrueVSize();
  MFEM_VERIFY(prod_n == static_cast<std::int64_t>(n_total),
              "Product of 1D control-point counts does not equal TrueVSize.");

  BuildTensorToDofMap(fespace, kvs, factors.full_sizes_1d,
                      factors.tensor_to_tdof);

  Array<int> ess_tdof_list;
  fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  Array<bool> is_essential(n_total);
  is_essential = false;
  for (int i = 0; i < ess_tdof_list.Size(); i++) {
    const int tdof = ess_tdof_list[i];
    MFEM_VERIFY(tdof >= 0 && tdof < n_total,
                "Essential true DoF index out of range.");
    is_essential[tdof] = true;
  }

  Array<bool> left_ess(dim);
  Array<bool> right_ess(dim);
  left_ess = true;
  right_ess = true;

  for (int tidx = 0; tidx < n_total; tidx++) {
    if (!is_essential[factors.tensor_to_tdof[tidx]]) {
      for (int d = 0; d < dim; d++) {
        const int i_d = (tidx / strides[d]) % factors.full_sizes_1d[d];
        if (i_d == 0) {
          left_ess[d] = false;
        }
        if (i_d == factors.full_sizes_1d[d] - 1) {
          right_ess[d] = false;
        }
      }
    }
  }

  for (int tidx = 0; tidx < n_total; tidx++) {
    bool in_tp_free = true;
    for (int d = 0; d < dim; d++) {
      const int i_d = (tidx / strides[d]) % factors.full_sizes_1d[d];
      if ((left_ess[d] && i_d == 0) ||
          (right_ess[d] && i_d == factors.full_sizes_1d[d] - 1)) {
        in_tp_free = false;
        break;
      }
    }
    const bool is_true_free = !is_essential[factors.tensor_to_tdof[tidx]];
    MFEM_VERIFY(in_tp_free == is_true_free,
                "Essential boundary conditions do not form a "
                "tensor-product free-DoF set.");
  }

  factors.free_dofs_1d.clear();
  factors.free_dofs_1d.resize(dim);
  factors.K.clear();
  factors.K.resize(dim);
  factors.M.clear();
  factors.M.resize(dim);

  for (int d = 0; d < dim; d++) {
    const int n_d = factors.full_sizes_1d[d];
    const int start = left_ess[d] ? 1 : 0;
    const int end = right_ess[d] ? (n_d - 1) : n_d;

    factors.free_dofs_1d[d].SetSize(0);
    if (start < end) {
      factors.free_dofs_1d[d].Reserve(end - start);
      for (int i = start; i < end; i++) {
        factors.free_dofs_1d[d].Append(i);
      }
    }

    AssembleUnivariateDirection(kvs[d], factors.free_dofs_1d[d], factors.K[d],
                                factors.M[d]);
  }
}

} // namespace mfem
