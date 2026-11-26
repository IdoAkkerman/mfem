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

#include "kronecker_preconditioner.hpp"
#include "../general/error.hpp"

namespace mfem {

KroneckerPreconditioner::KroneckerPreconditioner(
    int full_size, const Array<int> &ess_tdof_list_,
    const Array<int> &full_sizes_1d,
    const std::vector<Array<int>> &free_dofs_1d,
    const Array<int> &tensor_to_tdof,
    const FastDiagonalizationSolver &fd_solver_)
    : Solver(full_size), fd_solver(fd_solver_), ess_tdof_list(ess_tdof_list_) {
  const int dim = full_sizes_1d.Size();
  MFEM_VERIFY(dim >= 1 && dim == static_cast<int>(free_dofs_1d.size()),
              "Dimension mismatch between full_sizes_1d and free_dofs_1d.");
  MFEM_VERIFY(full_size >= 0, "full_size must be non-negative.");
  MFEM_VERIFY(tensor_to_tdof.Size() == full_size,
              "tensor_to_tdof.Size() must equal full_size.");

  long long full_prod = 1;
  long long free_prod = 1;
  for (int d = 0; d < dim; d++) {
    const int n_full_d = full_sizes_1d[d];
    const int n_free_d = free_dofs_1d[d].Size();
    MFEM_VERIFY(n_full_d > 0, "1D full sizes must be positive.");
    MFEM_VERIFY(n_free_d >= 0 && n_free_d <= n_full_d,
                "1D free DoF count out of valid range.");

    for (int j = 0; j < n_free_d; j++) {
      const int idx = free_dofs_1d[d][j];
      MFEM_VERIFY(idx >= 0 && idx < n_full_d,
                  "1D free DoF index out of bounds.");
    }

    full_prod *= n_full_d;
    free_prod *= n_free_d;
  }

  MFEM_VERIFY(full_prod == static_cast<long long>(full_size),
              "Product of full_sizes_1d does not equal full_size.");
  MFEM_VERIFY(fd_solver.Width() == fd_solver.Height(),
              "fd_solver must be square.");
  MFEM_VERIFY(free_prod == static_cast<long long>(fd_solver.Width()),
              "Product of free_dofs_1d sizes does not equal "
              "fd_solver.Width().");
  MFEM_VERIFY(fd_solver.Width() + ess_tdof_list.Size() == full_size,
              "fd_solver.Width() + ess_tdof_list.Size() must equal "
              "full_size.");

  for (int i = 0; i < ess_tdof_list.Size(); i++) {
    const int idx = ess_tdof_list[i];
    MFEM_VERIFY(idx >= 0 && idx < full_size,
                "Essential DoF index out of bounds.");
  }

  const int n_free = static_cast<int>(free_prod);
  free_map.SetSize(n_free);
  x_free.SetSize(n_free);
  y_free.SetSize(n_free);

  if (n_free > 0) {
    Array<int> strides(dim);
    strides[0] = 1;
    for (int d = 1; d < dim; d++) {
      strides[d] = strides[d - 1] * full_sizes_1d[d - 1];
    }

    Array<int> multi_idx(dim);
    multi_idx = 0;

    for (int k = 0; k < n_free; k++) {
      int flat_idx = 0;
      for (int d = 0; d < dim; d++) {
        flat_idx += free_dofs_1d[d][multi_idx[d]] * strides[d];
      }

      const int tdof = tensor_to_tdof[flat_idx];
      MFEM_VERIFY(tdof >= 0 && tdof < full_size,
                  "Mapped free DoF index out of bounds.");
      free_map[k] = tdof;

      for (int d = 0; d < dim; d++) {
        multi_idx[d]++;
        if (multi_idx[d] < free_dofs_1d[d].Size()) {
          break;
        }
        multi_idx[d] = 0;
      }
    }
  }

  Array<int> seen(full_size);
  seen = 0;
  for (int k = 0; k < n_free; k++) {
    const int idx = free_map[k];
    MFEM_VERIFY(seen[idx] == 0,
                "Duplicate index in KroneckerPreconditioner free_map.");
    seen[idx] = 1;
  }
  for (int i = 0; i < ess_tdof_list.Size(); i++) {
    const int idx = ess_tdof_list[i];
    MFEM_VERIFY(seen[idx] == 0,
                "ess_tdof_list overlaps free_map or contains duplicates.");
    seen[idx] = 1;
  }
}

void KroneckerPreconditioner::SetOperator(const Operator &op) {
  MFEM_VERIFY(op.Height() == height && op.Width() == width,
              "Operator dimensions do not match KroneckerPreconditioner.");
}

void KroneckerPreconditioner::Mult(const Vector &x, Vector &y) const {
  MFEM_VERIFY(x.Size() == width && y.Size() == height,
              "Incompatible vector sizes in KroneckerPreconditioner::Mult.");

  const real_t *x_data = x.HostRead();
  const int *fmap = free_map.HostRead();
  const int n_free = free_map.Size();

  real_t *xf_data = x_free.HostWrite();
  for (int k = 0; k < n_free; k++) {
    xf_data[k] = x_data[fmap[k]];
  }

  if (n_free > 0) {
    fd_solver.Mult(x_free, y_free);
  }

  const real_t *yf_data = y_free.HostRead();
  real_t *y_data = (&x == &y) ? y.HostReadWrite() : y.HostWrite();
  for (int k = 0; k < n_free; k++) {
    y_data[fmap[k]] = yf_data[k];
  }

  const int *ess = ess_tdof_list.HostRead();
  const int n_ess = ess_tdof_list.Size();
  for (int i = 0; i < n_ess; i++) {
    const int idx = ess[i];
    y_data[idx] = x_data[idx];
  }
}

} // namespace mfem
