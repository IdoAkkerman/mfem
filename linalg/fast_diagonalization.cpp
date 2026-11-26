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

#include "fast_diagonalization.hpp"
#include "../general/error.hpp"
#include "kronecker.hpp"

#include <cstddef>
#include <limits>

namespace mfem {

namespace {

int ComputeTotalSize(const std::vector<DenseMatrix> &U_1d,
                     const std::vector<Vector> &lambda_1d) {
  const std::size_t dim = U_1d.size();
  MFEM_VERIFY(dim >= 1 && dim == lambda_1d.size(),
              "FastDiagonalizationSolver requires non-empty U_1d and "
              "lambda_1d of equal size.");
  MFEM_VERIFY(dim <= 4,
              "FastDiagonalizationSolver supports dimensions 1 to 4.");

  const std::size_t max_int =
      static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t total = 1;
  for (std::size_t d = 0; d < dim; d++) {
    const int n_d = U_1d[d].Height();
    MFEM_VERIFY(n_d > 0 && n_d == U_1d[d].Width(),
                "Each U_1d[d] must be a non-empty square DenseMatrix.");
    MFEM_VERIFY(n_d == lambda_1d[d].Size(),
                "Each U_1d[d].Height() must match lambda_1d[d].Size().");
    MFEM_VERIFY(total <= max_int / static_cast<std::size_t>(n_d),
                "Total FastDiagonalizationSolver size overflows int.");
    total *= static_cast<std::size_t>(n_d);
  }
  return static_cast<int>(total);
}

} // namespace

FastDiagonalizationSolver::FastDiagonalizationSolver(
    const std::vector<DenseMatrix> &U_1d, const std::vector<Vector> &lambda_1d)
    : Solver(ComputeTotalSize(U_1d, lambda_1d), false), inv_eig_sum(height),
      work(height) {
  const std::size_t dim = U_1d.size();
  switch (dim) {
  case 1:
    kron_U = std::make_unique<KroneckerProductOperator<1>>(U_1d[0]);
    break;
  case 2:
    kron_U = std::make_unique<KroneckerProductOperator<2>>(U_1d[0], U_1d[1]);
    break;
  case 3:
    kron_U = std::make_unique<KroneckerProductOperator<3>>(U_1d[0], U_1d[1],
                                                           U_1d[2]);
    break;
  case 4:
    kron_U = std::make_unique<KroneckerProductOperator<4>>(U_1d[0], U_1d[1],
                                                           U_1d[2], U_1d[3]);
    break;
  default:
    MFEM_ABORT("FastDiagonalizationSolver supports dimensions 1 to 4.");
    break;
  }

  real_t *inv_data = inv_eig_sum.HostWrite();
  int stride = lambda_1d[0].Size();
  const real_t *lam0 = lambda_1d[0].HostRead();
  for (int i = 0; i < stride; i++) {
    inv_data[i] = lam0[i];
  }

  for (std::size_t d = 1; d < dim; d++) {
    const int n_d = lambda_1d[d].Size();
    const real_t *lam_d = lambda_1d[d].HostRead();
    for (int i_d = n_d - 1; i_d >= 0; i_d--) {
      const real_t lam = lam_d[i_d];
      const int offset = i_d * stride;
      for (int k = 0; k < stride; k++) {
        inv_data[offset + k] = inv_data[k] + lam;
      }
    }
    stride *= n_d;
  }

  for (int i = 0; i < height; i++) {
    const real_t denom = inv_data[i];
    MFEM_VERIFY(denom > real_t(0.0),
                "FastDiagonalizationSolver requires strictly positive "
                "eigenvalue sums.");
    inv_data[i] = real_t(1.0) / denom;
  }
}

void FastDiagonalizationSolver::SetOperator(const Operator &op) {
  MFEM_VERIFY(op.Height() == height && op.Width() == width,
              "Operator dimensions do not match FastDiagonalizationSolver.");
}

void FastDiagonalizationSolver::Mult(const Vector &x, Vector &y) const {
  MFEM_VERIFY(x.Size() == width && y.Size() == height,
              "Incompatible vector sizes in "
              "FastDiagonalizationSolver::Mult.");

  kron_U->MultTranspose(x, work);

  real_t *w_data = work.HostReadWrite();
  const real_t *inv_data = inv_eig_sum.HostRead();
  for (int i = 0; i < height; i++) {
    w_data[i] *= inv_data[i];
  }

  kron_U->Mult(work, y);
}

} // namespace mfem
