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

#ifndef MFEM_UNIVARIATE_ASSEMBLY
#define MFEM_UNIVARIATE_ASSEMBLY

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "../linalg/densemat.hpp"
#include "fespace.hpp"

#include <vector>

namespace mfem {

/// Univariate 1D B-spline stiffness and mass factors for a single-patch
/// tensor-product NURBS space with tensor-product essential BCs.
struct UnivariateFactors {
  /// 1D SPD stiffness matrices (free DoFs).
  std::vector<DenseMatrix> K;
  /// 1D SPD mass matrices (free DoFs).
  std::vector<DenseMatrix> M;
  /// Full 1D DoF counts (n_0, ..., n_{d-1}).
  Array<int> full_sizes_1d;
  /// Free 1D DoF indices per direction.
  std::vector<Array<int>> free_dofs_1d;
  /// Maps flat i_0 + n_0*(i_1 + ...) -> MFEM tdof
  Array<int> tensor_to_tdof;
};

/// Assemble 1D parametric B-spline stiffness and mass matrices on the free
/// DoFs of each direction for a single-patch tensor-product NURBS space.
void AssembleUnivariateLaplace(const FiniteElementSpace &fespace,
                               const Array<int> &ess_bdr,
                               UnivariateFactors &factors);

} // namespace mfem

#endif
