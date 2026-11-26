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

#ifndef MFEM_KRONECKER_PRECONDITIONER
#define MFEM_KRONECKER_PRECONDITIONER

#include "../config/config.hpp"
#include "../general/array.hpp"
#include "fast_diagonalization.hpp"
#include "operator.hpp"

#include <vector>

namespace mfem {

/// Preconditioner bridging a full true-DoF vector (with essential Dirichlet
/// DoFs kept as identity rows) to a free-DoF FastDiagonalizationSolver.
class KroneckerPreconditioner : public Solver {
protected:
  const FastDiagonalizationSolver &fd_solver;
  Array<int> ess_tdof_list;
  Array<int> free_map;
  mutable Vector x_free;
  mutable Vector y_free;

public:
  /// Construct a KroneckerPreconditioner over a tensor-product grid.
  KroneckerPreconditioner(int full_size, const Array<int> &ess_tdof_list,
                          const Array<int> &full_sizes_1d,
                          const std::vector<Array<int>> &free_dofs_1d,
                          const Array<int> &tensor_to_tdof,
                          const FastDiagonalizationSolver &fd_solver);

  /// Verify that @a op matches the preconditioner dimensions.
  void SetOperator(const Operator &op) override;

  /// Apply the preconditioner: solve on free DoFs, identity on essential DoFs.
  void Mult(const Vector &x, Vector &y) const override;
};

} // namespace mfem

#endif
