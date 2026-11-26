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

#ifndef MFEM_FAST_DIAGONALIZATION
#define MFEM_FAST_DIAGONALIZATION

#include "densemat.hpp"
#include "operator.hpp"

#include <memory>
#include <vector>

namespace mfem {

/// Fast Diagonalization solver for Kronecker-sum operators.
class FastDiagonalizationSolver : public Solver {
protected:
  std::unique_ptr<Operator> kron_U;
  Vector inv_eig_sum;
  mutable Vector work;

public:
  /// Construct from 1D eigenvector matrices @a U_1d and eigenvalue vectors
  /// @a lambda_1d (direction 0 varies fastest in memory).
  FastDiagonalizationSolver(const std::vector<DenseMatrix> &U_1d,
                            const std::vector<Vector> &lambda_1d);

  /// Verify that @a op matches the solver dimensions.
  void SetOperator(const Operator &op) override;

  /// Apply the Fast Diagonalization inverse action y = A^{-1} x.
  void Mult(const Vector &x, Vector &y) const override;
};

} // namespace mfem

#endif // MFEM_FAST_DIAGONALIZATION
