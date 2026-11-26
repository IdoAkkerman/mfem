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

#ifndef MFEM_DENSEMAT_GEN_EIGENSYSTEM
#define MFEM_DENSEMAT_GEN_EIGENSYSTEM

#include "../config/config.hpp"
#include "densemat.hpp"
#include "vector.hpp"

namespace mfem {

/// Generalized eigensystem solver for symmetric-definite dense matrices.
///
/// Solves K U = M U \Lambda with M-orthonormal eigenvectors U^T M U = I and
/// U^T K U = \Lambda, where K is symmetric and M is symmetric positive
/// definite.
class DenseMatrixGenEigensystem {
protected:
  DenseMatrix K_copy;
  DenseMatrix M_copy;
  Vector eigenvalues;
  DenseMatrix eigenvectors;
  bool evaluated;

public:
  /// Construct the generalized eigensystem for symmetric K and SPD M.
  DenseMatrixGenEigensystem(const DenseMatrix &K, const DenseMatrix &M);

  /// Compute the eigenvalues (in ascending order) and M-orthonormal
  /// eigenvectors using LAPACK.
  void Eval();

  /// Return the computed eigenvalues in ascending order.
  const Vector &Eigenvalues() const;

  /// Return the computed M-orthonormal eigenvectors stored as columns.
  const DenseMatrix &Eigenvectors() const;
};

} // namespace mfem

#endif
