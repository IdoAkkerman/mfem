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

#include "densemat_gen_eigensystem.hpp"
#include "../general/error.hpp"

#ifdef MFEM_USE_LAPACK
extern "C" void dsygv_(int *itype, char *jobz, char *uplo, int *n, double *a,
                       int *lda, double *b, int *ldb, double *w, double *work,
                       int *lwork, int *info);
extern "C" void ssygv_(int *itype, char *jobz, char *uplo, int *n, float *a,
                       int *lda, float *b, int *ldb, float *w, float *work,
                       int *lwork, int *info);
#endif

namespace mfem {

DenseMatrixGenEigensystem::DenseMatrixGenEigensystem(const DenseMatrix &K,
                                                     const DenseMatrix &M)
    : evaluated(false) {
  MFEM_VERIFY(K.Height() == K.Width(),
              "DenseMatrixGenEigensystem requires a square matrix K.");
  MFEM_VERIFY(M.Height() == M.Width(),
              "DenseMatrixGenEigensystem requires a square matrix M.");
  MFEM_VERIFY(K.Height() > 0,
              "DenseMatrixGenEigensystem requires non-empty matrices.");
  MFEM_VERIFY(K.Height() == M.Height(),
              "DenseMatrixGenEigensystem requires K and M to have "
              "identical dimensions.");

  K_copy = K;
  M_copy = M;
}

void DenseMatrixGenEigensystem::Eval() {
#ifdef MFEM_USE_LAPACK
  evaluated = false;

  int itype = 1;
  char jobz = 'V';
  char uplo = 'U';
  int n = K_copy.Height();
  int lda = n;
  int ldb = n;
  int info = 0;

  eigenvectors = K_copy;
  DenseMatrix m_work(M_copy);
  eigenvalues.SetSize(n);
  real_t *w_data = eigenvalues.HostWrite();

  int lwork = -1;
  real_t work_query = 0.0;

  auto call_sygv = [&](real_t *work_ptr) {
#ifdef MFEM_USE_SINGLE
    ssygv_(&itype, &jobz, &uplo, &n, eigenvectors.Data(), &lda, m_work.Data(),
           &ldb, w_data, work_ptr, &lwork, &info);
#else
    dsygv_(&itype, &jobz, &uplo, &n, eigenvectors.Data(), &lda, m_work.Data(),
           &ldb, w_data, work_ptr, &lwork, &info);
#endif
  };

  call_sygv(&work_query);
  MFEM_VERIFY(info == 0,
              "DenseMatrixGenEigensystem::Eval(): LAPACK sygv workspace "
              "query failed with info = "
                  << info << ".");

  lwork = static_cast<int>(work_query);
  if (lwork < 1) {
    lwork = 1;
  }

  Vector work(lwork);
  call_sygv(work.HostWrite());
  MFEM_VERIFY(info == 0,
              "DenseMatrixGenEigensystem::Eval(): LAPACK sygv failed with "
              "info = "
                  << info << ".");

  evaluated = true;
#else
  MFEM_ABORT("DenseMatrixGenEigensystem requires MFEM_USE_LAPACK.");
#endif
}

const Vector &DenseMatrixGenEigensystem::Eigenvalues() const {
  MFEM_VERIFY(evaluated,
              "DenseMatrixGenEigensystem::Eval() must be called before "
              "Eigenvalues().");
  return eigenvalues;
}

const DenseMatrix &DenseMatrixGenEigensystem::Eigenvectors() const {
  MFEM_VERIFY(evaluated,
              "DenseMatrixGenEigensystem::Eval() must be called before "
              "Eigenvectors().");
  return eigenvectors;
}

} // namespace mfem
