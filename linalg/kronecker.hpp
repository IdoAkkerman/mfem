#ifndef MFEM_KRONECKER
#define MFEM_KRONECKER

#include "../config/config.hpp"
#include "../general/error.hpp"
#include "densemat.hpp"
#include "operator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#define MFEM_KRON_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define MFEM_KRON_RESTRICT __restrict
#else
#define MFEM_KRON_RESTRICT
#endif

#ifdef MFEM_USE_LAPACK
extern "C" {
#ifdef MFEM_USE_SINGLE
void sgemm_(char *, char *, int *, int *, int *, float *, float *, int *,
            float *, int *, float *, float *, int *);
#else
void dgemm_(char *, char *, int *, int *, int *, double *, double *, int *,
            double *, int *, double *, double *, int *);
#endif
}
#endif

namespace mfem {

/// Matrix-free Kronecker product operator (D_{N-1} \otimes ... \otimes D_0).
template <size_t N> class KroneckerProductOperator : public Operator {
  static_assert(N > 0,
                "KroneckerProductOperator requires at least one matrix.");

public:
  /// Inner contraction backend selection.
  enum class ContractionBackend {
    Auto,      ///< Hybrid: BLAS for larger GEMMs, BlockedCPP for small.
    BLAS,      ///< Force BLAS GEMM when MFEM_USE_LAPACK is enabled.
    BlockedCPP ///< Pure C++ L1-tiled, register-blocked SIMD kernel.
  };

protected:
  std::array<DenseMatrix, N> matrices;
  std::array<DenseMatrix, N> matrices_t;

public:
  /// Passing (D_0, ..., D_{N-1}) computes (D_{N-1} \otimes ... \otimes D_0) x.
  template <typename... T>
  KroneckerProductOperator(const T &...mats)
      : Operator((mats.Height() * ...), (mats.Width() * ...)),
        matrices{mats...}, backend(ContractionBackend::Auto),
        use_even_odd(true) {
    static_assert(
        (std::is_same_v<T, DenseMatrix> && ...),
        "KroneckerProductOperator only supports DenseMatrix arguments.");
    static_assert(N == sizeof...(T),
                  "Operator dimension and matrix count must match.");

    for (size_t k = 0; k < N; ++k) {
      MFEM_VERIFY(matrices[k].Height() > 0 && matrices[k].Width() > 0,
                  "Kronecker factors must have positive dimensions.");
      matrices_t[k].Transpose(matrices[k]);
      eo_fwd[k] = BuildEvenOdd(matrices[k]);
      eo_adj[k] = BuildEvenOdd(matrices_t[k]);
    }

    mult_plan = BuildPlan(matrices);
    trans_plan = BuildPlan(matrices_t);

    const size_t peak = std::max(mult_plan.max_inter, trans_plan.max_inter);
    if (N >= 2 && peak > 0) {
      buf_a.SetSize(static_cast<int>(peak));
    }
    if (N >= 3 && peak > 0) {
      buf_b.SetSize(static_cast<int>(peak));
    }
  }

  /// Select contraction kernel backend (Auto, BLAS, or BlockedCPP).
  void SetBackend(ContractionBackend b) { backend = b; }

  /// Return currently active contraction backend.
  ContractionBackend GetBackend() const { return backend; }

  /// Enable or disable automatic even-odd symmetry acceleration.
  void SetUseEvenOdd(bool enable) { use_even_odd = enable; }

  /// Return whether automatic even-odd symmetry acceleration is enabled.
  bool GetUseEvenOdd() const { return use_even_odd; }

  /// Compute y = (D_{N-1} \otimes ... \otimes D_0) x.
  void Mult(const Vector &x, Vector &y) const override {
    MFEM_VERIFY(x.Size() == width,
                "Input vector size and operator width must match.");
    MFEM_VERIFY(y.Size() == height,
                "Output vector size and operator height must match.");
    ApplyPlan(mult_plan, matrices, eo_fwd, x, y);
  }

  /// Compute y = (D_{N-1}^T \otimes ... \otimes D_0^T) x.
  void MultTranspose(const Vector &x, Vector &y) const override {
    MFEM_VERIFY(x.Size() == height,
                "Input vector size and operator height must match.");
    MFEM_VERIFY(y.Size() == width,
                "Output vector size and operator width must match.");
    ApplyPlan(trans_plan, matrices_t, eo_adj, x, y);
  }

  enum class SymmetryType { None, RowEvenOdd, ColEvenOdd, Centrosymmetric };
  /// Return detected even-odd symmetry type for factor `mode`.
  SymmetryType GetSymmetryType(size_t mode, bool transpose = false) const {
    MFEM_VERIFY(mode < N, "Mode index out of range.");
    return transpose ? eo_adj[mode].type : eo_fwd[mode].type;
  }

private:
  struct ModeStep {
    size_t mode = 0;
    size_t L = 1;
    size_t R = 1;
    size_t m = 0;
    size_t n = 0;
  };

  struct ContractionPlan {
    std::array<ModeStep, N> steps;
    size_t max_inter = 0;
  };

  struct EvenOddFactor {
    SymmetryType type = SymmetryType::None;
    DenseMatrix a_even;
    DenseMatrix a_odd;
    std::vector<size_t> even_idx;
    std::vector<size_t> odd_idx;
  };

  ContractionBackend backend;
  bool use_even_odd;
  ContractionPlan mult_plan;
  ContractionPlan trans_plan;
  std::array<EvenOddFactor, N> eo_fwd;
  std::array<EvenOddFactor, N> eo_adj;

  mutable Vector buf_a;
  mutable Vector buf_b;
  mutable Vector eo_in;
  mutable Vector eo_out;

  static ContractionPlan BuildPlan(const std::array<DenseMatrix, N> &mats) {
    std::array<size_t, N> order;
    std::array<size_t, N> cur_dims;
    for (size_t k = 0; k < N; ++k) {
      order[k] = k;
      cur_dims[k] = static_cast<size_t>(mats[k].Width());
    }

    // Contract shrinking modes first to minimize FLOPs and scratch memory.
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      const double ma = static_cast<double>(mats[a].Height());
      const double na = static_cast<double>(mats[a].Width());
      const double mb = static_cast<double>(mats[b].Height());
      const double nb = static_cast<double>(mats[b].Width());
      return (ma - na) / (ma * na) < (mb - nb) / (mb * nb);
    });

    ContractionPlan plan;
    plan.max_inter = 0;

    for (size_t s = 0; s < N; ++s) {
      const size_t mode = order[s];
      size_t L = 1;
      for (size_t p = 0; p < mode; ++p) {
        L *= cur_dims[p];
      }
      size_t R = 1;
      for (size_t p = mode + 1; p < N; ++p) {
        R *= cur_dims[p];
      }

      const size_t m = static_cast<size_t>(mats[mode].Height());
      const size_t n = static_cast<size_t>(mats[mode].Width());
      plan.steps[s] = ModeStep{mode, L, R, m, n};
      cur_dims[mode] = m;

      if (s + 1 < N) {
        plan.max_inter = std::max(plan.max_inter, L * m * R);
      }
    }
    return plan;
  }

  static EvenOddFactor BuildEvenOdd(const DenseMatrix &A) {
    EvenOddFactor eo;
    const size_t m = static_cast<size_t>(A.Height());
    const size_t n = static_cast<size_t>(A.Width());
    constexpr size_t kMinDim = 8;
    if (m < kMinDim || n < kMinDim) {
      return eo;
    }

    const real_t *a_ptr = A.HostRead();
    real_t max_abs = 0.0;
    const size_t total = m * n;
    for (size_t idx = 0; idx < total; ++idx) {
      const real_t v = std::abs(a_ptr[idx]);
      if (v > max_abs) {
        max_abs = v;
      }
    }
    if (max_abs == real_t(0)) {
      return eo;
    }

    const real_t tol =
        real_t(1000.0) * std::numeric_limits<real_t>::epsilon() * max_abs;
    const size_t m_o = m / 2;
    const size_t m_e = (m + 1) / 2;
    const size_t n_o = n / 2;
    const size_t n_e = (n + 1) / 2;

    auto at = [&](size_t i, size_t j) -> real_t { return a_ptr[i + j * m]; };

    // 1. Check row even-odd symmetry (e.g., transposed symmetric eigvecs).
    {
      bool ok = true;
      std::vector<size_t> ev;
      std::vector<size_t> od;
      ev.reserve(m_e);
      od.reserve(m_o);
      for (size_t i = 0; i < m && ok; ++i) {
        bool is_even = true;
        bool is_odd = (n % 2 == 0) || (std::abs(at(i, n_o)) <= tol);
        for (size_t j = 0; j < n_o; ++j) {
          const real_t a = at(i, j);
          const real_t b = at(i, n - 1 - j);
          if (std::abs(a - b) > tol) {
            is_even = false;
          }
          if (std::abs(a + b) > tol) {
            is_odd = false;
          }
          if (!is_even && !is_odd) {
            break;
          }
        }
        if (is_even) {
          ev.push_back(i);
        } else if (is_odd) {
          od.push_back(i);
        } else {
          ok = false;
        }
      }
      if (ok) {
        eo.type = SymmetryType::RowEvenOdd;
        eo.even_idx = std::move(ev);
        eo.odd_idx = std::move(od);
        const size_t pe = eo.even_idx.size();
        const size_t po = eo.odd_idx.size();
        if (pe > 0) {
          eo.a_even.SetSize(static_cast<int>(pe), static_cast<int>(n_e));
          real_t *dst = eo.a_even.HostWrite();
          for (size_t j = 0; j < n_e; ++j) {
            for (size_t r = 0; r < pe; ++r) {
              dst[r + j * pe] = at(eo.even_idx[r], j);
            }
          }
        }
        if (po > 0) {
          eo.a_odd.SetSize(static_cast<int>(po), static_cast<int>(n_o));
          real_t *dst = eo.a_odd.HostWrite();
          for (size_t j = 0; j < n_o; ++j) {
            for (size_t r = 0; r < po; ++r) {
              dst[r + j * po] = at(eo.odd_idx[r], j);
            }
          }
        }
        return eo;
      }
    }

    // 2. Check column even-odd symmetry (e.g., symmetric eigenvectors).
    {
      bool ok = true;
      std::vector<size_t> ev;
      std::vector<size_t> od;
      ev.reserve(n_e);
      od.reserve(n_o);
      for (size_t j = 0; j < n && ok; ++j) {
        bool is_even = true;
        bool is_odd = (m % 2 == 0) || (std::abs(at(m_o, j)) <= tol);
        for (size_t i = 0; i < m_o; ++i) {
          const real_t a = at(i, j);
          const real_t b = at(m - 1 - i, j);
          if (std::abs(a - b) > tol) {
            is_even = false;
          }
          if (std::abs(a + b) > tol) {
            is_odd = false;
          }
          if (!is_even && !is_odd) {
            break;
          }
        }
        if (is_even) {
          ev.push_back(j);
        } else if (is_odd) {
          od.push_back(j);
        } else {
          ok = false;
        }
      }
      if (ok) {
        eo.type = SymmetryType::ColEvenOdd;
        eo.even_idx = std::move(ev);
        eo.odd_idx = std::move(od);
        const size_t qe = eo.even_idx.size();
        const size_t qo = eo.odd_idx.size();
        if (qe > 0) {
          eo.a_even.SetSize(static_cast<int>(m_e), static_cast<int>(qe));
          real_t *dst = eo.a_even.HostWrite();
          for (size_t c = 0; c < qe; ++c) {
            const size_t col = eo.even_idx[c];
            for (size_t i = 0; i < m_e; ++i) {
              dst[i + c * m_e] = at(i, col);
            }
          }
        }
        if (qo > 0) {
          eo.a_odd.SetSize(static_cast<int>(m_o), static_cast<int>(qo));
          real_t *dst = eo.a_odd.HostWrite();
          for (size_t c = 0; c < qo; ++c) {
            const size_t col = eo.odd_idx[c];
            for (size_t i = 0; i < m_o; ++i) {
              dst[i + c * m_o] = at(i, col);
            }
          }
        }
        return eo;
      }
    }

    // 3. Check centrosymmetry: A(m - 1 - i, n - 1 - j) == A(i, j).
    {
      bool ok = true;
      for (size_t j = 0; j < n && ok; ++j) {
        for (size_t i = 0; i < m_e; ++i) {
          const real_t a = at(i, j);
          const real_t b = at(m - 1 - i, n - 1 - j);
          if (std::abs(a - b) > tol) {
            ok = false;
            break;
          }
        }
      }
      if (ok) {
        eo.type = SymmetryType::Centrosymmetric;
        eo.a_even.SetSize(static_cast<int>(m_e), static_cast<int>(n_e));
        eo.a_odd.SetSize(static_cast<int>(m_o), static_cast<int>(n_o));
        real_t *ae = eo.a_even.HostWrite();
        real_t *ao = eo.a_odd.HostWrite();
        for (size_t j = 0; j < n_o; ++j) {
          for (size_t i = 0; i < m_e; ++i) {
            ae[i + j * m_e] = real_t(0.5) * (at(i, j) + at(i, n - 1 - j));
          }
          for (size_t i = 0; i < m_o; ++i) {
            ao[i + j * m_o] = real_t(0.5) * (at(i, j) - at(i, n - 1 - j));
          }
        }
        if (n % 2 == 1) {
          for (size_t i = 0; i < m_e; ++i) {
            ae[i + n_o * m_e] = at(i, n_o);
          }
        }
        return eo;
      }
    }

    return eo;
  }

  void ApplyPlan(const ContractionPlan &plan,
                 const std::array<DenseMatrix, N> &mats,
                 const std::array<EvenOddFactor, N> &eos, const Vector &x,
                 Vector &y) const {
    if (N == 1) {
      const real_t *in_ptr = x.HostRead();
      const real_t *y_raw = y.GetData();
      const bool aliases =
          (&x == &y) || (y_raw != nullptr && y_raw < in_ptr + x.Size() &&
                         in_ptr < y_raw + y.Size());
      if (aliases) {
        if (buf_a.Size() < x.Size()) {
          buf_a.SetSize(x.Size());
        }
        std::copy_n(in_ptr, x.Size(), buf_a.HostWrite());
        in_ptr = buf_a.HostRead();
      }
      real_t *out_ptr = y.HostWrite();
      const ModeStep &st = plan.steps[0];
      ApplyMode(st, mats[st.mode], eos[st.mode], in_ptr, out_ptr);
      return;
    }

    // Zero-copy boundary modes: x -> buf_a -> ... -> y.
    const real_t *cur_in = x.HostRead();
    real_t *cur_out = buf_a.HostWrite();
    {
      const ModeStep &st = plan.steps[0];
      ApplyMode(st, mats[st.mode], eos[st.mode], cur_in, cur_out);
    }

    bool a_is_src = true;
    for (size_t s = 1; s + 1 < N; ++s) {
      cur_in = a_is_src ? buf_a.HostRead() : buf_b.HostRead();
      cur_out = a_is_src ? buf_b.HostWrite() : buf_a.HostWrite();
      const ModeStep &st = plan.steps[s];
      ApplyMode(st, mats[st.mode], eos[st.mode], cur_in, cur_out);
      a_is_src = !a_is_src;
    }

    cur_in = a_is_src ? buf_a.HostRead() : buf_b.HostRead();
    cur_out = y.HostWrite();
    {
      const ModeStep &st = plan.steps[N - 1];
      ApplyMode(st, mats[st.mode], eos[st.mode], cur_in, cur_out);
    }
  }

  void ApplyMode(const ModeStep &st, const DenseMatrix &A,
                 const EvenOddFactor &eo, const real_t *in_ptr,
                 real_t *out_ptr) const {
    if (use_even_odd && eo.type != SymmetryType::None) {
      ApplyEvenOdd(st.L, st.R, st.m, st.n, eo, in_ptr, out_ptr);
    } else {
      ContractDense(st.L, st.R, st.m, st.n, A.HostRead(), in_ptr, out_ptr);
    }
  }

  bool ShouldUseBLAS(size_t L, size_t R, size_t m, size_t n) const {
#ifdef MFEM_USE_LAPACK
    const size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
    if (L > max_int || R > max_int || m > max_int || n > max_int) {
      return false;
    }
    if (backend == ContractionBackend::BLAS) {
      return true;
    }
    if (backend == ContractionBackend::BlockedCPP) {
      return false;
    }
    // Auto: use BLAS when a single GEMM amortizes call/packing overhead.
    const size_t gm = (L == 1) ? m : L;
    const size_t gn = (L == 1) ? R : m;
    const size_t gk = n;
    constexpr size_t kMinSide = 8;
    constexpr size_t kMinWork = 4096;
    return (gm >= kMinSide && gn >= kMinSide && gk >= kMinSide &&
            gm * gn * gk >= kMinWork);
#else
    (void)L;
    (void)R;
    (void)m;
    (void)n;
    return false;
#endif
  }

  void ContractDense(size_t L, size_t R, size_t m, size_t n,
                     const real_t *a_ptr, const real_t *x_ptr,
                     real_t *y_ptr) const {
    if (ShouldUseBLAS(L, R, m, n)) {
      ContractBLAS(L, R, m, n, a_ptr, x_ptr, y_ptr);
    } else {
      ContractBlockedCPP(L, R, m, n, a_ptr, x_ptr, y_ptr);
    }
  }

  static void CallGemm(char transa, char transb, int M, int N_col, int K,
                       const real_t *A, int lda, const real_t *B, int ldb,
                       real_t *C, int ldc) {
#ifdef MFEM_USE_LAPACK
    real_t alpha = real_t(1);
    real_t beta = real_t(0);
#ifdef MFEM_USE_SINGLE
    sgemm_(&transa, &transb, &M, &N_col, &K, &alpha, const_cast<float *>(A),
           &lda, const_cast<float *>(B), &ldb, &beta, C, &ldc);
#else
    dgemm_(&transa, &transb, &M, &N_col, &K, &alpha, const_cast<double *>(A),
           &lda, const_cast<double *>(B), &ldb, &beta, C, &ldc);
#endif
#else
    (void)transa;
    (void)transb;
    (void)M;
    (void)N_col;
    (void)K;
    (void)A;
    (void)lda;
    (void)B;
    (void)ldb;
    (void)C;
    (void)ldc;
#endif
  }

  static void ContractBLAS(size_t L, size_t R, size_t m, size_t n,
                           const real_t *a_ptr, const real_t *x_ptr,
                           real_t *y_ptr) {
    if (L == 1) {
      // Y(m x R) = A(m x n) * X(n x R)
      CallGemm('N', 'N', static_cast<int>(m), static_cast<int>(R),
               static_cast<int>(n), a_ptr, static_cast<int>(m), x_ptr,
               static_cast<int>(n), y_ptr, static_cast<int>(m));
      return;
    }

    // For each slice r: Y_r(L x m) = X_r(L x n) * A(m x n)^T
    const int Li = static_cast<int>(L);
    const int mi = static_cast<int>(m);
    const int ni = static_cast<int>(n);
    const size_t x_stride = L * n;
    const size_t y_stride = L * m;
    for (size_t r = 0; r < R; ++r) {
      CallGemm('N', 'T', Li, mi, ni, x_ptr + r * x_stride, Li, a_ptr, mi,
               y_ptr + r * y_stride, Li);
    }
  }

  static void ContractBlockedCPP(size_t L, size_t R, size_t m, size_t n,
                                 const real_t *MFEM_KRON_RESTRICT a_ptr,
                                 const real_t *MFEM_KRON_RESTRICT x_ptr,
                                 real_t *MFEM_KRON_RESTRICT y_ptr) {
    constexpr size_t kMBlk = 64;
    constexpr size_t kKBlk = 32;
    constexpr size_t kRBlk = 32;
    constexpr size_t kLBlk = 64;

    if (L == 1) {
      for (size_t r0 = 0; r0 < R; r0 += kRBlk) {
        const size_t r1 = std::min(r0 + kRBlk, R);
        for (size_t i0 = 0; i0 < m; i0 += kMBlk) {
          const size_t i1 = std::min(i0 + kMBlk, m);
          for (size_t j0 = 0; j0 < n; j0 += kKBlk) {
            const size_t j1 = std::min(j0 + kKBlk, n);
            size_t r = r0;
            for (; r + 4 <= r1; r += 4) {
              real_t *MFEM_KRON_RESTRICT y0 = y_ptr + (r + 0) * m;
              real_t *MFEM_KRON_RESTRICT y1 = y_ptr + (r + 1) * m;
              real_t *MFEM_KRON_RESTRICT y2 = y_ptr + (r + 2) * m;
              real_t *MFEM_KRON_RESTRICT y3 = y_ptr + (r + 3) * m;
              if (j0 == 0) {
                for (size_t i = i0; i < i1; ++i) {
                  y0[i] = real_t(0);
                  y1[i] = real_t(0);
                  y2[i] = real_t(0);
                  y3[i] = real_t(0);
                }
              }
              for (size_t j = j0; j < j1; ++j) {
                const real_t *MFEM_KRON_RESTRICT ac = a_ptr + j * m;
                const real_t x0 = x_ptr[j + (r + 0) * n];
                const real_t x1 = x_ptr[j + (r + 1) * n];
                const real_t x2 = x_ptr[j + (r + 2) * n];
                const real_t x3 = x_ptr[j + (r + 3) * n];
                for (size_t i = i0; i < i1; ++i) {
                  const real_t a = ac[i];
                  y0[i] += a * x0;
                  y1[i] += a * x1;
                  y2[i] += a * x2;
                  y3[i] += a * x3;
                }
              }
            }
            for (; r < r1; ++r) {
              real_t *MFEM_KRON_RESTRICT y0 = y_ptr + r * m;
              if (j0 == 0) {
                for (size_t i = i0; i < i1; ++i) {
                  y0[i] = real_t(0);
                }
              }
              for (size_t j = j0; j < j1; ++j) {
                const real_t *MFEM_KRON_RESTRICT ac = a_ptr + j * m;
                const real_t x0 = x_ptr[j + r * n];
                for (size_t i = i0; i < i1; ++i) {
                  y0[i] += ac[i] * x0;
                }
              }
            }
          }
        }
      }
      return;
    }

    // L > 1: unit-stride inner loop over l, 4-way register-blocked over i.
    const size_t x_stride = L * n;
    const size_t y_stride = L * m;
    for (size_t r = 0; r < R; ++r) {
      const real_t *MFEM_KRON_RESTRICT xr = x_ptr + r * x_stride;
      real_t *MFEM_KRON_RESTRICT yr = y_ptr + r * y_stride;
      for (size_t i0 = 0; i0 < m; i0 += kRBlk) {
        const size_t i1 = std::min(i0 + kRBlk, m);
        for (size_t l0 = 0; l0 < L; l0 += kLBlk) {
          const size_t l1 = std::min(l0 + kLBlk, L);
          for (size_t j0 = 0; j0 < n; j0 += kKBlk) {
            const size_t j1 = std::min(j0 + kKBlk, n);
            size_t i = i0;
            for (; i + 4 <= i1; i += 4) {
              real_t *MFEM_KRON_RESTRICT y0 = yr + (i + 0) * L;
              real_t *MFEM_KRON_RESTRICT y1 = yr + (i + 1) * L;
              real_t *MFEM_KRON_RESTRICT y2 = yr + (i + 2) * L;
              real_t *MFEM_KRON_RESTRICT y3 = yr + (i + 3) * L;
              if (j0 == 0) {
                for (size_t l = l0; l < l1; ++l) {
                  y0[l] = real_t(0);
                  y1[l] = real_t(0);
                  y2[l] = real_t(0);
                  y3[l] = real_t(0);
                }
              }
              for (size_t j = j0; j < j1; ++j) {
                const real_t *MFEM_KRON_RESTRICT xc = xr + j * L;
                const real_t a0 = a_ptr[(i + 0) + j * m];
                const real_t a1 = a_ptr[(i + 1) + j * m];
                const real_t a2 = a_ptr[(i + 2) + j * m];
                const real_t a3 = a_ptr[(i + 3) + j * m];
                for (size_t l = l0; l < l1; ++l) {
                  const real_t x = xc[l];
                  y0[l] += x * a0;
                  y1[l] += x * a1;
                  y2[l] += x * a2;
                  y3[l] += x * a3;
                }
              }
            }
            for (; i < i1; ++i) {
              real_t *MFEM_KRON_RESTRICT y0 = yr + i * L;
              if (j0 == 0) {
                for (size_t l = l0; l < l1; ++l) {
                  y0[l] = real_t(0);
                }
              }
              for (size_t j = j0; j < j1; ++j) {
                const real_t *MFEM_KRON_RESTRICT xc = xr + j * L;
                const real_t a0 = a_ptr[i + j * m];
                for (size_t l = l0; l < l1; ++l) {
                  y0[l] += xc[l] * a0;
                }
              }
            }
          }
        }
      }
    }
  }

  void EnsureEOScratch(size_t in_size, size_t out_size) const {
    if (static_cast<size_t>(eo_in.Size()) < in_size) {
      eo_in.SetSize(static_cast<int>(in_size));
    }
    if (static_cast<size_t>(eo_out.Size()) < out_size) {
      eo_out.SetSize(static_cast<int>(out_size));
    }
  }

  void ApplyEvenOdd(size_t L, size_t R, size_t m, size_t n,
                    const EvenOddFactor &eo, const real_t *x_ptr,
                    real_t *y_ptr) const {
    const size_t m_o = m / 2;
    const size_t m_e = (m + 1) / 2;
    const size_t n_o = n / 2;
    const size_t n_e = (n + 1) / 2;

    if (L == 1) {
      constexpr size_t kRChunk = 128;
      const size_t max_rc = std::min(R, kRChunk);
      EnsureEOScratch(n * max_rc, m * max_rc);

      real_t *in_buf = eo_in.HostWrite();
      real_t *out_buf = eo_out.HostWrite();

      for (size_t r0 = 0; r0 < R; r0 += kRChunk) {
        const size_t rc = std::min(kRChunk, R - r0);
        const real_t *xb = x_ptr + r0 * n;
        real_t *yb = y_ptr + r0 * m;

        if (eo.type == SymmetryType::RowEvenOdd) {
          const size_t pe = eo.even_idx.size();
          const size_t po = eo.odd_idx.size();
          real_t *xe = in_buf;
          real_t *xo = in_buf + n_e * rc;
          real_t *ye = out_buf;
          real_t *yo = out_buf + pe * rc;

          for (size_t c = 0; c < rc; ++c) {
            const real_t *xc = xb + c * n;
            real_t *xec = xe + c * n_e;
            real_t *xoc = xo + c * n_o;
            for (size_t j = 0; j < n_o; ++j) {
              const real_t a = xc[j];
              const real_t b = xc[n - 1 - j];
              xec[j] = a + b;
              xoc[j] = a - b;
            }
            if (n % 2 == 1) {
              xec[n_o] = xc[n_o];
            }
          }

          if (pe > 0) {
            ContractDense(1, rc, pe, n_e, eo.a_even.HostRead(), xe, ye);
          }
          if (po > 0) {
            ContractDense(1, rc, po, n_o, eo.a_odd.HostRead(), xo, yo);
          }

          for (size_t c = 0; c < rc; ++c) {
            real_t *yc = yb + c * m;
            const real_t *yec = ye + c * pe;
            const real_t *yoc = yo + c * po;
            for (size_t r = 0; r < pe; ++r) {
              yc[eo.even_idx[r]] = yec[r];
            }
            for (size_t r = 0; r < po; ++r) {
              yc[eo.odd_idx[r]] = yoc[r];
            }
          }
        } else if (eo.type == SymmetryType::ColEvenOdd) {
          const size_t qe = eo.even_idx.size();
          const size_t qo = eo.odd_idx.size();
          real_t *xe = in_buf;
          real_t *xo = in_buf + qe * rc;
          real_t *ye = out_buf;
          real_t *yo = out_buf + m_e * rc;

          for (size_t c = 0; c < rc; ++c) {
            const real_t *xc = xb + c * n;
            real_t *xec = xe + c * qe;
            real_t *xoc = xo + c * qo;
            for (size_t k = 0; k < qe; ++k) {
              xec[k] = xc[eo.even_idx[k]];
            }
            for (size_t k = 0; k < qo; ++k) {
              xoc[k] = xc[eo.odd_idx[k]];
            }
          }

          if (qe > 0) {
            ContractDense(1, rc, m_e, qe, eo.a_even.HostRead(), xe, ye);
          } else {
            std::fill_n(ye, m_e * rc, real_t(0));
          }
          if (qo > 0) {
            ContractDense(1, rc, m_o, qo, eo.a_odd.HostRead(), xo, yo);
          } else {
            std::fill_n(yo, m_o * rc, real_t(0));
          }

          for (size_t c = 0; c < rc; ++c) {
            real_t *yc = yb + c * m;
            const real_t *yec = ye + c * m_e;
            const real_t *yoc = yo + c * m_o;
            for (size_t i = 0; i < m_o; ++i) {
              const real_t e = yec[i];
              const real_t o = yoc[i];
              yc[i] = e + o;
              yc[m - 1 - i] = e - o;
            }
            if (m % 2 == 1) {
              yc[m_o] = yec[m_o];
            }
          }
        } else // Centrosymmetric
        {
          real_t *xe = in_buf;
          real_t *xo = in_buf + n_e * rc;
          real_t *ye = out_buf;
          real_t *yo = out_buf + m_e * rc;

          for (size_t c = 0; c < rc; ++c) {
            const real_t *xc = xb + c * n;
            real_t *xec = xe + c * n_e;
            real_t *xoc = xo + c * n_o;
            for (size_t j = 0; j < n_o; ++j) {
              const real_t a = xc[j];
              const real_t b = xc[n - 1 - j];
              xec[j] = a + b;
              xoc[j] = a - b;
            }
            if (n % 2 == 1) {
              xec[n_o] = xc[n_o];
            }
          }

          ContractDense(1, rc, m_e, n_e, eo.a_even.HostRead(), xe, ye);
          ContractDense(1, rc, m_o, n_o, eo.a_odd.HostRead(), xo, yo);

          for (size_t c = 0; c < rc; ++c) {
            real_t *yc = yb + c * m;
            const real_t *yec = ye + c * m_e;
            const real_t *yoc = yo + c * m_o;
            for (size_t i = 0; i < m_o; ++i) {
              const real_t e = yec[i];
              const real_t o = yoc[i];
              yc[i] = e + o;
              yc[m - 1 - i] = e - o;
            }
            if (m % 2 == 1) {
              yc[m_o] = yec[m_o];
            }
          }
        }
      }
      return;
    }

    // L > 1: process slice-by-slice to keep fold/unfold buffers cache-hot.
    EnsureEOScratch(L * n, L * m);
    real_t *in_buf = eo_in.HostWrite();
    real_t *out_buf = eo_out.HostWrite();
    const size_t x_stride = L * n;
    const size_t y_stride = L * m;

    for (size_t r = 0; r < R; ++r) {
      const real_t *xr = x_ptr + r * x_stride;
      real_t *yr = y_ptr + r * y_stride;

      if (eo.type == SymmetryType::RowEvenOdd) {
        const size_t pe = eo.even_idx.size();
        const size_t po = eo.odd_idx.size();
        real_t *xe = in_buf;
        real_t *xo = in_buf + L * n_e;
        real_t *ye = out_buf;
        real_t *yo = out_buf + L * pe;

        for (size_t j = 0; j < n_o; ++j) {
          const real_t *xa = xr + j * L;
          const real_t *xb = xr + (n - 1 - j) * L;
          real_t *xec = xe + j * L;
          real_t *xoc = xo + j * L;
          for (size_t l = 0; l < L; ++l) {
            xec[l] = xa[l] + xb[l];
            xoc[l] = xa[l] - xb[l];
          }
        }
        if (n % 2 == 1) {
          std::copy_n(xr + n_o * L, L, xe + n_o * L);
        }

        if (pe > 0) {
          ContractDense(L, 1, pe, n_e, eo.a_even.HostRead(), xe, ye);
        }
        if (po > 0) {
          ContractDense(L, 1, po, n_o, eo.a_odd.HostRead(), xo, yo);
        }

        for (size_t k = 0; k < pe; ++k) {
          std::copy_n(ye + k * L, L, yr + eo.even_idx[k] * L);
        }
        for (size_t k = 0; k < po; ++k) {
          std::copy_n(yo + k * L, L, yr + eo.odd_idx[k] * L);
        }
      } else if (eo.type == SymmetryType::ColEvenOdd) {
        const size_t qe = eo.even_idx.size();
        const size_t qo = eo.odd_idx.size();
        real_t *xe = in_buf;
        real_t *xo = in_buf + L * qe;
        real_t *ye = out_buf;
        real_t *yo = out_buf + L * m_e;

        for (size_t k = 0; k < qe; ++k) {
          std::copy_n(xr + eo.even_idx[k] * L, L, xe + k * L);
        }
        for (size_t k = 0; k < qo; ++k) {
          std::copy_n(xr + eo.odd_idx[k] * L, L, xo + k * L);
        }

        if (qe > 0) {
          ContractDense(L, 1, m_e, qe, eo.a_even.HostRead(), xe, ye);
        } else {
          std::fill_n(ye, L * m_e, real_t(0));
        }
        if (qo > 0) {
          ContractDense(L, 1, m_o, qo, eo.a_odd.HostRead(), xo, yo);
        } else {
          std::fill_n(yo, L * m_o, real_t(0));
        }

        for (size_t i = 0; i < m_o; ++i) {
          const real_t *yec = ye + i * L;
          const real_t *yoc = yo + i * L;
          real_t *ya = yr + i * L;
          real_t *yb = yr + (m - 1 - i) * L;
          for (size_t l = 0; l < L; ++l) {
            ya[l] = yec[l] + yoc[l];
            yb[l] = yec[l] - yoc[l];
          }
        }
        if (m % 2 == 1) {
          std::copy_n(ye + m_o * L, L, yr + m_o * L);
        }
      } else // Centrosymmetric
      {
        real_t *xe = in_buf;
        real_t *xo = in_buf + L * n_e;
        real_t *ye = out_buf;
        real_t *yo = out_buf + L * m_e;

        for (size_t j = 0; j < n_o; ++j) {
          const real_t *xa = xr + j * L;
          const real_t *xb = xr + (n - 1 - j) * L;
          real_t *xec = xe + j * L;
          real_t *xoc = xo + j * L;
          for (size_t l = 0; l < L; ++l) {
            xec[l] = xa[l] + xb[l];
            xoc[l] = xa[l] - xb[l];
          }
        }
        if (n % 2 == 1) {
          std::copy_n(xr + n_o * L, L, xe + n_o * L);
        }

        ContractDense(L, 1, m_e, n_e, eo.a_even.HostRead(), xe, ye);
        ContractDense(L, 1, m_o, n_o, eo.a_odd.HostRead(), xo, yo);

        for (size_t i = 0; i < m_o; ++i) {
          const real_t *yec = ye + i * L;
          const real_t *yoc = yo + i * L;
          real_t *ya = yr + i * L;
          real_t *yb = yr + (m - 1 - i) * L;
          for (size_t l = 0; l < L; ++l) {
            ya[l] = yec[l] + yoc[l];
            yb[l] = yec[l] - yoc[l];
          }
        }
        if (m % 2 == 1) {
          std::copy_n(ye + m_o * L, L, yr + m_o * L);
        }
      }
    }
  }
};

template <typename... T>
KroneckerProductOperator(const T &...)
    -> KroneckerProductOperator<sizeof...(T)>;

} // namespace mfem

#undef MFEM_KRON_RESTRICT

#endif // MFEM_KRONECKER
