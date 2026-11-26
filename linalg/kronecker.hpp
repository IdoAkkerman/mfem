
#include "linalg/densemat.hpp"
#include "operator.hpp"
#include <assert.h>
#include <cstdio>
#include <numeric>

namespace mfem {

template <size_t N> class KroneckerProductOperator : public Operator {
  static_assert(N > 0,
                "KroneckerProductOperator requires at least one matrix.");

protected:
  std::array<mfem::DenseMatrix, N> matrices;

public:
  // Passing (D_1, ..., D_N) computes (D_N ⊗ ... ⊗ D_1) x
  template <typename... T>
  KroneckerProductOperator(const T &...matrices)
      : matrices{matrices...},
        X(std::make_unique<Vector>(
            (std::max(matrices.Height(), matrices.Width()) *
             ...))), // Note: this overallocates in general.
        X_prime(std::make_unique<Vector>(
            (std::max(matrices.Height(), matrices.Width()) * ...))),
        Operator((matrices.Height() * ...), (matrices.Width() * ...)) {
    static_assert(
        (std::is_same_v<T, mfem::DenseMatrix> && ...),
        "KroneckerProductOperator only supports DenseMatrix type arguments.");
    static_assert(N == sizeof...(T),
                  "Operator size and matrix count is not the same.");
  }

  void Mult(const Vector &x, Vector &y) const override {
    MFEM_ASSERT(x.Size() == this->Width(),
                "Vector size and Matrix width must match!");
    MFEM_ASSERT(y.Size() == this->Height(),
                "Output vector size and Matrix height must match!");
    memcpy(X_prime->GetData(), x.GetData(), x.Size() * sizeof(real_t));
    for (size_t mode = 0; mode < N; ++mode) {
      std::swap(X, X_prime);
      modewise_contract_tensor(mode);
    }
    memcpy(y.GetData(), X_prime->GetData(), y.Size() * sizeof(real_t));
  };
  /// Compute y = (D_{N-1}^T \otimes ... \otimes D_0^T) x.
  void MultTranspose(const Vector &x, Vector &y) const override {
    MFEM_VERIFY(x.Size() == height,
                "Input vector size and operator height must match.");
    MFEM_VERIFY(y.Size() == width,
                "Output vector size and operator width must match.");

    std::copy_n(x.HostRead(), x.Size(), X_prime->HostWrite());
    for (size_t mode = 0; mode < N; ++mode) {
      std::swap(X, X_prime);
      modewise_contract_tensor_transpose(mode);
    }
    std::copy_n(X_prime->HostRead(), y.Size(), y.HostWrite());
  }

private:
  mutable std::unique_ptr<Vector> X;
  mutable std::unique_ptr<Vector> X_prime;
  /*
   * Invariant: X contains X to be used to create the next X'.
   * Output in X_prime
   * */
  void modewise_contract_tensor(size_t mode) const {
    const size_t L{
        std::accumulate(matrices.begin(), matrices.begin() + mode, size_t{1},
                        [](size_t acc, const DenseMatrix &m) {
                          return acc *= static_cast<size_t>(m.Height());
                        })};
    const size_t R{
        std::accumulate(matrices.begin() + mode + 1, matrices.end(), size_t{1},
                        [](size_t acc, const DenseMatrix &m) {
                          return acc *= static_cast<size_t>(m.Width());
                        })};
    const DenseMatrix &D_k = matrices[mode];
    const size_t m_k{static_cast<size_t>(D_k.Height())};
    const size_t n_k{static_cast<size_t>(D_k.Width())};

    // TODO: Optimize this
    for (size_t l{0}; l < L; ++l) {
      for (size_t r{0}; r < R; ++r) {
        for (size_t i{0}; i < m_k; ++i) {
          real_t acc{0};
          for (size_t j{0}; j < n_k; ++j) {
            acc += D_k.Elem(i, j) * X->Elem(l + j * L + r * L * n_k);
          }
          X_prime->Elem(l + i * L + r * L * m_k) = acc;
        }
      }
    }
  }
  void modewise_contract_tensor_transpose(size_t mode) const {
    const size_t L =
        std::accumulate(matrices.begin(), matrices.begin() + mode, size_t{1},
                        [](size_t acc, const DenseMatrix &m) {
                          return acc * static_cast<size_t>(m.Width());
                        });
    const size_t R =
        std::accumulate(matrices.begin() + mode + 1, matrices.end(), size_t{1},
                        [](size_t acc, const DenseMatrix &m) {
                          return acc * static_cast<size_t>(m.Height());
                        });

    const DenseMatrix &D_k = matrices[mode];
    const size_t m_k = static_cast<size_t>(D_k.Height());
    const size_t n_k = static_cast<size_t>(D_k.Width());

    const real_t *x_data = X->HostRead();
    real_t *xp_data = X_prime->HostWrite();

    for (size_t r = 0; r < R; ++r) {
      const size_t x_r_offset = r * L * m_k;
      const size_t xp_r_offset = r * L * n_k;
      for (size_t i = 0; i < n_k; ++i) {
        const int col = static_cast<int>(i);
        const size_t xp_ri_offset = xp_r_offset + i * L;
        for (size_t l = 0; l < L; ++l) {
          real_t acc = 0.0;
          for (size_t j = 0; j < m_k; ++j) {
            const real_t d_ji = D_k.Elem(static_cast<int>(j), col);
            acc += d_ji * x_data[l + j * L + x_r_offset];
          }
          xp_data[l + xp_ri_offset] = acc;
        }
      }
    }
  }
};

template <typename... T>
KroneckerProductOperator(const T &...)
    -> KroneckerProductOperator<sizeof...(T)>;
} // namespace mfem
