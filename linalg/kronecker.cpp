#include "config/config.hpp"
#include "general/error.hpp"
#include "linalg/densemat.hpp"
#include "operator.hpp"
#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <numeric>
#include <type_traits>

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
};

template <typename... T>
KroneckerProductOperator(const T &...)
    -> KroneckerProductOperator<sizeof...(T)>;
} // namespace mfem
