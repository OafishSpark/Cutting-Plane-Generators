#ifndef ESOLVER_SPARSE_COL_MATRIX_H
#define ESOLVER_SPARSE_COL_MATRIX_H

#include <cstdint>

#include "sparse_vector.h"

class SparseColMatrix {
public:
    explicit SparseColMatrix(size_t n_rows = 0, size_t n_cols = 0);

    std::vector<Scalar> operator*(const std::vector<Scalar> &x) const;

    void Resize(size_t n_rows, size_t n_cols) {
        n_rows_ = n_rows;
        n_cols_ = n_cols;
        cols_.resize(n_cols);
    }

    size_t GetNRows() const { return n_rows_; }

    size_t GetNCols() const { return n_cols_; }

    SparseVector &GetCol(Index i_col) { return cols_[i_col]; }
    const SparseVector &GetCol(Index i_col) const { return cols_[i_col]; }

    inline std::vector<SparseVector>& GetCols() { return cols_; }

private:
    size_t n_rows_;
    size_t n_cols_;
    std::vector<SparseVector> cols_;
};

std::vector<Scalar> operator*(const Scalar *x, const SparseColMatrix &m);

std::vector<Scalar> operator*(const std::vector<Scalar> &x, const SparseColMatrix &m);

std::vector<Scalar> operator*(const SparseVector &x, const SparseColMatrix &m);

#endif  // ESOLVER_SPARSE_COL_MATRIX_H
