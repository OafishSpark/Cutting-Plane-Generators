#ifndef ESOLVER_SPARSE_ROW_MATRIX_H
#define ESOLVER_SPARSE_ROW_MATRIX_H

#include <cstdint>

#include "sparse_vector.h"

class SparseRowMatrix {
public:
    explicit SparseRowMatrix(size_t n_rows = 0, size_t n_cols = 0);

    std::vector<Scalar> operator*(const std::vector<Scalar> &x) const;

    void Resize(size_t n_rows, size_t n_cols) {
        n_rows_ = n_rows;
        n_cols_ = n_cols;
        rows_.resize(n_rows);
    }

    size_t GetNRows() const { return n_rows_; }

    size_t GetNCols() const { return n_cols_; }

    SparseVector &GetRow(Index i_row) { return rows_[i_row]; }
    const SparseVector &GetRow(Index i_row) const { return rows_[i_row]; }

    inline std::vector<SparseVector>& GetRows() { return rows_; }

private:
    size_t n_rows_;
    size_t n_cols_;
    std::vector<SparseVector> rows_;
};

#endif  // ESOLVER_SPARSE_ROW_MATRIX_H
