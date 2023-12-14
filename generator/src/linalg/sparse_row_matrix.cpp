#include "../../include/linalg/sparse_row_matrix.h"

#include <cassert>

SparseRowMatrix::SparseRowMatrix(size_t n_rows, size_t n_cols) : n_rows_(n_rows), n_cols_(n_cols) {
    rows_.resize(n_rows_);
    for (Index i_row = 0; i_row < Index(n_rows_); i_row++) {
        rows_[i_row] = SparseVector();
    }
}

std::vector<Scalar> SparseRowMatrix::operator*(const std::vector<Scalar> &x) const {
    assert(n_cols_ == x.size());
    std::vector<Scalar> res(n_rows_, Scalar(0));
    for (Index i_row = 0; i_row < Index(n_rows_); i_row++) {
        for (Index i=0; i<rows_[i_row].size(); i++) {
            res[i_row] += rows_[i_row].values[i] * x[rows_[i_row].index[i]];
        }
    }
    return res;
}
