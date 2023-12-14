#include "../../include/linalg/sparse_col_matrix.h"

#include <cassert>

SparseColMatrix::SparseColMatrix(size_t n_rows, size_t n_cols) : n_rows_(n_rows), n_cols_(n_cols) {
    cols_.resize(n_cols);
    for (Index i_col = 0; i_col < Index(n_cols); i_col++) {
        cols_[i_col] = SparseVector();
    }
}

std::vector<Scalar> SparseColMatrix::operator*(const std::vector<Scalar> &x) const {
    assert(n_cols_ == x.size());
    std::vector<Scalar> res(n_rows_, Scalar(0));
    for (Index i_col = 0; i_col < Index(n_cols_); i_col++) {
        if (x[i_col] == Scalar(0)) {
            continue;
        }
        for (Index i=0; i<Index(cols_[i_col].size()); i++) {
            res[cols_[i_col].index[i]] += x[i_col] * cols_[i_col].values[i];
        }
    }
    return res;
}

std::vector<Scalar> operator*(const Scalar *x, const SparseColMatrix &m) {
    std::vector<Scalar> res(m.GetNCols(), Scalar(0));
    for (Index i_col = 0; i_col < Index(m.GetNCols()); i_col++) {
        for (Index i=0; i<Index(m.GetCol(i_col).size()); i++) {
            res[i_col] += m.GetCol(i_col).values[i] * x[m.GetCol(i_col).index[i]];
        }
    }
    return res;
}

std::vector<Scalar> operator*(const std::vector<Scalar> &x, const SparseColMatrix &m) {
    return x.data() * m;
}

std::vector<Scalar> operator*(const SparseVector &x, const SparseColMatrix &m) {
    std::vector<Scalar> res(m.GetNCols(), Scalar(0));
    for (Index i_col = 0; i_col < Index(m.GetNCols()); i_col++) {
        res[i_col] = x * m.GetCol(i_col);
    }
    return res;
}
