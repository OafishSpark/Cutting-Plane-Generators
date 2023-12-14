#include "../../include/linalg/matrix.h"

// #include <lapacke.h>

#include <cassert>
#include <cstring>

#include "linalg_utils.h"
#include "lup.h"
#include "sparse_col_matrix.h"
#include "sparse_vector.h"

Matrix::Matrix(size_t rows_cnt_, size_t cols_cnt_)
    : rows_cnt(rows_cnt_), cols_cnt(cols_cnt_), data_(rows_cnt_ * cols_cnt_, Scalar(0)) {}

Matrix::Matrix(const Matrix &other) = default;

// Remember: Matrix is in RowMajor format
Matrix::Matrix(SparseColMatrix &A)
    : rows_cnt(A.GetNRows()), cols_cnt(A.GetNCols()), data_(rows_cnt * cols_cnt) {
    for (Index i_col = 0; i_col < Index(cols_cnt); ++i_col) {
        auto &col = A.GetCol(i_col);
        for (Index i = 0; i < col.size(); i++) {
            data_[col.index[i] * cols_cnt + i_col] = col.values[i];
        }
    }
}

void Matrix::ResizeAsZero(size_t rows_cnt_, size_t cols_cnt_) {
    if (rows_cnt_ * cols_cnt_ > rows_cnt * cols_cnt) {
        data_.resize(rows_cnt_ * cols_cnt_);
    }
    std::fill(data_.begin(), data_.end(), Scalar(0));
    rows_cnt = rows_cnt_;
    cols_cnt = cols_cnt_;
}

const std::vector<Scalar> &Matrix::GetData() const { return data_; }

std::vector<Scalar> &Matrix::GetMutableData() { return data_; }

std::vector<Scalar> Matrix::GetRow(Index row) const {
    std::vector<Scalar> x(cols_cnt);
    memcpy(x.data(), data_.data() + row * cols_cnt, cols_cnt * sizeof(Scalar));
    return x;
}

Scalar *Matrix::GetRowPointer(Index row) const {
    assert(row < Index(rows_cnt));
    return const_cast<Scalar *>(data_.data() + row * cols_cnt);
}

std::vector<Scalar> Matrix::GetCol(Index col) const {
    std::vector<Scalar> x(rows_cnt);
    for (Index ic = 0; ic < Index(rows_cnt); ic++) {
        x[ic] = operator()(ic, col);
    }
    return x;
}

Scalar &Matrix::operator()(Index row, Index col) {
    assert(row >= 0 && row < Index(rows_cnt));
    assert(col >= 0 && col < Index(cols_cnt));
    return data_[row * cols_cnt + col];
}

Scalar Matrix::operator()(Index row, Index col) const {
    assert(row >= 0 && row < Index(rows_cnt));
    assert(col >= 0 && col < Index(cols_cnt));
    return data_[row * cols_cnt + col];
}

std::vector<Scalar> Matrix::operator*(const std::vector<Scalar> &x) const {
    assert(x.size() == cols_cnt);
    return operator*(x.data());
}

std::vector<Scalar> Matrix::operator*(const Scalar *x) const {
    std::vector<Scalar> res(rows_cnt, Scalar(0));
    for (Index ic = 0; ic < Index(rows_cnt); ++ic) {
        dot(Index(cols_cnt), GetRowPointer(ic), x, &res[ic]);
    }
    return res;
}

std::vector<Scalar> Matrix::operator*(const SparseVector &x) const {
    std::vector<Scalar> res(rows_cnt, Scalar(0));
    for (Index ic = 0; ic < Index(rows_cnt); ++ic) {
        for (Index i = 0; i < x.size(); ++i) {
            res[ic] += operator()(ic, x.index[i]) * x.values[i];
        }
    }
    return res;
}

std::vector<Scalar> operator*(std::vector<Scalar> &x, Matrix &m) {
    assert(x.size() == m.rows_cnt);
    return operator*(x.data(), m);
}

std::vector<Scalar> operator*(const Scalar *x, Matrix &m) {
    std::vector<Scalar> res(m.cols_cnt, Scalar(0));
    auto *x_p = const_cast<Scalar *>(x);
    gemv(m.rows_cnt, m.cols_cnt, Scalar(1), &m(0, 0), x_p, Scalar(0), &res[0]);
    return res;
}

Matrix Matrix::GetInverse() const {
    assert(cols_cnt == rows_cnt);
    Matrix inv(*this);
    inv.Inverse();
    return inv;
}

void Matrix::Inverse() {
    assert(cols_cnt == rows_cnt);

    auto n = Index(cols_cnt);
    std::vector<Index> p(n);

    LUP lup;
    lup.lu(n, n, data_.data(), n, p.data());
    lup.inverse(n, data_.data(), n, p.data());
}
