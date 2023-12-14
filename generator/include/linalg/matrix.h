#ifndef ESOLVER_MATRIX_H
#define ESOLVER_MATRIX_H

#include <cstdint>
#include <vector>

#include "sparse_col_matrix.h"

class Matrix {
   public:
    // zero matrix
    explicit Matrix(size_t rows_cnt_ = 0, size_t cols_cnt_ = 0);

    explicit Matrix(SparseColMatrix &);

    Matrix(const Matrix &);

    void ResizeAsZero(size_t rows_cnt_, size_t cols_cnt_);

    const std::vector<Scalar> &GetData() const;
    std::vector<Scalar> &GetMutableData();

    std::vector<Scalar> GetRow(Index row) const;
    Scalar *GetRowPointer(Index row) const;
    std::vector<Scalar> GetCol(Index col) const;

    Scalar &operator()(Index row, Index col);
    inline Scalar operator()(Index row, Index col) const;

    std::vector<Scalar> operator*(const std::vector<Scalar> &x) const;
    std::vector<Scalar> operator*(const Scalar *x) const;
    std::vector<Scalar> operator*(const SparseVector &x) const;

    void Inverse();

    Matrix GetInverse() const;

    size_t rows_cnt;
    size_t cols_cnt;

   private:
    std::vector<Scalar> data_;
};

std::vector<Scalar> operator*(std::vector<Scalar> &x, Matrix &m);
std::vector<Scalar> operator*(const Scalar *x, Matrix &m);

#endif  // ESOLVER_MATRIX_H
