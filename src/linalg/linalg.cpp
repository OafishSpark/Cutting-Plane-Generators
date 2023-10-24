#include "../../headers/linalg/linalg.h"


// Sparse vector part

SparseVector::SparseVector(std::vector<Scalar>& dense_vector) {
    len_ = dense_vector.size();
    Index nonzero_ind = 0;
    elements_.resize(len_, Elem({0.0, 0}));
    for (Index iv = 0; iv < len_; iv++) {
        Scalar elem = dense_vector[iv]; 
        if (abs(elem) >= Zero_Epsilon) {
            elements_.at(nonzero_ind).val_ = elem;
            elements_.at(nonzero_ind).ind_ = iv;
            nonzero_ind += 1;
        }
    }
    elements_.resize(nonzero_ind);
}

SparseVector::SparseVector(Index len) {
    len_ = len;
}


void SparseVector::Print() {
    for (auto elem : elements_) {
        std::cout << "(" << elem.ind_ << ", " << elem.val_ << ")";       
    }
    std::cout << " len=" << len_ << "\n";
}

// SparseMatrix part

SparseColMatrix::SparseColMatrix(DenseMatrix& dense_matrix) {
    rows_ = dense_matrix.size();
    assert(rows_ > 0);
    cols_ = dense_matrix.back().size();
    assert(cols_ > 0);
    for (const auto& row: dense_matrix) {
        assert(row.size() == cols_);
    }
    matrix_.resize(cols_, SparseVector(rows_));
    for (Index iv = 0; iv < cols_; ++iv) {
        matrix_[iv].elements_.resize(rows_, SparseVector::Elem({0.0, 0}));
        Index nonzero_ind = 0;
        for (Index jv = 0; jv < rows_; ++jv) {
            Scalar elem = dense_matrix[jv][iv]; 
            if (abs(elem) >= Zero_Epsilon) {
                matrix_[iv].elements_.at(nonzero_ind).val_ = elem;
                matrix_[iv].elements_.at(nonzero_ind).ind_ = jv;
                nonzero_ind += 1;
            }
        }
        matrix_[iv].elements_.resize(nonzero_ind);
    }
}

void SparseColMatrix::Print() {
    for (Index iv = 0; iv < cols_; ++iv) {
        std::cout << "Column " << iv << ": ";
        matrix_[iv].Print();
    }
    std::cout << " cols=" << cols_ << " rows=" << rows_ << "\n";
}
