#include "../../headers/linalg/linalg.h"


// Sparse vector part

SparseVector::SparseVector(std::vector<Scalar>& dense_vector) {
    len_ = dense_vector.size();
    int nonzero_ind = 0;
    elements_.resize(len_, Elem({0.0, 0}));
    for (int iv = 0; iv < len_; iv++) {
        Scalar elem = dense_vector[iv]; 
        if (Abs(elem) >= kZero_Epsilon) {
            elements_.at(nonzero_ind).val_ = elem;
            elements_.at(nonzero_ind).ind_ = iv;
            nonzero_ind += 1;
        }
    }
    elements_.resize(nonzero_ind);
}

SparseVector::SparseVector(int len) {
    len_ = len;
}

void SparseVector::Print() {
    for (auto elem : elements_) {
        std::cout << "(" << elem.ind_ << ", " << elem.val_ << ")";       
    }
    std::cout << " len=" << len_ << "\n";
}

std::string SparseVector::PrintInFile() {
    std::string answer = "";
    for (auto elem : elements_) {
        answer += " " + std::to_string(elem.ind_) + "," + std::to_string(elem.val_);       
    }
    return answer;
}

Scalar operator*(std::vector<Scalar> dns_vec, SparseVector sprs_vec) {
    assert(dns_vec.size() == sprs_vec.len_);
    Scalar answer = 0;
    for (const auto& elem: sprs_vec.elements_) {
        answer += elem.val_ * dns_vec.at(elem.ind_);
    }
    return answer;
}

Scalar operator*(SparseVector vec1, SparseVector vec2) {
    assert(vec1.len_ == vec2.len_);
    Scalar answer = 0;
    auto iter1 = vec1.elements_.begin();
    auto iter2 = vec2.elements_.begin();
    while ((iter1 != vec1.elements_.end()) && (iter2 != vec2.elements_.end())) {
        if (iter1->ind_ == iter2->ind_) {
            answer += iter1->val_ * iter2->val_;
            iter1++;
            iter2++;
            continue;
        } else if (iter1->ind_ > iter2->ind_) {
            iter2++;
        } else {
            iter1++;
        }
    } 
    return answer;
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
    for (int iv = 0; iv < cols_; ++iv) {
        matrix_[iv].elements_.resize(rows_, SparseVector::Elem({0.0, 0}));
        int nonzero_ind = 0;
        for (int jv = 0; jv < rows_; ++jv) {
            Scalar elem = dense_matrix[jv][iv]; 
            if (Abs(elem) >= kZero_Epsilon) {
                matrix_[iv].elements_.at(nonzero_ind).val_ = elem;
                matrix_[iv].elements_.at(nonzero_ind).ind_ = jv;
                nonzero_ind += 1;
            }
        }
        matrix_[iv].elements_.resize(nonzero_ind);
    }
}

void SparseColMatrix::Print() {
    for (int iv = 0; iv < cols_; ++iv) {
        std::cout << "Column " << iv << ": ";
        matrix_[iv].Print();
    }
    std::cout << " cols=" << cols_ << " rows=" << rows_ << "\n";
}
