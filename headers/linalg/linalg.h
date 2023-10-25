#ifndef LINALG_H
#define LINALG_H

#include "../utils/utils.h"

#include <vector>
#include <iostream>
#include <cassert>


class SparseVector {
public:
    struct Elem {
        Scalar val_;
        int ind_;
    };

    std::vector<Elem> elements_;
    int len_;

    SparseVector(std::vector<Scalar>& dense_vector);
    SparseVector(int len);

    void Print();

    SparseVector::Elem operator[](int ind) { return elements_[ind]; }
};

Scalar operator*(std::vector<Scalar>, SparseVector);

Scalar operator*(SparseVector, SparseVector);


class SparseColMatrix {
public:
    std::vector<SparseVector> matrix_;
    int rows_;
    int cols_;

    SparseColMatrix(DenseMatrix& dense_matrix);

    SparseColMatrix(){}

    void Print();

    SparseVector operator[](int ind) {
        assert(ind < cols_);
        return matrix_[ind];
    }
};


class RHS {
public:
    struct Elem {
        char type_;
        Scalar val_;
    };

    std::vector<Elem> rhs_;

    RHS::Elem operator[](int ind) { return rhs_[ind]; }
};


class Variables {
public:
    struct Elem {
        bool is_int_;
        Scalar bnd_lo_;
        Scalar bnd_up_;
    };

    std::vector<Elem> vars_;

    Variables::Elem operator[](int ind) { return vars_[ind]; }
};

#endif //LINALG_H
