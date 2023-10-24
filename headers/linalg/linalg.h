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
        Index ind_;
    };

    std::vector<Elem> elements_;
    Index len_;

    SparseVector(std::vector<Scalar>& dense_vector);
    SparseVector(Index len);

    void Print();
};


class SparseColMatrix {
public:
    std::vector<SparseVector> matrix_;
    Index rows_;
    Index cols_;

    SparseColMatrix(DenseMatrix& dense_matrix);

    SparseColMatrix(){}

    void Print();
};


class RHS {
public:
    struct Elem {
        char type_;
        Scalar val_;
    };

    std::vector<Elem> rhs_;
};


class Variables {
public:
    struct Elem {
        bool is_int_;
        Scalar bnd_lo_;
        Scalar bnd_up_;
    };

    std::vector<Elem> vars_;
};

#endif //LINALG_H
