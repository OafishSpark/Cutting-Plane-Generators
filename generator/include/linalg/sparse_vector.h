#ifndef ESOLVER_SPARSE_VECTOR_H
#define ESOLVER_SPARSE_VECTOR_H

#include <cassert>
#include <cmath>
#include <cstdint>
#include <vector>

#include "scalars/scalar.h"

class SparseVector {
   public:
    SparseVector();

    explicit SparseVector(const std::vector<Scalar> &x_dense);

    SparseVector(const Scalar *x_dense, size_t size);

    inline Index size() const {
        assert(index.size() == values.size());
        return Index(index.size());
    }

    inline void Resize(size_t n) {
        index.resize(n);
        values.resize(n);
    }
    inline void Reserve(size_t n) {
        index.reserve(n);
        values.reserve(n);
    }
    inline void PushBack(const Index i, const Scalar v) {
        index.push_back(i);
        values.push_back(v);
    }

    std::vector<Index> index;
    std::vector<Scalar> values;

    SparseVector operator-() const;

    SparseVector operator*(Scalar x) const;
    Scalar operator*(const SparseVector &x) const;
    Scalar operator*(const Scalar *x) const;
    Scalar operator*(const std::vector<Scalar> &x) const;

    SparseVector operator+(const SparseVector &x) const;
    SparseVector operator-(const SparseVector &x) const;

    void operator*=(Scalar x);

    size_t RemoveElements(const std::vector<Index> &i);

    void Sort();
};

#endif  // ESOLVER_SPARSE_VECTOR_H
