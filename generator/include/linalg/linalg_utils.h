#ifndef ESOLVER_LINALG_UTILS_H
#define ESOLVER_LINALG_UTILS_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <type_traits>
#include <vector>

#include "../scalars/scalar.h"

#define NBMAX 4096  // max size of row block

template <typename T>
std::vector<Index> sort_indexes(const std::vector<T> &v) {
    std::vector<Index> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });
    return idx;
}

template <class InputIt_, class OutputIt_, class Scalar_>
OutputIt_ move_and_scal(InputIt_ first1, InputIt_ last1, OutputIt_ first2, Scalar_ scal) {
    for (; first1 != last1; ++first1, ++first2) {
        *first2 = std::move(*first1);
        *first1 = scal;
    }

    return first2;
}

template <class ForwardIt1_, class ForwardIt2_>
ForwardIt2_ swap_ranges_forward(ForwardIt1_ first1, ForwardIt1_ last1, ForwardIt2_ first2) {
    for (; first1 != last1; ++first1, ++first2) std::swap(*first1, *first2);

    return first2;
}

template <class ForwardIt1_, class ForwardIt2_>
ForwardIt2_ swap_ranges_backward(ForwardIt1_ first1, ForwardIt1_ last1, ForwardIt2_ last2) {
    while (first1 != last1) std::swap(*(--last1), *(--last2));

    return last2;
}

template <class ForwardIt1_, class ForwardIt2_, class Scalar_>
ForwardIt2_ add_ranges_scalar(ForwardIt1_ first1, ForwardIt1_ last1, ForwardIt2_ first2,
                              Scalar_ scal) {
    for (; first1 != last1; ++first1, ++first2) (*first1) += scal * (*first2);

    return first2;
}

template <class ForwardIt_, class Scalar_>
ForwardIt_ scale_range(ForwardIt_ first, ForwardIt_ last, Scalar_ scal) {
    for (; first != last; ++first) (*first) *= scal;

    return first;
}

template <class ForwardIt_, class Scalar_>
ForwardIt_ div_scale_range(ForwardIt_ first, ForwardIt_ last, Scalar_ scal) {
    for (; first != last; ++first) (*first) /= scal;

    return first;
}

template <class InputIt, class T, class UnaryPredicate>
void fill_if(InputIt first, InputIt last, const T &value, UnaryPredicate pred) {
    for (; first != last; ++first)
        if (pred(*first)) *first = value;
}

// x^T*y
void dot(size_t, const Scalar *, const Scalar *, Scalar *);
void dot(const std::vector<Scalar> &, const std::vector<Scalar> &, Scalar *);

// y = a*A*x+b*y
// gemv(m,n,a,A,x,b,y) , A \in R^{m x n}, x \in R^n, y \in R^m
// A in colmajor format
void gemv(Index, Index, Scalar, Scalar *, Scalar *, Scalar, Scalar *);

void norm2(size_t, const Scalar *, Scalar *);

// swap lines in col major matrix
void swap_matrix_lines(Index n, Scalar *A, Index lda, Index k1, Index k2, Index *p);

#endif  // ESOLVER_LINALG_UTILS_H
