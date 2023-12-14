#include <cassert>

#include "../../include/linalg/linalg_utils.h"
#include "raw/params.h"

void dot(size_t n, const Scalar *x, const Scalar *y, Scalar *result) {
    Scalar dot_0, dot_1, dot_2, dot_3;
    dot_0 = dot_1 = dot_2 = dot_3 = Scalar(0);

    size_t n1 = (n & Index(-4 * 8));
    for (Index k = 0; k < Index(n1); k += 4 * 8) {
        dot_0 += x[k] * y[k] + x[k + 1] * y[k + 1] + x[k + 2] * y[k + 2] + x[k + 3] * y[k + 3] +
                 x[k + 4] * y[k + 4] + x[k + 5] * y[k + 5] + x[k + 6] * y[k + 6] +
                 x[k + 7] * y[k + 7];
        dot_1 += x[k + 8] * y[k + 8] + x[k + 9] * y[k + 9] + x[k + 10] * y[k + 10] +
                 x[k + 11] * y[k + 11] + x[k + 12] * y[k + 12] + x[k + 13] * y[k + 13] +
                 x[k + 14] * y[k + 14] + x[k + 15] * y[k + 15];
        dot_2 += x[k + 16] * y[k + 16] + x[k + 17] * y[k + 17] + x[k + 18] * y[k + 18] +
                 x[k + 19] * y[k + 19] + x[k + 20] * y[k + 20] + x[k + 21] * y[k + 21] +
                 x[k + 22] * y[k + 22] + x[k + 23] * y[k + 23];
        dot_3 += x[k + 24] * y[k + 24] + x[k + 25] * y[k + 25] + x[k + 26] * y[k + 26] +
                 x[k + 27] * y[k + 27] + x[k + 28] * y[k + 28] + x[k + 29] * y[k + 29] +
                 x[k + 30] * y[k + 30] + x[k + 31] * y[k + 31];
    }
    size_t n2 = (n & Index(-8));
    for (Index k(n1); k < Index(n2); k += 8) {
        dot_0 += x[k] * y[k] + x[k + 1] * y[k + 1] + x[k + 2] * y[k + 2] + x[k + 3] * y[k + 3] +
                 x[k + 4] * y[k + 4] + x[k + 5] * y[k + 5] + x[k + 6] * y[k + 6] +
                 x[k + 7] * y[k + 7];
    }
    *result = dot_0 + dot_1 + dot_2 + dot_3;
    for (Index k(n2); k < Index(n); ++k) *result += x[k] * y[k];
}

void dot(const std::vector<Scalar> &x, const std::vector<Scalar> &y, Scalar *result) {
    assert(x.size() == y.size());
    return dot(x.size(), x.data(), y.data(), result);
}

void norm2(size_t n, const Scalar *x, Scalar *result) { dot(n, x, x, result); }

static void gemv_kernel_4x4(Index n, Scalar **ac, const Scalar *x, Scalar *y) {
    Scalar *a0 = ac[0];
    Scalar *a1 = ac[1];
    Scalar *a2 = ac[2];
    Scalar *a3 = ac[3];

    for (Index i = 0; i < n; i += 4) {
        y[i] += a0[i] * x[0] + a1[i] * x[1] + a2[i] * x[2] + a3[i] * x[3];
        y[i + 1] += a0[i + 1] * x[0] + a1[i + 1] * x[1] + a2[i + 1] * x[2] + a3[i + 1] * x[3];
        y[i + 2] += a0[i + 2] * x[0] + a1[i + 2] * x[1] + a2[i + 2] * x[2] + a3[i + 2] * x[3];
        y[i + 3] += a0[i + 3] * x[0] + a1[i + 3] * x[1] + a2[i + 3] * x[2] + a3[i + 3] * x[3];
    }
}

static void gemv_kernel_4x1(Index n, const Scalar *ac, const Scalar *x, Scalar *y) {
    for (Index i = 0; i < n; i += 4) {
        y[i] += ac[i] * x[0];
        y[i + 1] += ac[i + 1] * x[0];
        y[i + 2] += ac[i + 2] * x[0];
        y[i + 3] += ac[i + 3] * x[0];
    }
}

static void scale_y_4(Index n, Scalar b, Scalar *y) {
    for (Index i = 0; i < n; i += 4) {
        y[i] *= b;
        y[i + 1] *= b;
        y[i + 2] *= b;
        y[i + 3] *= b;
    }
}

void gemv(Index m, Index n, Scalar a, Scalar *A, Scalar *x, Scalar b, Scalar *y) {
    // n/4 - number of 4 cols blocks in A
    Index n1 = n >> 2;
    // reminder block with cols < 4
    Index n2 = n % Index(4);

    // max number of rows mod 16
    Index m1 = m - (m % Index(16));
    // max number of rows mod 16 in reminder mod NBMAX (>=1024)
    Index m2 = (m % NBMAX) - (m % Index(16));
    // A cols pointers
    Scalar *ac[4];
    // x buffer
    Scalar x_b[4];
    // pointer to *x
    Scalar *x_p;
    // pointer to *A
    Scalar *A_p;

    auto nb = Index(NBMAX);
    while (nb == Index(NBMAX)) {
        m1 -= nb;
        if (m1 < Index(0)) {
            if (m2 == Index(0)) break;
            nb = m2;
        }

        scale_y_4(nb, b, y);

        x_p = x;
        A_p = A;

        for (Index i = 0; i < n1; ++i) {
            // get 4 cols of A
            ac[0] = A_p;
            ac[1] = ac[0] + m;
            ac[2] = ac[1] + m;
            ac[3] = ac[2] + m;

            x_b[0] = a * (*(x_p++));
            x_b[1] = a * (*(x_p++));
            x_b[2] = a * (*(x_p++));
            x_b[3] = a * (*(x_p++));
            gemv_kernel_4x4(nb, ac, x_b, y);
            A_p += m << 2;
        }

        for (Index i = 0; i < n2; ++i) {
            x_b[0] = a * (*(x_p++));
            gemv_kernel_4x1(nb, A_p, x_b, y);
            A_p += m;
        }
        A += nb;
        y += nb;
    }

    Index j = 0;
    while (j < (m % Index(16))) {
        A_p = A;
        x_p = x;
        Scalar tmp(0);
        for (Index i = 0; i < n; ++i) {
            tmp += (*A_p) * (*x_p);
            A_p += m;
            x_p++;
        }
        *y *= b;
        *y += a * tmp;
        y++;
        A++;
        j++;
    }
}

// swap lines in col major matrix
void swap_matrix_lines(Index n, Scalar *A, Index lda, Index k1, Index k2, Index *p) {
    const Index n16 = (n >> split_step) << split_step;
    if (n16 != 0)
        for (Index j = 0; j < n16; j += split_size) {
            for (Index i = k1; i < k2; ++i) {
                Index *ip = p + i;
                if (*ip != i)
                    for (Scalar *Acol = A + j * lda; Acol != A + (j + split_size) * lda;
                         Acol += lda) {
                        std::swap(*(Acol + i), *(Acol + *ip));
                    }
            }
        }
    if (n16 != n)
        for (Index i = k1; i < k2; ++i) {
            Index *ip = p + i;
            if (*ip != i)
                for (Scalar *Acol = A + n16 * lda; Acol != A + n * lda; Acol += lda)
                    std::swap(*(Acol + i), *(Acol + *ip));
        }
}
