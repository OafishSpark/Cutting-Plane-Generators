#include "../../include/linalg/lup.h"

//#include <cblas.h>
//#include <lapacke.h>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <numeric>

#include "linalg_utils.h"
#include "raw/params.h"

Index LUP::unblocked_lu(Index m, Index n, Scalar *A, Index lda, Index *p) {
    Scalar abs_value;
    Scalar *Acol = A;

    for (Index j = 0; j < n; ++j, ++p, Acol += lda) {
        Scalar piv_value = abs(*(Acol + j));
        *p = j;
        // find the pivot
        for (Index k = j + 1; k < m; ++k) {
            if ((abs_value = abs(*(Acol + k))) > piv_value) {
                piv_value = abs_value;
                *p = k;
            }
        }

        if (piv_value != Scalar(0)) {
            // swap rows of matrix A
            Scalar *Acol2 = A;
            if (*p > j) {
                for (; Acol2 != A + n * lda; Acol2 += lda) std::swap(*(Acol2 + j), *(Acol2 + *p));
            }

            // scal rows of A
            Index l = j + 1;
            scale_range(Acol + l, Acol + m, Scalar(1) / *(Acol + j));
            Acol2 = Acol + lda;
            for (Index k = l; k < n; ++k, Acol2 += lda) {
                Scalar alpha = *(Acol2 + j);
                add_ranges_scalar(Acol2 + l, Acol2 + m, Acol + l, -alpha);
            }

        } else {
            return j;
        }
    }

    return 0;
}

Index LUP::blocked_lu(Index m, Index n, Scalar *A, Index lda, Index *p) {
    if (!blocked_ || (n <= split_size)) {
        // Unblocked LU
        return unblocked_lu(m, n, A, lda, p);
    }

    // Splitting
    const Index n1 = matrix_split(n);
    const Index n2 = n - n1;
    const Index m2 = m - n1;

    // split A = [ A_11 A_12 ]
    //           [ A_21 A_22 ]
    Scalar *const A_11 = A;
    Scalar *const A_12 = A + n1 * lda;
    Scalar *const A_21 = A + n1;
    Scalar *const A_22 = A + n1 * lda + n1;

    // split p = [ p_1 ]
    //           [ p_2 ]
    Index *const p_1 = p;
    Index *const p_2 = p + n1;

    // apply lu to [A_11]
    //             [A_21]
    Index info = blocked_lu(m, n1, A_11, lda, p_1);
    if (info != 0) return info;
    // apply row permutation to A_12
    swap_matrix_lines(n2, A_12, lda, 0, n1, p_1);

    // solve lower triangular sistem: A_12 = A_11 \ A_12
    // solve_lu_l(n1, n2, A_11, lda, A_12, lda);
    // and
    // update: A_22 = A_22 - A_21 * A_12
    for (Index j = 0; j < n2; ++j) {
        Index jcol = j * lda;
        Scalar *A_11col = A_11;
        Scalar *A_22col = A_22 + jcol;
        Scalar *A_12col = A_12 + jcol;
        for (Index k = 0; k < n1; ++k, A_11col += lda) {
            Scalar *A_12el = A_12col + k;
            if (*(A_12el) != Scalar(0)) {
                // A_12 = A_11 \ A_12
                add_ranges_scalar(A_12col + k + 1, A_12col + n1, A_11col + k + 1, -(*A_12el));
                // A_22 = A_22 - A_21 * A_12
                Scalar *A_21col = A_21 + k * lda;
                for (Index i = 0; i < m2; ++i) *(A_22col + i) -= *(A_21col + i) * *(A_12el);
            }
        }
    }

    // apply lu to A_22
    info = blocked_lu(m2, n2, A_22, lda, p_2);
    if (info != 0) info += n1;
    // apply row permutation to A_21
    swap_matrix_lines(n1, A_21, lda, 0, n2, p_2);
    // shift pivots
    for (Index i = 0; i < n2; ++i) *(p_2 + i) += n1;

    return info;
}

Index LUP::unblocked_inv_u(Index n, Scalar *A, Index lda) {
    // inverse in place of upper triangular matrix
    Scalar *Acol = A;
    for (Index j = 0; j < n; ++j, Acol += lda) {
        *(Acol + j) = Scalar(1) / *(Acol + j);
        Scalar *A_k_col = A;
        for (Index k = 0; k < j; ++k, A_k_col += lda) {
            Scalar *Ael = Acol + k;
            if (*Ael != Scalar(0)) {
                add_ranges_scalar(Acol, Acol + k, A_k_col, *Ael);
                *Ael *= *(A + k * lda + k);
            }
        }
        scale_range(Acol, Acol + j, -*(Acol + j));
    }

    return 0;
}

Index LUP::blocked_inv_u(Index n, Scalar *A, Index lda) {
    if (!blocked_ || (n <= split_size)) {
        // Unblocked inversion
        return unblocked_inv_u(n, A, lda);
    }

    // Splitting
    const Index n1 = matrix_split(n);
    const Index n2 = n - n1;

    // split U = [ U_11 U_12 ]
    //           [   0  U_22 ]
    Scalar *const U_11 = A;
    Scalar *const U_12 = A + n1 * lda;
    Scalar *const U_22 = A + n1 * lda + n1;

    // apply inversion to U_11
    Index info = blocked_inv_u(n1, U_11, lda);
    if (info != 0) return info;

    // update: U_12 = - U_11 * U_12
    // and
    // solve upper triangular sistem: U_12  = U_12 * U_22^{-1}
    //
    // Equivalent BLAS functions:
    // cblas_dtrmm(CblasColMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,n1,n2,-1.0,U_11,lda,U_12,lda);
    // cblas_dtrsm(CblasColMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,n1,n2,1.0,U_22,lda,U_12,lda);
    //
    for (Index j = 0; j < n2; ++j) {
        Index jcol = j * lda;
        Scalar *U_11_k_col = U_11;
        Scalar *U_12_j_col = U_12 + jcol;
        // U_12 = - U_11 * U_12
        for (Index k = 0; k < n1; ++k, U_11_k_col += lda) {
            Scalar *U_12el = U_12_j_col + k;
            if (*U_12el != Scalar(0)) {
                *U_12el *= Scalar(-1);
                add_ranges_scalar(U_12_j_col, U_12_j_col + k, U_11_k_col, *U_12el);
                *U_12el *= *(U_11_k_col + k);
            }
        }

        // U_12 = U_12 * U_22^{-1}
        Scalar *U_22_j_col = U_22 + jcol;
        Scalar *U_12_k_col = U_12;
        for (Index k = 0; k < j; ++k, U_12_k_col += lda) {
            Scalar *U_22_kj_el = U_22_j_col + k;
            if (*U_22_kj_el != Scalar(0)) {
                add_ranges_scalar(U_12_j_col, U_12_j_col + n1, U_12_k_col, -*(U_22_kj_el));
            }
        }
        if (*(U_22_j_col + j) != Scalar(1)) {
            Scalar alpha = Scalar(1) / *(U_22_j_col + j);
            scale_range(U_12_j_col, U_12_j_col + n1, alpha);
        }
    }

    // apply inversion to U_22
    info = blocked_inv_u(n2, U_22, lda);
    if (info != 0) info += n1;

    return info;
}

void LUP::solve_lu_l(Index m, Index n, const Scalar *A, Index lda, Scalar *B, Index ldb) {
    if (n == 0) return;
    for (Scalar *Bcol = B; Bcol != B + n * ldb; Bcol += ldb) {
        const Scalar *Acol = A;
        for (Index k = 0; k < m; ++k, Acol += lda) {
            Scalar *alpha = (Bcol + k);
            if (*alpha != Scalar(0)) {
                add_ranges_scalar(Bcol + k + 1, Bcol + m, Acol + k + 1, -(*alpha));
            }
        }
    }
}

void LUP::solve_lu_u(Index m, Index n, const Scalar *A, Index lda, Scalar *B, Index ldb) {
    if (n == 0) return;
    for (Scalar *Bcol = B; Bcol != B + n * ldb; Bcol += ldb) {
        const Scalar *Acol = A + (m - 1) * lda;
        for (Index k = m - 1; k >= 0; --k, Acol -= lda) {
            Scalar *alpha = Bcol + k;
            if (*alpha != Scalar(0)) {
                *alpha /= *(Acol + k);
                add_ranges_scalar(Bcol + k + 1, Bcol + m, Acol + k + 1, -(*alpha));
            }
        }
    }
}

Index LUP::lu(Index m, Index n, Scalar *A, Index lda, Index *p) {
    Index info = 0;
    // Check arguments
    if (m < 0)
        info = -1;
    else if (n < 0)
        info = -2;
    else if (lda < std::max(Index(1), m))
        info = -4;
    if (info != 0) return info;

    if (m == 0 || n == 0) return info;

    const Index sn = m <= n ? m : n;
    info = blocked_lu(m, sn, A, lda, p);

    // Right remainder
    if (m < n) {
        // get reminder cols number
        const Index rn = n - m;

        // split matrix A = [A_L A_R]
        const Scalar *const A_L = A;
        Scalar *const A_R = A + lda * m;

        // apply row permutation to A_R
        swap_matrix_lines(rn, A_R, lda, 0, m, p);
        // solve lower triangular sistem: A_R = A_L \ A_R
        solve_lu_l(m, rn, A_L, lda, A_R, lda);
    }

    // Remove almost zero elems:
    // for (Scalar *Acol = A; Acol != A + n * lda; Acol += lda)
    //     fill_if(Acol, Acol + m, Scalar(0),
    //             [](Scalar x) { return abs(x) < Params::GetZeroEps() * Params::GetZeroEps(); });

    return info;
}

Index LUP::inverse(Index n, Scalar *A, Index lda, const Index *p) {
    Index info = 0;
    // Check arguments
    if (n < 0)
        info = -1;
    else if (lda < std::max(Index(1), n))
        info = -3;
    if (info != 0) return info;

    if (n == 0) return info;

    info = blocked_inv_u(n, A, lda);
    if (info > 0) return info;

    // Solve the equation A^{-1}*L = U^{-1} for A^{-1}
    if (!blocked_ || (n <= split_size)) {
        // unblocked
        work_ = new Scalar[n];
        Scalar *Acol = A + (n - 2) * lda;

        for (Index j = n - 2; j >= 0; --j, Acol -= lda) {
            std::copy(Acol + j + 1, Acol + n, work_ + j + 1);
            std::fill(Acol + j + 1, Acol + n, Scalar(0));

            Scalar *A_k_col = Acol + lda;
            Scalar *work_k = work_ + j + 1;
            for (; work_k != work_ + n; A_k_col += lda, ++work_k) {
                if (*(work_k) != Scalar(0)) add_ranges_scalar(Acol, Acol + n, A_k_col, -(*work_k));
            }
        }
        delete[] work_;
    } else {
        // blocked
        work_ = new Scalar[split_size * n];

        Index nn = (((n - 1) >> split_step) << split_step);
        Index jb = n - nn;

        // solve for first block
        Scalar *Acol = A + nn * lda;
        Scalar *Wcol = work_;
        for (Index jj = nn; jj < nn + jb; ++jj, Acol += lda, Wcol += n) {
            std::copy(Acol + jj + 1, Acol + n, Wcol + jj + 1);
            std::fill(Acol + jj + 1, Acol + n, Scalar(0));
        }
        // Equivalent BLAS function:
        // cblas_dtrsm(CblasColMajor,CblasRight,CblasLower,CblasNoTrans,CblasUnit, n, jb, 1.0,
        // W+nn,n, A+nn*lda, lda);
        Acol = A + (n - 1) * lda;
        Wcol = work_ + (jb - 1) * n + nn;
        for (Index jj = jb - 1; jj >= 0; --jj, Acol -= lda, Wcol -= n) {
            Scalar *Acol_k = Acol + lda;
            for (Scalar *Wel = Wcol + jj + 1; Wel != Wcol + jb; ++Wel, Acol_k += lda) {
                if (*(Wel) != Scalar(0)) add_ranges_scalar(Acol, Acol + n, Acol_k, -(*Wel));
            }
        }

        for (Index j = nn - split_size; j >= 0; j -= split_size) {
            Acol = A + j * lda;
            Wcol = work_;
            for (Index jj = j; jj < j + split_size; ++jj, Acol += lda, Wcol += n) {
                std::copy(Acol + jj + 1, Acol + n, Wcol + jj + 1);
                std::fill(Acol + jj + 1, Acol + n, Scalar(0));
            }

            // Equivalent BLAS function:
            // cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans, n, split_size,
            // n-j-split_size,-1.0, A+(j+split_size)*lda, lda, W+j+split_size, n, 1.0, A+j*lda,
            // lda);
            Wcol = work_ + j + split_size;
            Acol = A + j * lda;
            for (; Acol != A + (j + split_size) * lda; Acol += lda, Wcol += n) {
                Scalar *Acol_k = A + (j + split_size) * lda;
                for (Scalar *Wel = Wcol; Wel != Wcol + n - j - split_size; ++Wel, Acol_k += lda) {
                    if (*Wel != Scalar(0)) add_ranges_scalar(Acol, Acol + n, Acol_k, -(*Wel));
                }
            }

            // Equivalent BLAS function:
            // cblas_dtrsm(CblasColMajor,CblasRight,CblasLower,CblasNoTrans,CblasUnit, n,
            // split_size, 1.0, W+j, n, A+j*lda, lda);
            Acol = A + (j + split_size - 1) * lda;
            Wcol = work_ + (split_size - 1) * n + j;
            for (Index jj = split_size - 1; jj >= 0; --jj, Acol -= lda, Wcol -= n) {
                Scalar *A_k_col = Acol + lda;
                for (Scalar *Wel = Wcol + jj + 1; Wel != Wcol + split_size; ++Wel, A_k_col += lda) {
                    if (*(Wel) != Scalar(0)) add_ranges_scalar(Acol, Acol + n, A_k_col, -(*Wel));
                }
            }
        }
        delete[] work_;
    }
    // permute columns
    for (Index j = n - 2; j >= 0; --j) {
        if (*(p + j) != j) {
            swap_ranges_forward(A + j * lda, A + (j + 1) * n, A + lda * *(p + j));
        }
    }
    // Remove almost zero elems:
    // for (Scalar *Acol = A; Acol != A + n * lda; Acol += lda)
    //     fill_if(Acol, Acol + n, Scalar(0),
    //             [](Scalar x) { return abs(x) < kEpsInteger*kEpsInteger; });
    return info;
}

