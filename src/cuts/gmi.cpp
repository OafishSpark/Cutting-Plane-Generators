#include "../../headers/cuts/cutter.h"


std::vector<Cutter::Cut> Cutter::GomoryMixedIntegerCutGenerator() {
     std::vector<Cutter::Cut> cuts;
    // compute which variables should be made negative
    std::vector<bool> if_neg_var;
    if_neg_var.resize(a_matrix_.rows_ + a_matrix_.cols_, false);
    // l / u -- shift_val: x --> x - l or x --> u - x
    std::vector<Scalar> shift_val(a_matrix_.rows_ + a_matrix_.cols_, 0.0);
    for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
        if ((vars_[iv].bnd_lo_ == -kInf) or (Abs(vars_[iv].bnd_up_ - sol_[iv])) < kZero_Epsilon) {
            assert(Abs(vars_[iv].bnd_up_) < kInf);
            if_neg_var[iv] = true;
            shift_val[iv] = vars_[iv].bnd_up_;
        } else {
            shift_val[iv] = vars_[iv].bnd_lo_;
        }
    }
    for (int iv = 0; iv < a_matrix_.rows_; ++iv) {
        assert(if_neg_var[if_neg_var.size()-iv-1] == false);
        if (rhs_[iv].type_ == 'U') {
            if_neg_var[if_neg_var.size()-iv-1] = true;
        } else {
            shift_val[shift_val.size()-iv-1] = -rhs_[iv].val_;
        }
    }
    // compute cuts only for basis integer variables
    for (const auto& basis_ind: basis_) {
        if (basis_ind < 0) {
            // if basis_ind < 0, then it refers to slack variable which is definitely not integer
            continue;
        }
        const auto& basis_var = vars_[basis_ind];
        if (basis_var.is_int_ == false) {
            continue;
        }
        // Transform to Ax = b, x >= 0
        Scalar b_i;
        if (!if_neg_var[basis_ind]) {
            b_i = sol_[basis_ind] - shift_val[basis_ind];
        } else {
            b_i = shift_val[basis_ind] - sol_[basis_ind];
        }
        // Variable relaxation value must be continuous
        if (Fraq(b_i) < Cutter::kAway) {
            continue;
        }
        // Compute gmi cut cx >= b in space related to Ax = b, x >= 0
        Scalar f_0 = b_i - Floor(b_i);
        std::vector<Scalar> cut(a_matrix_.cols_ + a_matrix_.rows_);
        cut.resize(a_matrix_.cols_ + a_matrix_.rows_, 0);
        Scalar gamma = f_0 * (1 - f_0);
        // Compute GMI cut coefficients for original variables
        for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
            if (std::find(basis_.begin(), basis_.end(), iv) != basis_.end()) {
                continue;
            }
            Scalar c_i = 0;
            Scalar t_ij = b_inv_[basis_ind] * a_matrix_[iv];
            if (if_neg_var[iv]) {
                t_ij = -t_ij;
            }
            if (vars_[iv].is_int_) {
                Scalar f_j = t_ij - floor(t_ij);
                if (f_j <= f_0) {
                    c_i = f_j * (1 - f_0);
                } else {
                    c_i = (1 - f_j) * f_0;
                }
            } else {
                if (t_ij >= 0) {
                    c_i = t_ij * (1 - f_0);
                } else {
                    c_i = -t_ij * f_0;
                }
            }
            if (Abs(c_i) < kZero_Epsilon) {
                continue;
            }
            cut[iv] = c_i;
        }
        // Compute GMI cut coefficients for slack variables
        for (int iv = 0; iv < a_matrix_.rows_; ++iv) {
            if (std::find(basis_.begin(), basis_.end(), -iv-1) != basis_.end()) {
                continue;
            }
            Scalar c_i = 0;
            Scalar t_ij = b_inv_[basis_ind][iv];
            if (if_neg_var[if_neg_var.size() - iv - 1]) {
                t_ij = -t_ij;
            }
            if (t_ij >= 0) {
                c_i = t_ij * (1 - f_0);
            } else {
                c_i = -t_ij * f_0;
            }
            if (Abs(c_i) < kZero_Epsilon) {
                continue;
            }
            cut[cut.size() - iv - 1] = c_i;
        }
        // Transform back
        for (int iv = 0; iv < cut.size(); ++iv) {
            if (Abs(cut[iv]) < kZero_Epsilon) {
                continue;
            }
            if (if_neg_var.at(iv)) {
                gamma += -cut[iv] * shift_val[iv];
                cut[iv] = -cut[iv];
            } else {
                gamma += cut[iv] * shift_val[iv];
            }
        }
        // Remove slack vars
        std::vector<Scalar> slack_cut(a_matrix_.rows_, 0);
        for (int iv = 0; iv < a_matrix_.rows_; ++iv) {
            slack_cut[iv] = cut[cut.size() - iv - 1];
        } 
        cut.resize(a_matrix_.cols_);
        for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
            cut[iv] -= slack_cut * a_matrix_[iv];
        }
        // Estimate it by finding Euclidean distance to relaxation solution
        SparseVector lhv(cut); 
        Scalar violation = sol_ * lhv - gamma;
        // if (violation > Scalar(0)) {
        //     continue;
        // }
        Scalar estimation = (violation * violation) / (lhv * lhv);

        // Add cut to temporary pool
        cuts.push_back(Cutter::Cut{lhv, gamma, estimation, false});
    }
    return cuts; 
}
