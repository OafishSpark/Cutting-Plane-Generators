#include "cuts/cutter.h"

// objective
//     minimize cx
// s.t.:
//     L <= Ax <= U
// bounds:
//     lo <= x <= up
// some parameters in x can be integer

// Compute if_neg_var flag and shift_val value  
// If (if_neg_var == true) then transformation:
//     x --> x* = - x + shift_val
// Back to x:
//     x* --> x = - x* + shift_val
// Else:
//     x --> x* = x - shift_val
// Back to x:
//     x* --> x = x* + shift_val

std::pair<std::vector<bool>, std::vector<Scalar>> Cutter::LinearTransform() {
    std::vector<bool> if_neg_var;
    if_neg_var.resize(a_matrix_scm.GetNRows() + a_matrix_scm.GetNCols(), false);
     
    std::vector<Scalar> shift_val(a_matrix_scm.GetNRows() + a_matrix_scm.GetNCols(), 0.0);

    // variables tranform
    // iv variable ~ iv index in if_neg_var and shift_val
    for (size_t iv = 0; iv < a_matrix_scm.GetNCols(); ++iv) {
        if ((variables_[iv].bnd_lo_ == -kInf) or (Abs(variables_[iv].bnd_up_ - rel_sol_[iv])) < kZeroEpsilon) {
            assert(Abs(variables_[iv].bnd_up_) < kInf);
            if_neg_var[iv] = true;
            // x --> up - x
            // shift_val = up
            shift_val[iv] = variables_[iv].bnd_up_;
        } else {
            // x --> x - lo
            // shift_val = lo
            shift_val[iv] = variables_[iv].bnd_lo_;
        }
    }

    // slack variables transform
    // jv slack variable ~ n_rows + n_cols - jv - 1 index in if_neg_var and shift_val
    for (size_t iv = 0; iv < a_matrix_scm.GetNRows(); ++iv) {
        assert(if_neg_var[if_neg_var.size() - iv - 1] == false);
        if (rhs_[iv].type_ == 'u') {
            // s --> -L - s
            if_neg_var[if_neg_var.size() - iv - 1] = true;
        }
        // s --> s + U
        // either shift_val = -U or shift_val = -L
        shift_val[shift_val.size() - iv - 1] = -rhs_[iv].val_;
    }

    return {if_neg_var, shift_val};
}


std::vector<Cutter::Cut> Cutter::GomoryMixedIntegerCutGenerator() {
    std::vector<Cutter::Cut> cuts;
    
    // compute which variables should be made negative
    std::vector<bool> if_neg_var;
    std::vector<Scalar> shift_val;
    std::tie(if_neg_var, shift_val) = LinearTransform();

    // compute cuts only for basis integer variables
    for (size_t kv = 0; kv < basis_.size(); ++kv) {
        int basis_ind = basis_[kv];
        
        // Check if we can create cut with such a basis variable
        if (basis_ind < 0) {
            // if basis_ind < 0, then it refers to slack variable which is definitely not integer
            continue;
        }
        const auto& basis_var = variables_[basis_ind];
        if (basis_var.is_int_ == false) {
            continue;
        }

        // Transform basis variable to x >= 0
        Scalar b_i;
        if (!if_neg_var[basis_ind]) {
            b_i = rel_sol_[basis_ind] - shift_val[basis_ind];
        } else {
            b_i = shift_val[basis_ind] - rel_sol_[basis_ind];
        }

        // Variable relaxation value must be continuous
        if (Fraq(b_i) < Cutter::kAway) {
            continue;
        }

        // Compute gmi cut cx >= gamma in space related to Ax = b, x >= 0
        Scalar f_0 = b_i - Floor(b_i);
        std::vector<Scalar> cut;
        cut.resize(a_matrix_scm.GetNCols() + a_matrix_scm.GetNRows(), 0);
        Scalar gamma = f_0 * (1 - f_0);
        
        // Compute GMI cut coefficients for original variables
        for (size_t iv = 0; iv < a_matrix_scm.GetNCols(); ++iv) {
            // All coeffs related to basis variables are equal to 0
            if (std::find(basis_.begin(), basis_.end(), iv) != basis_.end()) {
                continue;
            }
            Scalar c_i = 0;
            // Let T be the the simplex table. t_ij --- its element
            Scalar t_ij = SparseVector(b_inv_[kv]) * a_matrix_scm.GetCol(iv);
            if (if_neg_var[iv]) {
                t_ij = -t_ij;
            }
            if (variables_[iv].is_int_) {
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
            if (Abs(c_i) < kZeroEpsilon) {
                continue;
            }
            cut[iv] = c_i;
        }
        
        // Compute GMI cut coefficients for slack variables
        for (size_t iv = 0; iv < a_matrix_scm.GetNRows(); ++iv) {
            // the same movements as in the previous part, but using another tranformation
            if (std::find(basis_.begin(), basis_.end(), -iv - 1) != basis_.end()) {
                continue;
            }
            Scalar c_i = 0;
            Scalar t_ij = b_inv_[kv][iv];
            if (if_neg_var[if_neg_var.size() - iv - 1]) {
                t_ij = -t_ij;
            }
            if (t_ij >= 0) {
                c_i = t_ij * (1 - f_0);
            } else {
                c_i = -t_ij * f_0;
            }
            if (Abs(c_i) < kZeroEpsilon) {
                continue;
            }
            cut[cut.size() - iv - 1] = c_i;
        }

        // Transform back
        for (size_t iv = 0; iv < cut.size(); ++iv) {
            if (Abs(cut[iv]) < kZeroEpsilon) {
                continue;
            }
            if (if_neg_var.at(iv)) {
                gamma -= cut[iv] * shift_val[iv];
                cut[iv] = -cut[iv];
            } else {
                gamma += cut[iv] * shift_val[iv];
            }
        }

        // Remove slack vars
        std::vector<Scalar> slack_cut(a_matrix_scm.GetNRows(), 0);
        for (size_t iv = 0; iv < a_matrix_scm.GetNRows(); ++iv) {
            slack_cut[iv] = cut[cut.size() - iv - 1];
        }
        cut.resize(a_matrix_scm.GetNCols());
        SparseVector sparse_slack_cut = SparseVector(slack_cut);
        for (size_t iv = 0; iv < a_matrix_scm.GetNCols(); ++iv) {
            cut[iv] -= sparse_slack_cut * a_matrix_scm.GetCol(iv);
        }
        
        int nonzero_num{0};
        for (int iv{0}; iv < cut.size(); ++iv) {
            if (Abs(cut[iv]) <= kZeroEpsilon) {
                cut[iv] = Scalar(0);
                continue;                
            }
            ++nonzero_num;
        }
        assert(nonzero_num > 0);
        
        SparseVector lhv(cut);
        
        // Estimate it by finding Euclidean distance to relaxation solution project onto antigradient
        Scalar violation = (SparseVector(rel_sol_) * lhv - gamma);
        if (violation > Scalar(0)) {
            assert(false);
        }
        //Scalar estimation = (violation * violation) / (lhv * lhv);
        Scalar estimation = (((lhv * violation) * objective_) * ((lhv * violation) * objective_)) / (lhv * lhv) / (objective_ * objective_); 
        
        // Add cut to temporary pool
        cuts.push_back(Cutter::Cut{lhv, gamma, estimation, false});
    }
    return cuts;
}
