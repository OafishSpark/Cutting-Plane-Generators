#include "cuts/cutter.h"
#include "linalg/sparse_vector.h"
#include "scalars/scalar.h"
#include "utils/utils.h"


struct VarBnd {
    // See Achterberg constraint programming
    // related to j-th variable
    // x_j <= (>=) s * x_k + d, x_k is integer
    int k;
    Scalar s;
    Scalar d;
    // Lower 'l' or upper 'u'
    char mode;
    // | x_j - (s * s_k + d) |
    Scalar value;
};


struct Roww{
    // a_matrix_row
    SparseVector row;
    Scalar rhs;
    // violation of the row
    Scalar score;
};

bool CompareRowws(Roww& r1, Roww& r2) {
    return r1.score < r2.score;
}

inline Scalar Max(Scalar a, Scalar b) {
    return a > b ? a : b;
}

inline Scalar Min(Scalar a, Scalar b) {
    return a < b ? a : b;
}

std::vector<Cutter::Cut> Cutter::MixedIntegerRoundingCutGenerator() {
    std::vector<Cut> cuts;
    // Compute columns for aggregation and variable bounds
    std::vector<Roww> rows(n_rows_, Roww({SparseVector(), 0, 10}));
    std::vector<VarBnd> var_bnds;
    var_bnds.resize(n_cols_, {-1, 0, 0, 'l', -1});
    Index nv = 0;
    for (int iv{0}; iv < n_rows_; ++iv) {
        SparseVector row = a_matrix_srm_.GetRow(iv);
        // Check if row is variable bound
        if (row.index.size() == 2) {
            for (int jv{0}; jv < 2; ++jv) {
                if (!variables_[row.index[jv]].is_int_ || variables_[row.index[int(jv == 0)]].is_int_)
                    continue;
                Index fst = row.index[jv];
                Index scnd_jv = int(jv == 0);
                Index scnd = row.index[scnd_jv];
                Scalar d = rhs_[iv].val_ / row.values[scnd_jv];
                Scalar s = - row.values[jv] / row.values[scnd_jv];
                char mode = 'l';
                if (row.values[scnd_jv] < 0) {
                    if (rhs_[iv].type_ == 'l')
                        mode = 'u';
                } else {
                    mode = rhs_[iv].type_;
                }
                Scalar value = Abs(rel_sol_[scnd] - (s * rel_sol_[fst] + d));
                if (var_bnds[scnd].value == -1 or var_bnds[fst].value > value)
                    var_bnds[scnd] = VarBnd({fst, s, d, mode, value});
            }
        }
        // Density check
        Scalar density = Scalar(row.index.size()) / n_cols_;
        if (density > kMaxDensity)
            continue;
        // Slack check
        Scalar slack = row * SparseVector(rel_sol_) - rhs_[iv].val_;
        if (Abs(slack) > kMaxSlack)
            continue;
        // Compute score
        Scalar score = density + slack;
        // Add to row set
        Scalar rhs = rhs_[iv].val_;
        if (rhs_[iv].type_ == 'l') {
            rhs *= -1;
            row *= -1;
        }
        rows[nv++] = {row, rhs, score};
    }
    rows.resize(nv);
    std::sort(rows.begin(), rows.end(), CompareRowws);
    std::vector<Scalar> aggregaion_list(n_cols_, 0.0);
    for (int iv{0}; iv < n_cols_; ++iv)
        if (!variables_[iv].is_int_)
            aggregaion_list[iv] = Min(variables_[iv].bnd_up_ - rel_sol_[iv], rel_sol_[iv] - variables_[iv].bnd_lo_);
    // Calculating cuts
    for (int iv{0}; iv < nv; ++iv) {
        SparseVector row = rows[iv].row;
        Scalar rhs = rows[iv].rhs;
        // mir_cut * x >= gamma
        SparseVector mir_cut;
        Scalar gamma = 0;
        Scalar violation = 1;
        for (int jv{0}; jv < 6; ++jv) {
            std::vector<Scalar> dense_temp_row(n_cols_, 0.0);
            Scalar temp_b = rhs;
            std::vector<bool> if_neg_var(n_cols_, false);
            std::vector<Scalar> if_var_bnd(n_cols_, false);
            // Transform to x >= 0
            for (int kv{0}; kv < row.index.size(); ++kv) {
                int ind = row.index[kv];
                Scalar score_lo = rel_sol_[ind] - variables_[ind].bnd_lo_;
                Scalar score_up = variables_[ind].bnd_up_ - rel_sol_[ind];
                Scalar score_var = var_bnds[ind].value;
                if (score_var == -1)
                    score_var = kInf;
                if (score_var <= score_lo && score_var <= score_up) {
                    if_var_bnd[ind] = true;
                    if (var_bnds[ind].mode == 'u') {
                        dense_temp_row[ind] -= row.values[kv];
                    } else
                        dense_temp_row[ind] += row.values[kv];
                    dense_temp_row[var_bnds[ind].k] += row.values[kv] * var_bnds[ind].s;
                    temp_b -= row.values[kv] * var_bnds[ind].d;
                    continue;
                }
                if (score_up < score_lo) {
                    dense_temp_row[ind] -= row.values[kv];
                    if_neg_var[ind] = true;
                    temp_b -= row.values[kv] * variables_[ind].bnd_up_;
                } else {
                    dense_temp_row[ind] += row.values[kv];
                    temp_b -= row.values[kv] * variables_[ind].bnd_lo_;
                }
            }
            Scalar max_int_mod = 0.0;
            for (int kv{0}; kv < n_cols_; ++kv)
                if (variables_[kv].is_int_)
                    if (dense_temp_row[kv] > max_int_mod)
                        max_int_mod = dense_temp_row[kv];

            std::vector<Scalar> deltas = {1, 2, 1.0 / 2, 1.0 / 4, max_int_mod, max_int_mod / 2, max_int_mod / 4};
            for (auto delta : deltas) {
                // Build MIR cut
                std::vector<Scalar> dense_cut(dense_temp_row);
                for (auto& elem: dense_cut)
                    elem /= delta;
                Scalar temp_gamma = temp_b / delta;
                bool if_leq = (delta < 0);
                Scalar f_0 = temp_gamma - Floor(temp_gamma);
                if (Abs(f_0) < 0.99) {
                    for (int kv{0}; kv < n_cols_; ++kv) {
                        Scalar f_j = dense_cut[kv] - Floor(dense_cut[kv]);
                        if (variables_[kv].is_int_)
                            dense_cut[kv] = Floor(dense_cut[kv]) + Max(f_j - f_0, 0) / (1 - f_0);
                        else
                            dense_cut[kv] = Max(-dense_cut[kv], 0) / (1 - f_0);
                    }
                    temp_gamma = Floor(temp_gamma);
                    // Transform back
                    for (int kv{0}; kv < n_cols_; ++kv) {
                        if (if_var_bnd[kv])
                            continue;
                        if (Abs(dense_cut[kv]) < kZeroEpsilon) {
                            continue;
                        }
                        if (if_neg_var.at(kv)) {
                            temp_gamma -= dense_cut[kv] * variables_[kv].bnd_up_;
                            dense_cut[kv] = -dense_cut[kv];
                        } else {
                            temp_gamma += dense_cut[kv] * variables_[kv].bnd_lo_;
                        }
                    }
                    for (int kv{0}; kv < n_cols_; ++kv) {
                        if (!if_var_bnd[kv])
                            continue;
                        if (var_bnds[kv].mode == 'u') {
                            temp_gamma -= dense_cut[kv] * var_bnds[kv].d;
                            dense_cut[var_bnds[kv].k] += dense_cut[kv] * var_bnds[kv].s;
                            dense_cut[kv] = - dense_cut[kv];
                        } else {
                            temp_gamma -= dense_cut[kv] * var_bnds[kv].d;
                            dense_cut[var_bnds[kv].k] += dense_cut[kv] * var_bnds[kv].s;
                        }
                    }
                    // Compare to existing cut
                    SparseVector lhs = SparseVector(dense_cut);
                    if (if_leq) {
                        lhs *= -1;
                        temp_gamma *= -1;
                    }
                    Scalar temp_violation = lhs * SparseVector(rel_sol_) - temp_gamma;
                    if (temp_violation >= 0)
                        continue;
                    if (temp_violation > violation)
                        continue;
                    mir_cut = lhs;
                    violation = temp_violation;
                    gamma = temp_gamma;
                }
            }
            // Aggregate to other row
            bool done = false;
            int max_iter = 0;
            for (int hv{0}; hv < (rows.size() > max_iter ? max_iter : rows.size()); ++hv) {
                if (hv == iv)
                    continue;
                SparseVector row1 = rows[hv].row;
                Scalar b1 = rows[hv].rhs;
                for (int kv{0}; kv < n_cols_; ++kv) {
                    if (aggregaion_list[kv] == 0)
                        continue;
                    auto orig_iter = std::find(row.index.begin(), row.index.end(), kv);
                    if (orig_iter == row.index.end())
                        continue;
                    auto other_iter = std::find(row1.index.begin(), row1.index.end(), kv);
                    if (other_iter == row1.index.end())
                        continue;
                    int orig_ind = std::distance(row.index.begin(), orig_iter);
                    int other_ind = std::distance(row1.index.begin(), other_iter);
                    if (row.values[orig_ind] * row.values[other_ind] >= 0)
                        continue;
                    Scalar temp = Abs(row.values[orig_ind]) / Abs(row1.values[other_ind]);
                    row = row + row1 * temp;
                    rhs = rhs + b1 * temp;
                    done = true;
                    break;
                }
                if (done)
                    break;
            }
        }
        // estimate and add best cut
        if (mir_cut.size() > 0) {
            // Scalar estimation = (violation * violation) / (lhv * lhv);
            Scalar estimation = (((mir_cut * violation) * objective_) * ((mir_cut * violation) * objective_)) / (mir_cut * mir_cut) / (objective_ * objective_); 
            cuts.push_back({mir_cut, gamma, estimation, false});
        }
    }
    return cuts;
}
