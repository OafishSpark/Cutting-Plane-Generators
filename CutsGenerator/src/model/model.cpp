#include "model/model.h"
#include <lp_data/HConst.h>
#include "linalg/sparse_col_matrix.h"
#include "linalg/sparse_row_matrix.h"
#include "linalg/sparse_conversions.h"
#include "linalg/sparse_vector.h"
#include "utils/utils.h"



MyUberModel::MyUberModel(std::string filepath) {
    for (CoolSilenceMode silence; silence.Mode(); silence.Off()) {
        Highs solver;
        solver.readModel(filepath);

        const HighsLp& lp = solver.getLp();

        std::vector<std::vector<Scalar>> temp_a = DensifyColwiseHighsMatrix(lp.a_matrix_);
        a_matrix_scm_ = SparseColMatrix(temp_a);
        a_matrix_srm_ = Col2Row(a_matrix_scm_);
        
        assert(a_matrix_scm_.GetNRows() > 0);
        n_rows_ = a_matrix_scm_.GetNRows();
        n_cols_ = a_matrix_scm_.GetNCols();

        objective_ = SparseVector(lp.col_cost_);

        // variables
        variables_.vars_.resize(n_cols_);
        for (int iv{0}; iv < n_cols_; ++iv) {
            auto& var = variables_.vars_[iv];
            var.bnd_lo_ = Scalar(lp.col_lower_[iv]);
            var.bnd_up_ = Scalar(lp.col_upper_[iv]);
            var.is_int_ = (lp.integrality_[iv] == HighsVarType::kInteger);
        }

        // rhs
        rhs_.rhs_.resize(n_rows_);
        for (int iv{0}; iv < n_rows_; ++iv) {
            auto& rhs = rhs_.rhs_[iv];
            auto& lo = lp.row_lower_[iv];
            auto& up = lp.row_upper_[iv];
            if (lo == -kHighsInf) {
                // upper or equal inequality
                assert(up < kHighsInf);
                rhs.type_ = 'l';
                rhs.val_ = up;
            } else if (up == kHighsInf) {
                // lower or equal inequality
                assert(lo > -kHighsInf);
                rhs.type_ = 'u';
                rhs.val_ = lo;
            } else if (up == lo) {
                // equality
                assert(up < kHighsInf);
                rhs.type_ = 'e';
                rhs.val_ = up;
            } else {
                // rhs is boxed
                assert(false);
            }
        }
        
        // lp relaxation
        const HighsModel& model = solver.getModel();
        HighsModel lp_relaxation = CreateLpRelaxation(model);
        
        Highs relaxation_solver;
        relaxation_solver.setHighsOptionValue("presolve", "off");

        relaxation_solver.passModel(lp_relaxation);
        relaxation_solver.run();

        // relaxation basis
        basis_ = std::vector<Index>(n_rows_);
        relaxation_solver.getBasicVariables(basis_.data());

        // basis inverse matrix
        b_inv_ = std::vector<std::vector<Scalar>>(n_rows_, std::vector<Scalar>(n_rows_, 0));
        for (int iv{0}; iv < basis_.size(); ++iv) {
            relaxation_solver.getBasisInverseRow(iv, b_inv_[iv].data());
        }

        // relaxation solution
        const HighsSolution& rel_sol = relaxation_solver.getSolution();
        rel_sol_ = std::vector<Scalar>(rel_sol.col_value);
    }
}


bool MyUberModel::Compare(MyUberModel other) {
    // 1. n_rows
    if (n_rows_ != other.n_rows_) {
        std::cout << "#1 Number of rows is not equal" << std::endl;
        return false;
    }

    // 2. n_cols
    if (n_cols_ != other.n_cols_) {
        std::cout << "#2 Number of cols is not equal" << std::endl;
        return false;
    }

    // 3. a_matrix
    for (int iv{0}; iv < n_cols_; ++iv) {
        const auto& orig_col = a_matrix_scm_.GetCols()[iv];
        const auto& other_col = other.a_matrix_scm_.GetCols()[iv];
        if (orig_col.size() != other_col.size()) {
            std::cout << "#3 Column sizes in matrices A are different" << std::endl;
            return false;
        }
        for (int jv{0}; jv < orig_col.size(); ++jv) {
            if (orig_col.index[jv] != other_col.index[jv]) {
                std::cout << "#3 Indexes in " << iv << " sparse columns of matrices A are different" << std::endl;
                return false;
            }
            if (abs(orig_col.values[jv] - other_col.values[jv]) >= kZeroEpsilon) {
                std::cout << "#3 Values in " << iv << " sparse columns of matrices A are different" << std::endl;
                return false;
            }
        }
    }

    // 4. objective
    if (objective_.size() != other.objective_.size()) {
        std::cout << "#4 Column sizes in objectives are different" << std::endl;
        return false;
    }
    for (int iv{0}; iv < objective_.size(); ++iv) {
        if (objective_.index[iv] != other.objective_.index[iv]) {
            std::cout << "#4 Indexes in sparse columns of objectives are different" << std::endl;
            return false;
        }
        if (abs(objective_.values[iv] - other.objective_.values[iv]) >= kZeroEpsilon) {
            std::cout << "#4 Values in sparse columns of objectives are different" << std::endl;
            return false;
        }
    }

    // 5. Variables
    if (variables_.vars_.size() != other.variables_.vars_.size()) {
        std::cout << "#5 Variables numbers are different" << std::endl;
        return false;
    }
    for (int iv{0}; iv < variables_.vars_.size(); ++iv) {
        if (variables_[iv].is_int_ != other.variables_[iv].is_int_) {
            std::cout << "#5 Integrality of variables is different" << std::endl;
            return false;
        }
        if (abs(variables_[iv].bnd_lo_ - other.variables_[iv].bnd_lo_) >= kZeroEpsilon) {
            std::cout << "#5 Lower bounds of variables are different" << std::endl;
            return false;
        }
        if (abs(variables_[iv].bnd_up_ - other.variables_[iv].bnd_up_) >= kZeroEpsilon) {
            std::cout << "#5 Upper bounds of variables are different" << std::endl;
            return false;
        }
    }

    // 6. RHS
    if (rhs_.rhs_.size() != other.rhs_.rhs_.size()) {
        std::cout << "#6 RHS numbers are different" << std::endl;
        return false;
    }
    for (int iv{0}; iv < rhs_.rhs_.size(); ++iv) {
        if (rhs_[iv].type_ != other.rhs_[iv].type_) {
            std::cout << "#6 RHS types are different" << std::endl;
            return false;
        }
        if (abs(rhs_[iv].val_ - other.rhs_[iv].val_) >= kZeroEpsilon) {
            std::cout << "#6 RHS values are different" << std::endl;
            return false;
        }
    }

    // 8. Relaxation solution
    assert(rel_sol_.size() == n_cols_);
    assert(other.rel_sol_.size() == n_cols_);
    for (int iv{0}; iv < n_cols_; ++iv) {
        if (abs(rel_sol_[iv] - other.rel_sol_[iv]) >= kZeroEpsilon) {
            std::cout << "#8 Relaxation solution values are different"<< std::endl;
            return false;
        }
    }

    // 9. Basis
    if (basis_.size() != other.basis_.size()) {
        std::cout << "#9 Basis dimensions are different" << std::endl;
        return false;
    }
      for (int iv{0}; iv < n_rows_; ++iv) {
        if (std::find(other.basis_.begin(), other.basis_.end(), basis_[iv]) == other.basis_.end()) {
            std::cout << "#9 Basis values are diferent"<< std::endl;
            return false;
        }
    }

    // 7. B inverse (lol)
    assert(b_inv_.size() == n_rows_);
    assert(other.b_inv_.size() == n_rows_);
    for (int iv{0}; iv < n_rows_; ++iv) {
        for (int jv{0}; jv < n_rows_; ++jv) {
            Scalar other_result = Scalar(other.basis_[jv] == other.basis_[iv]);
            SparseVector col;
            if (other.basis_[jv] >= 0) {
                col = other.a_matrix_scm_.GetCol(other.basis_[jv]);
            } else {
                std::vector<Scalar> temp(n_rows_, 0);
                temp[-other.basis_[jv] - 1] = 1;
                col = SparseVector(temp);
            }
            Scalar other_res = SparseVector(other.b_inv_[iv]) * col;
            assert(abs(other_res - other_result) < kZeroEpsilon);
            Scalar orig_result = Scalar(basis_[jv] == basis_[iv]);
            if (basis_[jv] >= 0) {
                col = a_matrix_scm_.GetCol(basis_[jv]);
            } else {
                std::vector<Scalar> temp(n_rows_, 0);
                temp[-basis_[jv]-1] = 1;
                col = SparseVector(temp);
            }
            Scalar orig_res = SparseVector(b_inv_[iv]) * col;
            assert(abs(orig_res - orig_result) < kZeroEpsilon);
        }
        
    }

    // All fields are the same
    return true;
}
