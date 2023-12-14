#include "utils/highs_utils.h"
#include <lp_data/HConst.h>


std::vector<std::vector<double>> DensifyColwiseHighsMatrix(const HighsSparseMatrix& highs_matrix) {
    auto& num_col = highs_matrix.num_col_;
    auto& num_row = highs_matrix.num_row_;
    auto& start = highs_matrix.start_;
    auto& index = highs_matrix.index_;
    auto& value = highs_matrix.value_;
    
    std::vector<std::vector<double>> matrix(num_row, std::vector<double>(num_col, 0.0));
    
    for (int iv{0}; iv < start.size() - 1; ++iv) {
        for (int jv{start[iv]}; jv < start[iv + 1]; ++jv) {
            matrix[index[jv]][iv] = value[jv];
        }
    }
    
    return matrix;
}


HighsModel CreateLpRelaxation(const HighsModel& model) {
    HighsModel lp_relaxation(model);
    for (auto& integrality: lp_relaxation.lp_.integrality_) {
        integrality = HighsVarType::kContinuous;
    }
    return lp_relaxation;
}
