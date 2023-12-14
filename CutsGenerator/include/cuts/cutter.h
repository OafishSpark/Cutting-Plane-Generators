#ifndef CUTTER_H
#define CUTTER_H

#include "linalg/sparse_row_matrix.h"
#include "utils/utils.h"
#include "model/model.h"
#include "linalg/sparse_vector.h"
#include "linalg/sparse_col_matrix.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <tuple>


class Cutter {
    const size_t kNCuts = 500;
    const Scalar kAway = Scalar(1) / Scalar(100);
    const Scalar kMaxDensity = Scalar(5) / Scalar(100);
    const Scalar kMaxSlack = Scalar(10) / Scalar(100);

    // additional methods for GMI
    std::pair<std::vector<bool>, std::vector<Scalar>> LinearTransform();

public:
    struct Cut{
        SparseVector lhs;
        Scalar rhs;
        Scalar estimation;
        bool removed;
    };

    struct CutGenerator{
        std::string name_;
        std::vector<Cut> (Cutter::*cut_generator_)();
    };

    std::vector<CutGenerator> cut_generators_;

    unsigned int n_cols_;
    unsigned int n_rows_;
    SparseColMatrix a_matrix_scm;
    SparseRowMatrix a_matrix_srm_;
    SparseVector objective_;
    Variables variables_;
    RHS rhs_; 
    std::vector<std::vector<Scalar>> b_inv_;
    std::vector<Scalar> rel_sol_;
    std::vector<Index> basis_;

    Cutter(MyUberModel&, std::vector<CutGenerator> cut_generators);
    
    Cutter(MyUberModel&);

    // todo
    // Cutter(MilpModel&, DualRevisedSimplex&, LpResult&);

    void AddCuts();

    void WriteCutsInFile(std::vector<Cut>&);

    size_t RunGenerator(std::vector<Cut> (Cutter::*)());

    std::vector<Cut> GomoryMixedIntegerCutGenerator();

    std::vector<Cut> MixedIntegerRoundingCutGenerator();

    // Adds cuts to highs model
    void AddCutsToHighsModel(Highs& model);
};

#endif //CUTTER_H
