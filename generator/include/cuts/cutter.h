#ifndef CUTTER_H
#define CUTTER_H

#include "../utils/utils.h"
#include "../parser/parser.h"
#include "../linalg/linalg.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>


class Cutter {
    const Scalar kAway = 0.01;
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

    SparseColMatrix a_matrix_;
    DenseMatrix b_inv_;
    RHS rhs_;
    Variables vars_;
    std::vector<Scalar> sol_;
    std::vector<int> basis_;

    Cutter(Model&, std::vector<CutGenerator> cut_generators);
    
    Cutter(Model&);

    void AddCuts();

    void WriteCutsInFile(std::vector<Cut>&);

    size_t RunGenerator(std::vector<Cut> (Cutter::*)());

    std::vector<Cut> GomoryMixedIntegerCutGenerator();
};

#endif //CUTTER_H
