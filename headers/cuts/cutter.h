#ifndef CUTTER_H
#define CUTTER_H

#include "../utils/utils.h"
#include "../parser/parser.h"
#include "../linalg/linalg.h"

#include <vector>
#include <string>
#include <iostream>


class Cutter {
public:
    struct Cut{
        SparseVector lhs;
        Scalar rhs;
        Scalar estimation;
        bool removed;
    };

    struct CutGenerator{
        std::string name;
        std::vector<Cut> (Cutter::*cut_generator)();
    };

    std::vector<CutGenerator> cut_generators_;

    Model* model_;

    Cutter(Model*, std::vector<CutGenerator> cut_generators);
    
    Cutter(Model*);

    std::vector<Cut> GomoryMixedIntegerCutGenerator();
};

#endif //CUTTER_H
