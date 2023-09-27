#ifndef CUTTER_H
#define CUTTER_H

#include "../utils/utils.h"

#include <vector>
#include <string>

class Cutter {
    struct Cut{
        // SparseVector lhs;
        Scalar rhs;
        Scalar estimation;
        bool removed;
    };

    struct CutGenerator{
        std::string name;
        std::vector<Cut> (Cutter::*cut_generator)();
    };

    std::vector<CutGenerator> cut_generators_;
};

#endif CUTTER_H
