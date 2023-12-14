#ifndef PARSER_H
#define PARSER_H

#include "../utils/utils.h"
#include "../linalg/linalg.h"

#include <vector>
#include <fstream>
#include <string>
#include <cassert>
#include <iostream>
#include <sstream>


std::vector<std::vector<Scalar>> ReadBinv(std::string filepath);


class Model {
public:
    SparseColMatrix a_matrix_;
    DenseMatrix b_inv_;
    RHS rhs_;
    Variables vars_;
    std::vector<Scalar> sol_;
    std::vector<int> basis_;

    Model(std::string filepath);

    void Print();
};


#endif //PARSER_H
