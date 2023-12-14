#include "cuts/cutter.h"


class gmi_tests {
public:
    struct gmi_general_test {
        // input
        SparseColMatrix a_matrix_;
        DenseMatrix b_inv_;
        RHS rhs_;
        Variables vars_;
        std::vector<Scalar> sol_;
        std::vector<int> basis_;
        // output
        std::vector<Cutter::Cut> cuts;
    };
    
    struct linear_transform_test {
        // input
        SparseColMatrix a_matrix_;
        RHS rhs_;
        Variables vars_;
        // output
        std::vector<bool> if_neg;
        std::vector<Scalar> shift_val;
    };

};
