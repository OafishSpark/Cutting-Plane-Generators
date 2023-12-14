#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>

#include "linalg/sparse_col_matrix.h"
#include "linalg/sparse_row_matrix.h"
#include "utils/highs_utils.h"

#include <Highs.h>


// enum class Target {
//     kMin,
//     kMax,
//     kUnknown,
// };


// struct Objective {
//     Objective();
//     Objective(Target target_, std::vector<Scalar> coeffs_);
//     Target target;
//     std::vector<Scalar> coeffs;
// };


// enum class VarType { kFixed, kLower, kUpper, kBoxed, kFree, kUnknown };


// struct Bounds {
//     Bounds(Scalar le_ = Scalar(0), Scalar ri_ = kInf): le(le_), ri(ri_) {}
//     Scalar le;  // >= le
//     Scalar ri;  // <= ri
// };


// struct Variables {
//     std::vector<Bounds> bounds;
//     std::vector<VarType> types;
//     std::vector<bool> integrality;
//     void Resize(size_t n) {
//         bounds.resize(n);
//         types.resize(n);
//         integrality.resize(n);
//     }
// };


// class MilpModel {
//     inline SparseRowMatrix &GetSrm() { return srm_; }
//     inline const SparseRowMatrix &GetSrm() const { return srm_; }
//     inline SparseColMatrix &GetScm() { return scm_; }
//     inline const SparseColMatrix &GetScm() const { return scm_; }

//     inline Objective &GetObjective() { return objective_; };
//     inline const Objective &GetObjective() const { return objective_; };
//     inline void SetObjective(const Objective &objective) { objective_ = objective; };

//     inline Variables CopyAllVariables() const { return var_; }
//     inline void SetAllVariables(const Variables &var) { var_ = var; }
//     inline bool GetIntegrality(Index iv) const { return var_.integrality[iv]; }
//     inline void SetIntegrality(Index iv, bool is_int) { var_.integrality[iv] = is_int; }
//     inline const Bounds &GetVarBounds(Index iv) const { return var_.bounds[iv]; }
//     void SetVarBounds(Index iv, const Bounds &bounds);
//     inline VarType GetVarType(Index iv) const { return var_.types[iv]; }

//     inline std::vector<Bounds> CopyAllRhs() { return rhs_; }
//     inline Bounds &GetRhs(Index ic) { return rhs_[ic]; }
//     inline const Bounds &GetRhs(Index ic) const { return rhs_[ic]; }
//     void SetRhs(Index ic, const Bounds &bounds) { rhs_[ic] = bounds; }

//     inline size_t GetNVars() const { return srm_.GetNCols(); }
//     inline size_t GetNCons() const { return srm_.GetNRows(); }

//     // !!! in cuts work over getters

//     SparseRowMatrix srm_;
//     SparseColMatrix scm_;

//     Objective objective_;
//     Variables var_;

//     std::vector<Bounds> rhs_;
// };


// class LpResult {
//     std::vector<Scalar> x;
//     Scalar y;
// };


// class DualRevisedSimplex {
//     inline const std::vector<Index> &GetBasis() const { return basis_; }
//     inline const SparseColMatrix &GetB() const { return B_; }
    
//     // !!! in cuts work over getters
    
//     SparseColMatrix B_;
//     std::vector<Index> basis_;
// }; 

class RHS {
public:
    struct Elem {
        char type_;
        Scalar val_;
    };

    std::vector<Elem> rhs_;

    RHS::Elem operator[](int ind) { return rhs_[ind]; }
};


class Variables {
public:
    struct Elem {
        bool is_int_;
        Scalar bnd_lo_;
        Scalar bnd_up_;
    };

    std::vector<Elem> vars_;

    Variables::Elem operator[](int ind) { return vars_[ind]; }
};



class MyUberModel {
public:
    // task data
    unsigned int n_cols_;
    unsigned int n_rows_;
    SparseColMatrix a_matrix_scm_;
    SparseRowMatrix a_matrix_srm_;
    SparseVector objective_;
    Variables variables_;
    RHS rhs_; 
    
    // relaxation data 
    std::vector<std::vector<Scalar>> b_inv_;
    std::vector<Scalar> rel_sol_;
    std::vector<Index> basis_;

    //esolver interface
    // MilpModel milp_;
    // DualRevisedSimplex drs_;
    // LpResult lp_res_;

    // ToDo: rewrite thing
    MyUberModel(std::string);
    MyUberModel() {}

    bool Compare(MyUberModel);
};

#endif //MODEL_H
