#include "testing/tests.h"
#include <lp_data/HConst.h>
#include "linalg/sparse_col_matrix.h"
#include "linalg/sparse_vector.h"


MyUberModel TestingSystem::model_test_1() {
    MyUberModel m;
    
    m.n_cols_ = 2;
    m.n_rows_ = 2;
    
    std::vector<std::vector<Scalar>> temp_mat = {
        {1, 1},
        {-1, 1}
    };
    
    m.a_matrix_scm_ = SparseColMatrix(temp_mat);

    m.objective_ = SparseVector(std::vector<Scalar>({1, -1}));

    m.variables_.vars_.push_back({true, 0, 2});
    m.variables_.vars_.push_back({false, 0, kHighsInf});

    m.rhs_.rhs_.push_back({'l', 6});
    m.rhs_.rhs_.push_back({'u', 1});
    
    m.b_inv_.push_back({1, 0});
    m.b_inv_.push_back({-1, 1});

    m.basis_ = std::vector<Index>({1, -2});

    m.rel_sol_ = std::vector<Scalar>({0, 6});

    return m;
}


MyUberModel TestingSystem::model_test_2() {
    MyUberModel m;
    
    m.n_cols_ = 2;
    m.n_rows_ = 3;
    
    std::vector<std::vector<Scalar>> temp_mat = {
        {1, 1},
        {1, -1},
        {1, -3}
    };
    
    m.a_matrix_scm_ = SparseColMatrix(temp_mat);

    m.objective_ = SparseVector(std::vector<Scalar>({-0.1, -10}));

    m.variables_.vars_.push_back({false, 3.3, 4.4});
    m.variables_.vars_.push_back({true, 0, kHighsInf});

    m.rhs_.rhs_.push_back({'l', 6});
    m.rhs_.rhs_.push_back({'u', 0.5});
    m.rhs_.rhs_.push_back({'e', -2.5});
    
    m.b_inv_.push_back({3.0/4, 0, 1.0/4});
    m.b_inv_.push_back({1.0/4, 0, -1.0/4});
    m.b_inv_.push_back({-1.0/2, 1,  -1.0/2});

    m.basis_ = std::vector<Index>({0, 1, -2});

    m.rel_sol_ = std::vector<Scalar>({31.0/8, 17.0/8});

    return m;
}


MyUberModel TestingSystem::model_test_3() {
    MyUberModel m;
    
    m.n_cols_ = 3;
    m.n_rows_ = 2;
    
    std::vector<std::vector<Scalar>> temp_mat = {
        {-11, 13, 0},
        {0, 10, 0},
    };
    
    m.a_matrix_scm_ = SparseColMatrix(temp_mat);

    m.objective_ = SparseVector(std::vector<Scalar>({11, -13.1}));

    m.variables_.vars_.push_back({true, -10, 10});
    m.variables_.vars_.push_back({true, 0, kHighsInf});
    m.variables_.vars_.push_back({false, -1, 1});

    m.rhs_.rhs_.push_back({'l', -10});
    m.rhs_.rhs_.push_back({'l', 90});
    
    m.b_inv_.push_back({1.0/13, 0});
    m.b_inv_.push_back({-10.0/13, 1});

    m.basis_ = std::vector<Index>({ 1, -2});

    m.rel_sol_ = std::vector<Scalar>({10, 100.0/13, -1});

    return m;
}


TestingSystem::TestingSystem() {
    // general test #1
    ModelTest test1({"files/tests/matrix_parse_test_1.lp", &TestingSystem::model_test_1});
    model_tests.push_back(test1);
    // general test #2
    ModelTest test2({"files/tests/matrix_parse_test_2.lp", &TestingSystem::model_test_2});
    model_tests.push_back(test2);
    // general test #3
    ModelTest test3({"files/tests/matrix_parse_test_3.lp", &TestingSystem::model_test_3});
    model_tests.push_back(test3);
}


bool TestingSystem::RunModelTests() {
    for (const auto& test: model_tests) {
        MyUberModel test_model = (this->*test.test)();
        MyUberModel alg_model = MyUberModel(test.filepath);
        bool result = test_model.Compare(alg_model);
        if (!result) {
            return result;
        }
    }
    return true;
}
