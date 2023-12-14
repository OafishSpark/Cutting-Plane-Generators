#ifndef TESTS_H
#define TESTS_H

#include "model/model.h"


class TestingSystem {
    MyUberModel model_test_1();
    MyUberModel model_test_2();
    MyUberModel model_test_3();

    struct ModelTest {
        std::string filepath;
        MyUberModel (TestingSystem::*test)();
    };
    std::vector<ModelTest> model_tests;

public:
    TestingSystem();
    bool RunModelTests();
};

#endif //TESTS_H
