#include "./headers/cuts/cutter.h"
#include "./headers/parser/parser.h"
#include "./headers/utils/utils.h"
#include "./headers/linalg/linalg.h"

#include <iostream>
#include <string>
#include <cassert>
#include <vector>


int main(int argc, char* argv[]) {
    std::string filepath;
    if (argc > 1) {
        filepath = argv[1];
    } else {
        filepath = "files/data.txt";
    }
    std::cout << filepath << std::endl;

    std::vector<Scalar> da1 = {1, 2, 3, 4};
    SparseVector a1(da1);
    std::cout << a1 * a1 << std::endl;
    std::vector<Scalar> da2 = {0, 0, 0, 0};
    SparseVector a2(da2);
    std::cout << a1 * a2 << std::endl;
    std::vector<Scalar> da3 = {0, 2, 0, 4};
    SparseVector a3(da3);
    std::cout << a1 * a3 << std::endl;
    std::vector<Scalar> da4 = {1, 0, 0, 0};
    SparseVector a4(da4);
    std::cout << a3 * a4 << std::endl;
    std::cout << a4 * a4 << std::endl;
    

    Model model(filepath);

    Cutter c1(model);
    c1.AddCuts();

    std::cout << "Message!" << std::endl;
    return 0;
}
