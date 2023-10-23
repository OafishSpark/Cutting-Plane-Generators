#include "./headers/cuts/cutter.h"
#include "./headers/parser/parser.h"
#include "./headers/utils/utils.h"

#include <iostream>
#include <string>
#include <cassert>


int main(int argc, char* argv[]) {
    std::string filepath;
    if (argc > 1) {
        filepath = argv[1];
    } else {
        filepath = "files/b-inv.txt";
    }
    std::cout << filepath << std::endl;
    DenseMatrix b_inv = ReadBinv(filepath);
    for (const auto& dense_vector: b_inv) {
        for (const auto& elem: dense_vector) {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
    Cutter c1;
    std::cout << "Message!" << std::endl;
    return 0;
}
