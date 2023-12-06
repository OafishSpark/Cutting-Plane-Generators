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

    Model model(filepath);

    Cutter c1(model);
    c1.AddCuts();
    return 0;
}
