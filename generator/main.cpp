#include "include/cuts/cutter.h"
#include "include/parser/parser.h"
#include "include/utils/utils.h"
#include "include/linalg/linalg.h"

#include <iostream>
#include <string>
#include <cassert>
#include <vector>


int main(int argc, char* argv[]) {
    std::string filepath;
    if (argc > 1) {
        filepath = argv[1];
    } else {
        filepath = "../files/data.txt";
    }
    std::cout << filepath << std::endl;

    Model model(filepath);

    Cutter c1(model);
    c1.AddCuts();
    return 0;
}
