#include "../../headers/cuts/cutter.h"


Cutter::Cutter(Model* model, std::vector<Cutter::CutGenerator> cut_generators) {
    model_ = model;
    cut_generators_ = cut_generators;
}

Cutter::Cutter(Model* model) {
    model_ = model;
    std::cout << "Nothing happens\n";
}
