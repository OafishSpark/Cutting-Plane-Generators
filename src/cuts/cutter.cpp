#include "../../headers/cuts/cutter.h"


Cutter::Cutter(std::vector<Cutter::CutGenerator> cut_generators) {
    cut_generators_ = cut_generators;
}

Cutter::Cutter() {
    std::cout << "Nothing happens\n";
}
