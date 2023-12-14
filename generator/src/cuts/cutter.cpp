#include "../../include/cuts/cutter.h"


Cutter::Cutter(Model& model, std::vector<Cutter::CutGenerator> cut_generators) {
    a_matrix_ = model.a_matrix_;
    b_inv_ = model.b_inv_;
    rhs_ = model.rhs_;
    vars_ = model.vars_;
    sol_ = model.sol_;
    basis_ = model.basis_;
    cut_generators_ = cut_generators;
}

Cutter::Cutter(Model& model) {
    a_matrix_ = model.a_matrix_;
    b_inv_ = model.b_inv_;
    rhs_ = model.rhs_;
    vars_ = model.vars_;
    sol_ = model.sol_;
    basis_ = model.basis_;
    cut_generators_ = std::vector<Cutter::CutGenerator>({Cutter::CutGenerator({"GMI1", &Cutter::GomoryMixedIntegerCutGenerator})}); 
}

void Cutter::WriteCutsInFile(std::vector<Cut>& cuts) {
    std::string filepath = "files/cuts.txt";
    std::ofstream file(filepath);
    for (auto& cut: cuts) {
        file << cut.lhs.PrintInFile() << " ";
        file << ">= ";
        file << cut.rhs << std::endl;
    }
    file.close();
}

size_t Cutter::RunGenerator(std::vector<Cut> (Cutter::*cut_generator_)()) {
    std::vector<Cut> cuts = (this->*cut_generator_)();
    // for (auto& cut: cuts) {
    //     std::cout << cut.rhs; 
    //     cut.lhs.Print();
    //     std::cout << std::endl;
    // }
    WriteCutsInFile(cuts);
    return cuts.size();
}

void Cutter::AddCuts() {
    for (auto generator : cut_generators_) {
        size_t n_c = RunGenerator(generator.cut_generator_);
        std::cout << generator.name_ << " generated: " << n_c << std::endl;
    }
}
