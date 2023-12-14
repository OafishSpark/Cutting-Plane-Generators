#include "cuts/cutter.h"

Cutter::Cutter(MyUberModel& model, std::vector<Cutter::CutGenerator> cut_generators) {
    a_matrix_scm = model.a_matrix_scm_;
    a_matrix_srm_ = model.a_matrix_srm_;
    n_cols_ = model.n_cols_;
    n_rows_ = model.n_rows_;
    b_inv_ = model.b_inv_;
    rhs_ = model.rhs_;
    variables_ = model.variables_;
    rel_sol_ = model.rel_sol_;
    basis_ = model.basis_;
    cut_generators_ = cut_generators;
}

Cutter::Cutter(MyUberModel& model) {
    a_matrix_scm = model.a_matrix_scm_;
    a_matrix_srm_ = model.a_matrix_srm_;
    n_cols_ = model.n_cols_;
    n_rows_ = model.n_rows_;
    b_inv_ = model.b_inv_;
    rhs_ = model.rhs_;
    variables_ = model.variables_;
    rel_sol_ = model.rel_sol_;
    basis_ = model.basis_;
    cut_generators_ = std::vector<Cutter::CutGenerator>();
    cut_generators_.push_back(Cutter::CutGenerator({"GMI1", &Cutter::GomoryMixedIntegerCutGenerator}));
    cut_generators_.push_back(Cutter::CutGenerator({"MIR1", &Cutter::MixedIntegerRoundingCutGenerator}));
}

void Cutter::WriteCutsInFile(std::vector<Cut>& cuts) {
    std::string filepath = "files/temp-files/cuts.txt";
    std::ofstream file(filepath);
    for (auto& cut : cuts) {
        file << cut.lhs.ToStr() << " ";
        file << ">= ";
        file << cut.rhs << std::endl;
    }
    file.close();
}

bool CutCompare(const Cutter::Cut& fst, const Cutter::Cut& scnd) {
    return (fst.estimation < scnd.estimation);
}

size_t Cutter::RunGenerator(std::vector<Cut> (Cutter::*cut_generator_)()) {
    std::vector<Cut> cuts((this->*cut_generator_)());
    std::sort(cuts.begin(), cuts.end(), CutCompare);
    size_t t = cuts.size();
    if (t > kNCuts) {
        cuts.resize(kNCuts);
    }
    WriteCutsInFile(cuts);
    return cuts.size();
}

void Cutter::AddCuts() {
    for (auto generator : cut_generators_) {
        size_t n_c = RunGenerator(generator.cut_generator_);
        std::cout << generator.name_ << " generated: " << n_c << std::endl;
    }
}

void Cutter::AddCutsToHighsModel(Highs& model) {
    for (auto generator : cut_generators_) {
        std::vector<Cut> cuts((this->*generator.cut_generator_)());
        std::sort(cuts.begin(), cuts.end(), CutCompare);
        size_t t = cuts.size();
        if (t > kNCuts) cuts.resize(kNCuts);

        // add cuts to model
        for (int iv{0}; iv < cuts.size(); ++iv)
            model.addRow(cuts[iv].rhs, kHighsInf, cuts[iv].lhs.index.size(),
                         cuts[iv].lhs.index.data(), cuts[iv].lhs.values.data());

        std::cout << generator.name_ << " generated: " << cuts.size() << std::endl;
    }
}
