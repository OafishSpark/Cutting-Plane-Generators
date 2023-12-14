#include <lp_data/HStruct.h>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <stdio.h>

#include <Highs.h>

#include "model/model.h"
#include "utils/highs_utils.h"
#include "testing/tests.h"

#include "cuts/cutter.h"


void RunCutsGenerator(std::string input_filepath, std::string output_filepath) {
    // solver for original model
    MyUberModel milp(input_filepath);
    
    // cutting planes generator
    Cutter cut_generator(milp);
    
    for (CoolSilenceMode silence(false); silence.Mode(); silence.Off()) {
        Highs cutted_model;
        cutted_model.readModel(input_filepath);
        cut_generator.AddCutsToHighsModel(cutted_model);
        cutted_model.writeModel(output_filepath);
        
        return;

        cutted_model.run();
        std::vector<double> cutted_sol = cutted_model.getSolution().col_value;

        Highs model;

        model.setHighsOptionValue("mip_pool_soft_limit","0");

        model.readModel(input_filepath);
        model.run();
        std::vector<double> sol = model.getSolution().col_value;

        assert(cutted_sol.size() == sol.size());
        bool correctness = true;
        for (int iv{0}; iv < sol.size(); ++iv) {
            if (Abs(cutted_sol[iv] - sol[iv]) >= kZeroEpsilon) {
                std::cout << "\nCuts cut off the relaxation solution!!!\n";
                std::cout << "Error in index: " << iv << "\n";
                std::cout << "Cutted value is: " << cutted_sol[iv] << "\n";
                std::cout << "Original value is: " << sol[iv] << "\n";
                correctness = false;
                break;
            }
        }
        if (correctness) {
            std::cout << "Cuts do not cut off the relaxation solution!\n";
        }
    }
}


int main() { 
    // Todo: these strings as console parameters
    std::string input_filepath = "files/instances/noswot.mps";
    std::string output_filepath = "files/cutted-instances/noswot_cutted.lp";
    
    TestingSystem tests;
    if (tests.RunModelTests()) {
        std::cout << "All Tests Are Passed!\n";
    }
    
    RunCutsGenerator(input_filepath, output_filepath);

    return 0;
}
