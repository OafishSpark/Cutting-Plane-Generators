#ifndef HIGHS_UTILS_H
#define HIGHS_UTILS_H

#include <vector>
#include <stdio.h>

#include "Highs.h"


// for silencing printf in HiGHS
// see MyUberModel(std::string) constructor for usage
class CoolSilenceMode {
    bool mode;
    FILE* old_stdout;
public:  
    
    CoolSilenceMode(bool flag = true) {
        mode = true;
        old_stdout = stdout;
        if (flag) {
            stdout = fopen("/dev/null", "w");
        }
    }

    ~CoolSilenceMode() {
        mode = false;
        stdout = old_stdout;
    }

    void Off() {
        mode = false;
    }

    inline bool Mode() {
        return mode;
    } 
};


std::vector<std::vector<double>> DensifyColwiseHighsMatrix(const HighsSparseMatrix&);


HighsModel CreateLpRelaxation(const HighsModel&);


#endif //HIGHS_UTILS_H
