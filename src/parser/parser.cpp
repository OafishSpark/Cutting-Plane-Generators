#include "../../headers/parser/parser.h"


std::vector<std::vector<Scalar>> ReadBinv(std::string filepath)
{
    DenseMatrix answer;
    std::ifstream b_inv_file(filepath);
    std::string line;
    while (std::getline(b_inv_file, line)) {
        answer.push_back(std::vector<Scalar>());
        std::istringstream line_stream(line);
        std::string elem;
        while (std::getline(line_stream, elem, ' ')) {
            if (elem != "") {
                Scalar d_elem = std::stod(elem);
                if (abs(d_elem) < Zero_Epsilon) {
                    d_elem = 0.0;
                }
                answer.back().push_back(d_elem);
            }
            
        }
    }
    b_inv_file.close();
    return answer;
}