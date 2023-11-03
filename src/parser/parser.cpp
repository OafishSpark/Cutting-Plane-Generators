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
                if (Abs(d_elem) < kZero_Epsilon) {
                    d_elem = 0.0;
                }
                answer.back().push_back(d_elem);
            }
            
        }
    }
    b_inv_file.close();
    return answer;
}

Model::Model(std::string filepath) {
    std::ifstream file(filepath);
    std::string line;
    std::istringstream line_stream;
    std::string elem;
    while (std::getline(file, line)) {
        if (line == "") {
            break;
        } else if (line == "A") {
            std::getline(file, line);
            line_stream = std::istringstream(line);
            std::getline(line_stream, elem, ' ');
            int rows = std::stoi(elem);
            std::getline(line_stream, elem, ' ');
            int cols = std::stoi(elem);
            DenseMatrix a_temp;
            a_temp.resize(rows, std::vector<Scalar>(cols));
            for (int iv = 0; iv < rows; ++iv) {
                std::getline(file, line);
                line_stream = std::istringstream(line);
                for (int jv = 0; jv < cols; ++jv) {
                    std::getline(line_stream, elem, ' ');
                    Scalar temp = std::stod(elem);
                    if (Abs(temp) < kZero_Epsilon) {
                        temp = 0.0;
                    }
                    a_temp[iv][jv] = temp;
                }
            }
            a_matrix_ = SparseColMatrix(a_temp);
        } else if (line == "rhs_signs") {
            std::getline(file, line);
            line_stream = std::istringstream(line);
            rhs_ = RHS();
            rhs_.rhs_.resize(a_matrix_.rows_, RHS::Elem({'L', 0}));
            for (int iv = 0; iv < a_matrix_.rows_; ++iv) {
                std::getline(line_stream, elem, ' ');
                rhs_.rhs_[iv].type_ = elem[0];
            }
            std::getline(file, line);
            assert(line == "rhs_values");
            std::getline(file, line);
            line_stream = std::istringstream(line);
            for (int iv = 0; iv < a_matrix_.rows_; ++iv) {
                std::getline(line_stream, elem, ' ');
                Scalar temp = std::stod(elem);
                if (Abs(temp) < kZero_Epsilon) {
                    temp = 0.0;
                }
                rhs_.rhs_[iv].val_ = temp;
            }
        } else if (line == "is_integer") {
            std::getline(file, line);
            line_stream = std::istringstream(line);
            vars_ = Variables();
            vars_.vars_.resize(a_matrix_.cols_, Variables::Elem({false, 0.0, kInf}));
            for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
                std::getline(line_stream, elem, ' ');
                if (elem == "integral") {
                    vars_.vars_[iv].is_int_ = true;
                } else if (elem == "continuous") {
                    vars_.vars_[iv].is_int_ = false;
                } else {
                    assert(false);
                }
            }
            std::getline(file, line);
            assert(line == "bnd_lo");
            std::getline(file, line);
            line_stream = std::istringstream(line);
            for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
                std::getline(line_stream, elem, ' ');
                if (elem == "-inf") {
                    elem = "-1E10";
                }
                Scalar temp = std::stod(elem);
                if (Abs(temp) < kZero_Epsilon) {
                    temp = 0.0;
                }
                vars_.vars_[iv].bnd_lo_ = temp;
            }
            std::getline(file, line);
            assert(line == "bnd_up");
            std::getline(file, line);
            line_stream = std::istringstream(line);
            for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
                std::getline(line_stream, elem, ' ');
                if (elem == "inf") {
                    elem = "1E10";
                }
                Scalar temp = std::stod(elem);
                if (Abs(temp) < kZero_Epsilon) {
                    temp = 0.0;
                }
                vars_.vars_[iv].bnd_up_ = temp;
            }
        } else if (line == "basis_inds") {
            basis_ = std::vector<int>();
            basis_.resize(a_matrix_.rows_, 0);
            std::getline(file, line);
            line_stream = std::istringstream(line);
            for (int iv = 0; iv < a_matrix_.rows_; ++iv) {
                std::getline(line_stream, elem, ' ');
                basis_[iv] = std::stoi(elem);
            }
        } else if (line == "sol") {
            sol_ = std::vector<Scalar>();
            sol_.resize(a_matrix_.cols_, 0.0);
            std::getline(file, line);
            line_stream = std::istringstream(line);
            for (int iv = 0; iv < a_matrix_.cols_; ++iv) {
                std::getline(line_stream, elem, ' ');
                Scalar temp = std::stod(elem);
                if (Abs(temp) < kZero_Epsilon) {
                    temp = 0.0;
                }
                sol_[iv] = temp;
            }
        } else if (line == "b_inv") {
            std::getline(file, line);
            line_stream = std::istringstream(line);
            std::getline(line_stream, elem);
            int rows = std::stoi(elem);
            std::getline(line_stream, elem);
            int cols = std::stoi(elem);
            b_inv_ = DenseMatrix();
            b_inv_.resize(rows, std::vector<Scalar>(cols));
            for (int iv = 0; iv < rows; ++iv) {
                std::getline(file, line);
                line_stream = std::istringstream(line);
                for (int jv = 0; jv < cols; ++jv) {
                    std::getline(line_stream, elem, ' ');
                    Scalar temp = std::stod(elem);
                    if (Abs(temp) < kZero_Epsilon) {
                        temp = 0.0;
                    } 
                    b_inv_[iv][jv] = temp; 
                }
            }
        } else {
            if (line == "") {
                continue;
            }
            assert(false);
        }
    }
}

void Model::Print() {

}
