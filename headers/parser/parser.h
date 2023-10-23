#ifndef PARSER_H
#define PARSER_H

#include "../utils/utils.h"

#include <vector>
#include <fstream>
#include <string>
#include <cassert>
#include <sstream>


std::vector<std::vector<Scalar>> ReadBinv(std::string filepath);

#endif //PARSER_H
