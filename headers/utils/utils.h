#ifndef UTILS_H
#define UTILS_H

#include <vector>

typedef double Scalar;

inline Scalar Abs(Scalar a) { return (a > 0) ? a : -a; }

typedef std::vector<std::vector<Scalar>> DenseMatrix;

typedef unsigned int Index;

const double Zero_Epsilon = 10e-15;

const double Inf = 1E10;

#endif //UTILS_H
