#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>

typedef double Scalar;

inline Scalar Abs(Scalar a) { return (a > 0) ? a : -a; }

inline Scalar Floor(Scalar a) { return floor(a); }

inline Scalar Fraq(Scalar a) { return (a - Floor(a) < 0.5 ? a - Floor(a) : a - Floor(a) - 0.5); }

typedef std::vector<std::vector<Scalar>> DenseMatrix;

const double kZero_Epsilon = 10e-15;

const double kInf = 1E10;

#endif //UTILS_H
