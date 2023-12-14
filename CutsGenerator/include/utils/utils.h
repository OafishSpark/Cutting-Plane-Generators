#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>
#include <limits>

#include "scalars/scalar.h"

inline Scalar Abs(Scalar a) { return (a > 0) ? a : -a; }

inline Scalar Floor(Scalar a) { return std::floor(a); }

inline Scalar Fraq(Scalar a) { return (Abs(a) - Floor(Abs(a)) < 0.5 ? Abs(a) - Floor(Abs(a)) : Abs(a) - Floor(Abs(a)) - 0.5); }

typedef std::vector<std::vector<Scalar>> DenseMatrix;

constexpr Scalar kInf = std::numeric_limits<double>::infinity();
constexpr Scalar kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr Scalar k05 = Scalar(0.5);
constexpr Scalar kMaxAbs = Scalar(1.e10);
constexpr Scalar kEpsInteger = Scalar(1.e-5);
constexpr Scalar kZeroEpsilon = Scalar(1.e-10);

#endif //UTILS_H
