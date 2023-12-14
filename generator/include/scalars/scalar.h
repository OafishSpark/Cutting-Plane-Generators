#ifndef ESOLVER_SCALAR_H
#define ESOLVER_SCALAR_H

#include <string>

#include "scalars/rational.h"

using Index = int32_t;

#define SCALAR_DOUBLE 1
#define SCALAR_FLOAT 2
#define SCALAR_RATIONAL 3

#ifndef SCALAR
#define SCALAR SCALAR_DOUBLE
// #define SCALAR SCALAR_FLOAT
// #define SCALAR SCALAR_RATIONAL
#endif

#if SCALAR == SCALAR_DOUBLE
using Scalar = double;
#elif SCALAR == SCALAR_FLOAT
using Scalar = float;
#elif SCALAR == SCALAR_RATIONAL
using Scalar = Rational<int32_t>;
#endif

Scalar str2scalar(std::string &s, std::size_t *pos);
std::string ScalarToString(const Scalar& scalar, const int width, const int precision);

#endif  // ESOLVER_SCALAR_H
