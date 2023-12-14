#ifndef ESOLVER_SCALAR_H
#define ESOLVER_SCALAR_H

#include <string>

#include "scalars/rational.h"

using Index = int32_t;

#define SCALAR_DOUBLE 1

#ifndef SCALAR
#define SCALAR SCALAR_DOUBLE
#endif

#if SCALAR == SCALAR_DOUBLE
using Scalar = double;
#endif

Scalar str2scalar(std::string &s, std::size_t *pos);
std::string ScalarToString(const Scalar& scalar, const int width, const int precision);

#endif  // ESOLVER_SCALAR_H
