#ifndef ESOLVER_LINALG_RAW_PARAMS_H
#define ESOLVER_LINALG_RAW_PARAMS_H

#include "scalars/scalar.h"

constexpr Index split_step = Index(4);
constexpr Index split_size = Index(1) << split_step;
constexpr Index split_size_half = split_size >> 1;
constexpr Index split_size_onehalf = split_size + split_size_half;

constexpr Index matrix_split(Index n) {
    return (n >= split_size) ? ((n + split_size_half) / split_size) * split_size_half : n / 2;
};

#endif // ESOLVER_LINALG_RAW_PARAMS_H

