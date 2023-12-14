#ifndef ESOLVER_SPARSE_CONVERSIONS_H
#define ESOLVER_SPARSE_CONVERSIONS_H

#include "sparse_col_matrix.h"
#include "sparse_row_matrix.h"

SparseRowMatrix Col2Row(SparseColMatrix &scm);
SparseColMatrix Row2Col(SparseRowMatrix &srm);

#endif  // ESOLVER_SPARSE_CONVERSIONS_H
