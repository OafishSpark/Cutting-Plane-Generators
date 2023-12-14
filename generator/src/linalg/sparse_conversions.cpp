#include "../../include/linalg/sparse_conversions.h"

SparseRowMatrix Col2Row(SparseColMatrix &scm) {
    SparseRowMatrix srm(scm.GetNRows(), scm.GetNCols());
    for (Index i_col = 0; i_col < Index(scm.GetNCols()); i_col++) {
        auto &col = scm.GetCol(i_col);
        for (Index i = 0; i < Index(col.size()); i++) {
            srm.GetRow(col.index[i]).PushBack(i_col, col.values[i]);
        }
    }
    return srm;
}

SparseColMatrix Row2Col(SparseRowMatrix &srm) {
    SparseColMatrix scm(srm.GetNRows(), srm.GetNCols());
    for (Index i_row = 0; i_row < Index(srm.GetNRows()); i_row++) {
        auto &row = srm.GetRow(i_row);
        for (Index i = 0; i < Index(row.size()); i++) {
            scm.GetCol(row.index[i]).PushBack(i_row, row.values[i]);
        }
    }
    return scm;
}
