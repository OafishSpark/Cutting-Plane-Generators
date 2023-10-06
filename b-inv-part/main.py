from copy import copy

import numpy as np
from numpy.linalg import inv

from pysmps import smps_loader as mps
from pyscipopt import LP


def b_inverse(filepath: str):
    """ b_inverse(filepath) -> b_inv
        --- computes inverse basis matrix
    Input:
        filepath --- path to .lp file with task description
    Output:
        b_inv --- inverse basis matrix
    """
    # read .lp file
    lp = LP()
    lp.readLP(bytes(filepath, encoding='utf-8'))
    # solve LP relaxation and get basis indexes
    lp.solve()
    basis_indexes = lp.getBasisInds()
    # write and read with pysmps .mps file
    new_filepath = copy(filepath).replace(".lp", ".mps")
    lp.writeLP(bytes(new_filepath, encoding='utf-8'))
    columns = np.transpose(mps.load_mps(new_filepath)[7])
    # compute basis matrix
    k = len(basis_indexes)
    b_matrix = []
    for index in basis_indexes:
        if index >= 0:
            b_matrix.append(columns[index])
        else:
            column = [0] * k
            column[-index - 1] = 1
            b_matrix.append(column)
    b_matrix = np.array(b_matrix)
    print(b_matrix)
    b_inv = inv(np.transpose(b_matrix))
    return b_inv


if __name__ == "__main__":
    print(b_inverse("tests/task.lp"))
