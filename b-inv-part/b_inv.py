from copy import copy

import numpy as np
from numpy.linalg import inv

from pysmps import smps_loader as mps
from pyscipopt import LP


def read_lp_write_mps(filepath: str) -> (LP, np.array):
    # read .lp file
    lp = LP()
    lp.readLP(bytes(filepath, encoding='utf-8'))
    # write .mps file
    new_filepath = copy(filepath).replace(".lp", ".mps")
    lp.writeLP(bytes(new_filepath, encoding='utf-8'))
    # get instance matrix from .mps file
    a_matrix = mps.load_mps(new_filepath)[7]
    return lp, a_matrix


def compute_b_inverse(lp: LP, a_matrix: np.array) -> np.array:
    lp.solve()
    basis_indexes = lp.getBasisInds()
    # get columns of constraint matrix
    columns = np.transpose(a_matrix)
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
    b_inv = inv(np.transpose(b_matrix))
    return b_inv


def b_inverse(filepath: str) -> np.array:
    """Compute the inverse of basis matrix"""
    lp, a_matrix = read_lp_write_mps(filepath)
    return compute_b_inverse(lp, a_matrix)


def b_inverse_a(filepath: str) -> np.array:
    """Compute B^(-1) * A"""
    lp, a_matrix = read_lp_write_mps(filepath)
    return compute_b_inverse(lp, a_matrix) @ a_matrix


def write_matrix(matrix: np.array, name: str = "Instance") -> None:
    with open(name + ".txt", 'w') as f:
        answer = ""
        for line in matrix:
            answer += " ".join(list(map(str, line)))
            answer += '\n'
        f.write(answer)


def print_matrix(matrix: np.array) -> None:
    answer = ""
    for line in matrix:
        answer += " ".join(list(map(str, line)))
        answer += '\n'
    print(answer)
