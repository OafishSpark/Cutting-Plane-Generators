from copy import copy

import numpy as np
from numpy.linalg import inv

from pysmps import smps_loader as mps
from pyscipopt import LP, Model


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
    a_matrix = np.concatenate((a_matrix, np.eye(a_matrix.shape[0])), axis=1)
    return compute_b_inverse(lp, a_matrix) @ a_matrix


def write_matrix(matrix: np.array, name: str = "Instance") -> None:
    with open(name + ".txt", 'w') as f:
        answer = ""
        for line in matrix:
            answer += " ".join(list(map(str, line)))
            answer += '\n'
        f.write(answer)


def print_matrix(matrix: np.array) -> str:
    return '\n'.join([" ".join(list(map(str, line))) for line in matrix])


def get_prepared_for_gmi(filepath: str):
    # read .lp
    lp = LP()
    lp.readLP(bytes(filepath, encoding='utf-8'))
    milp = Model()
    milp.readProblem(filepath)
    # write .mps file
    new_filepath = copy(filepath).replace(".lp", ".mps")
    milp.writeProblem(new_filepath)
    # get data from .mps file
    parsed_mps = mps.load_mps(new_filepath)
    a_matrix = parsed_mps[7]
    is_integer = parsed_mps[4]
    rhs_signs = parsed_mps[5]
    rhs_values = [parsed_mps[9][elem] for elem in parsed_mps[8]][0]
    if (parsed_mps[10]):
        bounds_lo, bounds_up = parsed_mps[11][parsed_mps[10][0]].values()
    else:
        bounds_lo, bounds_up = [np.zeros((a_matrix.shape[1])), np.ones((a_matrix.shape[1])) * np.inf]
    lp.solve()
    solution = lp.getPrimal()
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
    return {
        "A": a_matrix, 
        "rhs_signs": rhs_signs,
        "rhs_values": rhs_values,
        "is_integer": is_integer,
        "bnd_lo": bounds_lo,
        "bnd_up": bounds_up,
        "basis_inds": basis_indexes,
        "sol": solution,
        "b_inv": b_inv
        }


def print_gmi_data(filepath: str):
    tmp_dict = get_prepared_for_gmi(filepath)
    answer = ""
    for key, value in zip(tmp_dict.keys(), tmp_dict.values()):
        answer += key + '\n'
        if key == "A" or key == "b_inv":
            answer += f'{value.shape[0]} {value.shape[1]}' + '\n'
            answer += print_matrix(value) + '\n'
        else:
            answer += " ".join(list(map(str, value))) + '\n'
    return answer
