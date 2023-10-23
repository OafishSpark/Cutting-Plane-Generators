from sys import argv

from b_inv import print_matrix, b_inverse_a, b_inverse, get_prepared_for_gmi, print_gmi_data

try:
    _, filepath = argv
except ValueError:
    filepath = "./b-inv-part/tests/task.lp"

if __name__ == "__main__":
    # print_matrix(b_inverse_a(filepath))
    # print_matrix(b_inverse(filepath))
    # print(get_prepared_for_gmi(filepath))
    print_gmi_data(filepath)
