from sys import argv
from os.path import dirname

from b_inv import print_matrix, b_inverse_a, b_inverse, get_prepared_for_gmi, print_gmi_data


try:
    _, filepath = argv
    print(filepath)
except ValueError:
    # filepath = "./b-inv-part/tests/task.lp"
    filepath = "./files/task.lp"


if __name__ == "__main__":
    # print_matrix(b_inverse_a(filepath))
    # print_matrix(b_inverse(filepath))
    # print(get_prepared_for_gmi(filepath))
    ans = print_gmi_data(filepath)
    new_path = dirname(filepath) + "/data.txt"
    with open(new_path, 'w') as f:
        f.write(ans)
