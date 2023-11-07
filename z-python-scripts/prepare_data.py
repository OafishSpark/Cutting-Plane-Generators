from sys import argv
from os.path import dirname

from b_inv import print_matrix, b_inverse_a, b_inverse, get_prepared_for_gmi, print_gmi_data


try:
    _, filepath = argv
    print(filepath)
except ValueError:
    filepath = "files/task.lp"


if __name__ == "__main__":
    ans = print_gmi_data(filepath)
    new_path = "files/data.txt"
    with open(new_path, 'w') as f:
        f.write(ans)
