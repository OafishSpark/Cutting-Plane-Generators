from sys import argv

from b_inv import print_matrix, b_inverse_a


_, filepath = argv

if __name__ == "__main__":
    print_matrix(b_inverse_a(filepath))
