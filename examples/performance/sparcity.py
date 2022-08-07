import timeit

import casadi as ca

matrix = ca.SX([[0, 1, 2, 3], [4, 5, 6, 7]])
sparsity = ca.Sparsity(2, 4, [0, 0, 0, 1, 1], [0])


def performance_get_true():
    matrix.get(True, sparsity)


def performance_get_false():
    matrix.get(False, sparsity)


def performance_project():
    ca.project(matrix, sparsity)


def performance_index():
    matrix[sparsity]


for function in [
    performance_get_true,
    performance_get_false,
    performance_project,
    performance_index,
]:
    print(
        f"{timeit.timeit(function, globals=globals(), number=100000)} : {function.__name__}"
    )
