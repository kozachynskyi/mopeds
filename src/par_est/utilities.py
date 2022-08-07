""" Here comes methods that are used by main library"""
import copyreg
import pickle
import tempfile
import webbrowser

import casadi as ca
import numpy as np
import pandas as pd


def reduce_MX(mx_object):
    return (str, (repr(mx_object),))


class MXPickler(pickle.Pickler):
    """This class allows to pickle casadi objects by replacing
    them with their __repr__ values. Source:
    https://docs.python.org/3/library/pickle.html#dispatch-tables
    """

    dispatch_table = copyreg.dispatch_table.copy()  # type: ignore
    dispatch_table[ca.MX] = reduce_MX
    dispatch_table[ca.Function] = reduce_MX


def show_html_from_dataframe(dataframe: pd.DataFrame):
    tmp_file = tempfile.NamedTemporaryFile("w", delete=False)
    min_max_html = (
        dataframe.style.highlight_max(color="green")
        .highlight_min(color="red")
        .highlight_null(null_color="yellow")
        .to_html()
    )
    tmp_file.writelines(min_max_html)
    descibe_html = dataframe.describe().style.to_html()
    tmp_file.writelines(descibe_html)
    tmp_file.seek(0)
    tmp_file.close()
    webbrowser.open("file://" + tmp_file.name)


def plot_array(array, xticks=None, yticks=None):
    """Plots given array in an observable way."""
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt

    div_norm = colors.TwoSlopeNorm(vcenter=0)
    plt.close()
    plt.imshow(array, cmap=cm.coolwarm, norm=div_norm)
    if xticks is not None:
        plt.xticks(range(0, array.shape[1]), xticks)
    if yticks is not None:
        plt.yticks(range(0, array.shape[0]), yticks)
    # plt.colorbar()
    plt.show()


def plot_arrays(arrays):
    """Same as plot_array but plots every array in a subplot."""
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt

    plt.close()
    num_plots = len(arrays)
    fig = plt.figure()

    count_plots = 1
    for array in arrays:
        fig.add_subplot(1, num_plots, count_plots).imshow(array, cmap=cm.Greens_r)
        count_plots = count_plots + 1

    plt.show()


def generate_hammersley(D, N):
    """generate returns a set N points for a D-dimensional Hammersley sequence
    the interval (0,1). Taken from Erik / not tested
    """

    S = np.empty(shape=[N, D], dtype=float)
    prime = prime_class()

    """last column is simple"""
    # S[:,-1] = np.arange(1.,N+1) / N - 1./(2*N)
    S[:, -1] = np.arange(0.0, N) / N

    pn = 2.0 * D
    """getting list of prime numbers"""
    p = prime.List(pn)

    while len(p) < D:
        pn = 2 * pn
        p = prime.List(pn)

    P = p[0 : D - 1]  # last dimension is already set
    """loop for dimensions"""
    for k in range(0, D - 1):
        pk = P[k]
        """loop for hammersley points"""
        for j in range(0, N):
            bj = j + 1
            """maximum for the devision of binary logarithms"""
            n = np.int(np.max([1, np.round(np.log2(bj + 1) / np.log2(pk))]))
            while pk**n <= bj:
                n = n + 1

            b = np.zeros(n)
            b[n - 1] = bj % pk
            while bj and n > 1:
                n = n - 1
                bj = np.floor_divide(bj, pk)
                b[n - 1] = bj % pk

            S[j, k] = np.sum(b[::-1] / pk ** np.arange(1.0, len(b) + 1))
    return S


def make_startpoints(bound0, N, sampling="lhs"):
    """bound0 the boundaries for all sampling points, where the number of tuples gives the number of dimensions D
    N is the number of sampling points
    bound0 = np.array([[0, 10],[0, 100]])
    output = B[num_of_samples, num_of_variables], example B[0] would return an array of variables guesses for all variables
    Taken from Erik, not tested"""
    import pyDOE

    D = len(
        bound0[
            :,
        ]
    )

    if sampling == "lhs":
        S = pyDOE.lhs(D, N)
    elif sampling == "hammersley":
        S = generate_hammersley(D, N)

    B = np.zeros(S.shape)

    for i in range(
        len(
            bound0[
                :,
            ]
        )
    ):
        B[:, i] = S[:, i] * (bound0[i, 1] - bound0[i, 0]) + bound0[i, 0]

    return B


class prime_class:
    """Taken from Erik, not tested"""

    def isEven(self, n):
        return n % 2 == 0

    def isPrime(self, n):
        Dmax = np.sqrt(n)
        if n == 2:
            return True
        if self.isEven(n):
            return False
        d = 3
        while n % d != 0 and d <= Dmax:
            d += 2
        return d > Dmax

    def List(self, nMax):
        return [n for n in np.arange(2, nMax) if self.isPrime(n)]
