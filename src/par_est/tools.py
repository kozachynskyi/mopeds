from matplotlib import pyplot as plt
import matplotlib.cm as cm


def plot_array(array):
    """ Plots given array in an observable way. """
    plt.close()
    plt.imshow(array, cmap=cm.Greens_r)
    plt.show()


def plot_arrays(arrays):
    """ Same as plot_array but plots every array in a subplot. """
    plt.close()
    num_plots = len(arrays)
    fig = plt.figure()

    count_plots = 1
    for array in arrays:
        fig.add_subplot(1, num_plots, count_plots).imshow(array, cmap=cm.Greens_r)
        count_plots = count_plots + 1

    plt.show()
