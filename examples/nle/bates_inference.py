"""This are examples from Bates book, showing how to plot
model uncertainty that is caused by parameter uncertainty"""

import matplotlib.pyplot as plt

import mopeds
import numpy as np


def plot_data(exp_data, label):
    plt.plot(exp_data["x"], exp_data["f"], "ko", label=label)


def plot_inference(exp_inference_results):
    plt.plot(
        exp_inference_results["x"],
        exp_inference_results["f"]["simulation"],
        "k-",
        label="Sim. data",
    )
    plt.plot(
        exp_inference_results["x"],
        exp_inference_results["f"]["lower bound"],
        "b:",
        label="Lower bound",
    )
    plt.plot(
        exp_inference_results["x"],
        exp_inference_results["f"]["upper bound"],
        "r:",
        label="Upper bound",
    )


def set_plot_legend(x_bounds, y_bounds, x_label, y_label):
    plt.legend()
    plt.xlim(*x_bounds)
    plt.xlabel(x_label)
    plt.ylim(y_bounds)
    plt.ylabel(y_label)
    plt.grid()
    plt.show()


def bod():
    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.bod_model()

    dict_of_params = {
        "theta1": 19.143,
        "theta2": 0.5311,
    }

    x_bounds = [0.0, 8.0]
    y_bounds = [-5.0, 30.0]
    dict_of_controls = {
        "x": [x_bounds[0], x_bounds[1], 90 + 1],
    }

    dict_of_responses = {
        "f": 1e1,
    }

    pe = mopeds.ParameterEstimationNLE(MODEL, EXP_DATA)
    exp_inference_results, exp_data, sim_data = pe.calculate_inference_bounds(
        dict_of_params, dict_of_responses, dict_of_controls
    )

    plot_data(exp_data, "Exp. data")
    plot_inference(exp_inference_results)
    set_plot_legend(x_bounds, y_bounds, "Time", "Oxygen Demand")

    dict_of_artificial_controls = {
        "x": [x_bounds[0], x_bounds[1], 10 + 1],
    }
    (
        artificial_inference_results,
        artificial_data,
        artificial_data_sim,
    ) = pe.calculate_inference_bounds(
        dict_of_params, dict_of_responses, dict_of_controls, dict_of_artificial_controls
    )

    plot_data(artificial_data, "Artificial data")
    plot_inference(artificial_inference_results)
    set_plot_legend(x_bounds, y_bounds, "Time", "Oxygen Demand")

    print(np.round(exp_inference_results["f"]["s"], 2))
    print(np.round(artificial_inference_results["f"]["s"], 2))


def puromycin():
    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.puromycin_model()

    data = EXP_DATA["Treated"]
    list_of_params = list(["theta1", "theta2"])

    dict_of_params = {}
    for param in list_of_params:
        dict_of_params[param] = float(data[0][param].value[0])

    x_bounds = [0.0, 1.2]
    y_bounds = [0.0, 250]

    dict_of_controls = {
        "x": [x_bounds[0], x_bounds[1], 120 + 1],
    }

    dict_of_responses = {
        "f": 5e1,
    }

    pe = mopeds.ParameterEstimationNLE(MODEL, data)
    exp_inference_results, exp_data, sim_data = pe.calculate_inference_bounds(
        dict_of_params, dict_of_responses, dict_of_controls
    )

    plot_data(exp_data, "Exp. data")
    plot_inference(exp_inference_results)
    set_plot_legend(x_bounds, y_bounds, "Concentration", "Velocity")

    dict_of_artificial_controls = {
        "x": [x_bounds[0], x_bounds[1], 10 + 1],
    }
    (
        artificial_inference_results,
        artificial_data,
        artificial_sim_data,
    ) = pe.calculate_inference_bounds(
        dict_of_params,
        dict_of_responses,
        dict_of_controls,
        dict_of_artificial_controls,
        rng=np.random.default_rng(42),
    )

    plot_data(artificial_data, "Artificial data")
    plot_inference(artificial_inference_results)
    set_plot_legend(x_bounds, y_bounds, "Concentration", "Velocity")

    print(np.round(exp_inference_results["f"]["s"], 2))
    print(np.round(artificial_inference_results["f"]["s"], 2))


if __name__ == "__main__":
    # bod()
    puromycin()
