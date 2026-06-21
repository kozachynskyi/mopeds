import copy

import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples

plt.ion()
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    variable_list, m, exp_data = mopeds.examples.free_fall_example()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(0, 30, 5)
    orig = copy.deepcopy(exp_data[0]["s"].dataframe)
    exp_data[0]["s"].dataframe.iloc[-1] = 100

    # Create simulation Object
    sim_idas = mopeds.Simulator(
        m,
        time_grid,
        variable_list,
        use_idas_constraints=True,
        simulate_jac=True,
        integrator_name="idas",
    )

    integrator_settings = {}
    integrator_settings["acados"] = {
        "integrator_type": "IRK",
        "collocation_type": "GAUSS_RADAU_IIA",
        "num_stages": 3,
        "num_steps": 10,
        "newton_tol": 1e-8,
        "newton_iter": 100,
        "code_reuse": True,
    }
    r = []

    pe = mopeds.ParameterEstimation(
        m,
        exp_data,
        simulator_name="acados",
        simulator_settings=integrator_settings,
    )
    par = pe.optimize(True, objective_function="fair")["x_dict"]

    res = pe.calculate_objective_and_residual(par)
    plt.plot(pe.list_simulators[0].time_grid_relative[1:], res["y"])
    plt.scatter(
        pe.list_simulators[0].time_grid_relative[1:],
        exp_data[0].dataframe["s"].iloc[1:],
    )
    plt.scatter(pe.list_simulators[0].time_grid_relative[1:], orig["s"].iloc[1:])
    plt.show()
    # print(sim_fixed.simulate_jac())
    # for sim in [sim_acados, sim_idas]:
    #     print(res_simple)
