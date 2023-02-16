import copy

import numpy as np
from matplotlib import pyplot as plt

import par_est
import par_est.examples

plt.ion()
import pandas as pd

if __name__ == "__main__":

    variable_list, m, exp_data = par_est.examples.free_fall_example()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(0, 30, 5)

    # Create simulation Object
    sim_idas = par_est.Simulator(
        m, time_grid, variable_list, use_idas_constraints=True, simulate_jac=True, integrator_name="idas",
    )


    integrator_settings = {}
    integrator_settings["acados"] = { 
        "integrator_type": "IRK",
        "collocation_type": "GAUSS_RADAU_IIA",
        "num_stages": 3,
        "num_steps": 10,
        "newton_tol": 1e-8,
        "newton_iter": 100,
        "code_reuse": False,
    }
    r = []

    for i in range(2):
        if i == 0:
            pe = par_est.ParameterEstimation(m, exp_data*3, simulator_name="acados", simulator_settings=integrator_settings)
            par = pe.optimize(True)["x_dict"]
        else:
            pe = par_est.ParameterEstimation(m, exp_data*3, simulator_name="idas")

        r.append(pe.calculate_objective_and_residual(par)["y"])
    print(r[0] - r[1])
    # print(sim_fixed.simulate_jac())
    # for sim in [sim_acados, sim_idas]:
    #     res_simple = sim.simulate_sym()
    #     print(res_simple)
