import copy

import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples

plt.ion()
import pandas as pd

if __name__ == "__main__":

    piecewiseswitch = False
    variable_list, m = mopeds.examples.cstr_dae(piecewiseswitch)

    T_in = variable_list["e0_T_in"]
    if isinstance(T_in, mopeds.VariableControlPiecewiseConstant):
        T_in.expand_horizon([4000, 7500], [283, 400])

    # Create time-grid. Zero should be first
    time_grid = np.linspace(0, 100, 2)
    # time_grid = np.insert(time_grid, 0, 0)

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # var_list_fixed["e0_U"].fixed = False
    var_list_fixed["e0_T"].fixed = False
    # Create simulation Object
    sim_fixed = mopeds.Simulator(
        m, time_grid, var_list_fixed, use_idas_constraints=True, simulate_jac=False
    )
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    print(sim_fixed.simulate_sym())
    print(sim_fixed.simulate_sym_unfixed([273]))
