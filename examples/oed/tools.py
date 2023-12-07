import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import copy

import mopeds

def unfix_parameters(varlist):
    unfix_variables = ["theta1", "theta2", "theta3","theta4"]
    varlist_oed = copy.deepcopy(varlist)

    for par_name in unfix_variables:
        varlist_oed[par_name].fixed = False
    return varlist_oed

def plot_res(pe, parameters):
    data = pe.calculate_objective_and_residual(parameters, "wls")
    y_all = data["y"]
    y_m = pe.array_data

    y_0 = pe.list_simulators[0]._initial_state
    x = pe.list_simulators[0].time_grid_relative

    y_all = np.insert(y_all, 0, y_0.T, 0)
    y_m = np.insert(y_m, 0, y_0.T, 0)
    fig, ax = plt.subplots(1, 2)
    for y, y_m, axis in zip(y_all.T, y_m.T, ax):
        axis.plot(x, y)
        axis.scatter(x, y_m)
    fig.show()

def plot_fig3(oed, obj_f_name):
    varlist = oed.list_input_varlist[0]
    control_bounds = {}
    for var_name in ["u1", "u2"]:
        lb = varlist[var_name].lower_bound
        ub = varlist[var_name].upper_bound
        num_p = 15
        control_bounds[var_name] = [lb, ub, num_p]

    grid, meshgrid = mopeds.tools.create_grid(control_bounds.values())

    obj_f = []
    for grid_i in grid:
        controls = {"u1": grid_i[0], "u2": grid_i[1]}
        res = oed.calculate_objective_and_jacobian(controls, obj_f_name)
        obj_f.append(res["f"])
    z = np.reshape(obj_f, meshgrid[0].shape)
    plt.contourf(meshgrid[0], meshgrid[1], z, locator=ticker.LogLocator())
    plt.colorbar()
    plt.show()
    print(obj_f)
