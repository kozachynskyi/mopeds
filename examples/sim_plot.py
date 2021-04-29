import copy
import numpy as np

import par_est
import par_est.examples

if __name__ == "__main__":

    piecewiseswitch = True
    variable_list, m = par_est.examples.cstr_dae(piecewiseswitch)
    if piecewiseswitch:
        T_in = variable_list["e0_T_in"]
        T_in.expand_horizon([4000, 7500], [283, 400])

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 40)
    time_grid = np.insert(time_grid, 0, 0)

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed, use_idas_constraints=True)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data(algebraic=True)
    res.plot_states(algebraic=True)
    # np.savetxt("exp.txt", res.toarray().T, delimiter="\t")
