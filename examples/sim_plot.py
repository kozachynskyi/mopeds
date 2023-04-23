import copy

import numpy as np
from matplotlib import pyplot as plt

import par_est
import par_est.examples

plt.ion()
import pandas as pd

if __name__ == "__main__":

    piecewiseswitch = False
    variable_list, m = par_est.examples.cstr_dae(piecewiseswitch)

    T_in = variable_list["e0_T_in"]
    if isinstance(T_in, par_est.VariableControlPiecewiseConstant):
        T_in.expand_horizon([4000, 7500], [283, 400])

    # Create time-grid. Zero should be first
    time_grid = np.linspace(1000, 10000, 2)
    time_grid = np.insert(time_grid, 0, 0)

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = par_est.Simulator(
        m, time_grid, var_list_fixed, use_idas_constraints=True, simulate_jac=True
    )
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate_sym()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data(algebraic=True)
    print(res.dataframe)
    # 1970-01-01 01:06:40.000000000

    a = sim_fixed.simulate_jac()["jac_xf_p"]
    b = pd.DataFrame(a.toarray(), index=list(sim_fixed.model.varlist_state.keys()))
    breakpoint()
    # res["e0_c_i1"].plot()
    # res.plot_new()
    # res_res.plot_new()
    # for var in res.values():
    #     a = var.value.equals(res_res[var.name].value)
    #     # breakpoint()
    #     # diff = var.value - res_res[var.name].value
    #     # b = diff.isin([0]).all(axis=None)
    #     if not a:
    #         aa = var.value
    #         bb = res_res[var.name].value
    #         breakpoint()
    # var.value.compare(res_res[var.name].value)
    # res.plot_states(algebraic=True, as_one_plot=False)
    # np.savetxt("exp.txt", res.toarray().T, delimiter="\t")
