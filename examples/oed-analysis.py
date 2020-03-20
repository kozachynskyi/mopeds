import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est
import cstr


if __name__ == "__main__":

    variable_list, m = cstr.initialize_problem()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    variable_list["e0_E_r1"].fixed = True
    variable_list["e0_E_r2"].fixed = True
    variable_list["e0_E_r3"].fixed = True
    variable_list["e0_k_pre_r1"].fixed = True
    variable_list["e0_k_pre_r2"].fixed = True
    variable_list["e0_k_pre_r3"].fixed = True
    variable_list["e0_U"].fixed = True
    variable_list["e0_c_p"].fixed = True
    # variable_list["e0_greek_Deltah_r1"].fixed = True
    # variable_list["e0_greek_Deltah_r2"].fixed = True
    # variable_list["e0_greek_Deltah_r3"].fixed = True

    # variable_list["e0_c_in_i1"].fixed = True
    # variable_list["e0_c_in_i2"].fixed = True
    # variable_list["e0_c_in_i3"].fixed = True
    # variable_list["e0_c_in_i4"].fixed = True
    # variable_list["e0_T_in"].fixed = True
    # variable_list["e0_T_j"].fixed = True
    # variable_list["e0_F"].fixed = True

    oed = par_est.OptimalExperimentalDesign(m, [variable_list], time_grid)
    # oed.optimize()
    # oed.optimize(False)

    sens_m, fim_m = oed.get_fim_matrix()
    sens_m = ca.fabs(sens_m.toarray())
    fim_m = fim_m.toarray()

    turn_off_states = np.array([1, 1, 1, 1, 1])
    scaling_states = [1, 0.01, 0.01, 0.01, 0.01]

    sc = np.diagflat(np.tile(turn_off_states, len(time_grid) - 1))
    scaling_states_full = np.diagflat(np.tile(turn_off_states / scaling_states, len(time_grid) - 1))

    # Generate scaling for parameters based on their values
    sc_params = []
    for var in variable_list.values():
        if isinstance(var, par_est.Parameter_variable):
            if var.fixed is False:
                sc_params.append(var.value)
    scaling_param_full = np.diagflat(sc_params)

    sens_scaled_full = scaling_states_full @ (sens_m @ scaling_param_full)
    # sc = np.tile(sc, 2)
    fig = plt.figure()
    fig.add_subplot(151).imshow(sens_m, cmap=cm.Greens_r)
    fig.add_subplot(152).imshow(ca.inv(fim_m), cmap=cm.Greens_r)
    fig.add_subplot(153).imshow(sens_scaled_full, cmap=cm.Greens_r)
    plt.show()
