import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est
import cstr

if __name__ == "__main__":

    variable_list, m = cstr.initialize_problem_ode()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    # Generate experimental data for parameter estimation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).generate_exp_data()

    # Replace empty state variables with results from simulation
    variable_list_optimizer = copy.deepcopy(variable_list)
    for key, var in var_list_exp.items():
        variable_list_optimizer[key] = var

    # variable_list_optimizer["e0_T"].value = par_est.Experimental_Data()

    # variable_list_optimizer["e0_E_r1"].fixed = True
    # variable_list_optimizer["e0_E_r2"].fixed = True
    # variable_list_optimizer["e0_E_r3"].fixed = True
    # variable_list_optimizer["e0_k_pre_r1"].fixed = True
    # variable_list_optimizer["e0_k_pre_r2"].fixed = True
    # variable_list_optimizer["e0_k_pre_r3"].fixed = True
    # variable_list_optimizer["e0_U"].fixed = True
    # variable_list_optimizer["e0_c_p"].fixed = True
    # variable_list_optimizer["e0_greek_Deltah_r1"].fixed = True
    # variable_list_optimizer["e0_greek_Deltah_r2"].fixed = True
    # variable_list_optimizer["e0_greek_Deltah_r3"].fixed = True

    # variable_list_optimizer["e0_c_in_i1"].fixed = True
    # variable_list_optimizer["e0_c_in_i2"].fixed = True
    # variable_list_optimizer["e0_c_in_i3"].fixed = True
    # variable_list_optimizer["e0_c_in_i4"].fixed = True
    # variable_list_optimizer["e0_T_in"].fixed = True
    # variable_list_optimizer["e0_T_j"].fixed = True
    # variable_list_optimizer["e0_F"].fixed = True

    pe = par_est.ParameterEstimation(m, [variable_list_optimizer, variable_list_optimizer])
    pe.optimize()
    pe.optimize(False)

    # variable_list_optimizer["e0_E_r1"].fixed = True
    # variable_list_optimizer["e0_E_r2"].fixed = True
    # variable_list_optimizer["e0_E_r3"].fixed = True
    # variable_list_optimizer["e0_k_pre_r1"].fixed = True
    # variable_list_optimizer["e0_k_pre_r2"].fixed = True
    # variable_list_optimizer["e0_k_pre_r3"].fixed = True
    # variable_list_optimizer["e0_U"].fixed = True
    # variable_list_optimizer["e0_c_p"].fixed = True
    variable_list_optimizer["e0_greek_Deltah_r1"].fixed = True
    variable_list_optimizer["e0_greek_Deltah_r2"].fixed = True
    variable_list_optimizer["e0_greek_Deltah_r3"].fixed = True

    # variable_list_optimizer["e0_c_in_i1"].fixed = True
    # variable_list_optimizer["e0_c_in_i2"].fixed = True
    # variable_list_optimizer["e0_c_in_i3"].fixed = True
    # variable_list_optimizer["e0_c_in_i4"].fixed = True
    # variable_list_optimizer["e0_T_in"].fixed = True
    # variable_list_optimizer["e0_T_j"].fixed = True
    # variable_list_optimizer["e0_F"].fixed = True

    oed = par_est.OptimalExperimentalDesign(m, [variable_list_optimizer], time_grid)
    oed.optimize()
    oed.optimize(False)
    # a = oed.get_fim_matrix()
    # b = a[0].toarray()
    # b = ca.fabs(b)
    # c = a[1].toarray()
    # turn_off_states = np.array([1, 1, 1, 1, 1])
    # sc_states = [1, 0.01, 0.01, 0.01, 0.01]
    # sc = np.diagflat(np.tile(turn_off_states, len(time_grid) - 1))
    # sc_full = np.diagflat(np.tile(turn_off_states / sc_states, len(time_grid) - 1))
    # # num_states = 5
    # # sc_full_params = np.tile(sc_params, ((len(time_grid) - 1) * num_states, 1)).T
    # sc_params = [5000000.0, 10000000.0, 500000.0, 1.4]
    # sc_full_params = np.diagflat(sc_params)
    # b_scaled = sc @ b
    # b_scaled_full = sc_full @ (b @ sc_full_params)
    # # sc = np.tile(sc, 2)
    # fig = plt.figure()
    # fig.add_subplot(151).imshow(b, cmap=cm.Greens_r)
    # fig.add_subplot(152).imshow(ca.inv(c), cmap=cm.Greens_r)
    # plt.show()
