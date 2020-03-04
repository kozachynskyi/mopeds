import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est


def initialize_problem():
    e0_greek_nu_i1_r1 = -1.0
    e0_greek_nu_i1_r2 = 1.0
    e0_greek_nu_i2_r2 = -1.0
    e0_greek_nu_i3_r1 = 1.0
    e0_greek_nu_i1_r3 = -1.0
    e0_greek_nu_i4_r3 = 1.0
    e0_greek_rho = 800.0
    e0_A = 1.0
    e0_R = 8.314
    e0_V = 1.0

    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.State_variable("e0_T", 273.0, 10))
    variable_list.add_variable(par_est.State_variable("e0_c_i1", 3.0, 20))
    variable_list.add_variable(par_est.State_variable("e0_c_i2", 10.0, 30))
    variable_list.add_variable(par_est.State_variable("e0_c_i3", 0.0, 40))
    variable_list.add_variable(par_est.State_variable("e0_c_i4", 0.0, 50))

    variable_list.add_variable(par_est.Parameter_variable("e0_E_r1", 9.6e4, 9.0e4, 10.0e4))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r2", 7.2e4, 6.8e4, 7.6e4))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r3", 6.9e4, 6.5e4, 7.3e4))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_pre_r1", 5.0e6, 4.5e6, 5.5e6))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_pre_r2", 1.0e7, 0.5e7, 1.5e7))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_pre_r3", 5.0e5, 4.5e5, 5.5e5))
    variable_list.add_variable(par_est.Parameter_variable("e0_U", 1.4, 1.0, 1.8))
    variable_list.add_variable(par_est.Parameter_variable("e0_c_p", 3.5, 3.0, 4.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_Deltah_r1", 4.5e-3, 4.0e-3, 5.0e-3))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_Deltah_r2", -5.5e-3, -5.0e-3, -6.0e-3))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_Deltah_r3", 4.5e-3, 4.0e-3, 5.0e-3))

    variable_list.add_variable(par_est.Control_variable("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.Control_variable("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.Control_variable("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.Control_variable("e0_c_in_i4", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.Control_variable("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.Control_variable("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.Control_variable("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))
    # fmt: on

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_T_in"].casadi_var - m._all_variables["e0_T"].casadi_var))) + (((m._all_variables["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (m._all_variables["e0_c_p"].casadi_var * e0_V))) * ((m._all_variables["e0_T_j"].casadi_var - m._all_variables["e0_T"].casadi_var)))) + (((-m._all_variables["e0_greek_Deltah_r1"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-m._all_variables["e0_greek_Deltah_r2"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-m._all_variables["e0_greek_Deltah_r3"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c1dot = ((((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i1"].casadi_var - m._all_variables["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c2dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i2"].casadi_var - m._all_variables["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c3dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i3"].casadi_var - m._all_variables["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c4dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i4"].casadi_var - m._all_variables["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    # fmt: on

    m.add_equations([tdot, c1dot, c2dot, c3dot, c4dot])

    return variable_list, m


variable_list, m = initialize_problem()

# Create time-grid. Zero should be first
time_grid = np.linspace(10, 10000, 40)
time_grid = np.insert(time_grid, 0, 0)

# Generate experimental data
var_list_fixed = copy.deepcopy(variable_list)
for var in var_list_fixed.values():
    var.fixed = True

sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
res = sim_fixed.generate_exp_data()
# res.plot_states()
# np.savetxt("exp.txt", res.toarray().T, delimiter="\t")

ss = par_est.Simulator(m, time_grid, variable_list)
var_list_exp = sim_fixed.generate_exp_data()

# Replace empty state variables with results from simulation
var_list2 = copy.deepcopy(variable_list)
for key, var in var_list_exp.items():
    var_list2[key] = var

# var_list2["e0_T"].value = par_est.Experimental_Data()
# var_list2["e0_c_i4"].value = par_est.Experimental_Data()

var_list2["e0_E_r1"].fixed = False
var_list2["e0_E_r1"].guess = var_list2["e0_E_r1"].lower_bound
var_list2["e0_E_r2"].fixed = False
var_list2["e0_E_r2"].guess = var_list2["e0_E_r2"].lower_bound
var_list2["e0_E_r3"].fixed = False
var_list2["e0_E_r3"].guess = var_list2["e0_E_r3"].lower_bound

var_list2["e0_k_pre_r1"].fixed = False
var_list2["e0_k_pre_r1"].guess = var_list2["e0_k_pre_r1"].lower_bound
var_list2["e0_k_pre_r2"].fixed = False
var_list2["e0_k_pre_r2"].guess = var_list2["e0_k_pre_r2"].lower_bound
var_list2["e0_k_pre_r3"].fixed = False
var_list2["e0_k_pre_r3"].guess = var_list2["e0_k_pre_r3"].lower_bound

var_list2["e0_U"].fixed = False
var_list2["e0_U"].guess = var_list2["e0_U"].lower_bound
var_list2["e0_c_p"].fixed = False
var_list2["e0_c_p"].guess = var_list2["e0_c_p"].lower_bound

# var_list2["e0_greek_Deltah_r1"].fixed = False
# var_list2["e0_greek_Deltah_r1"].guess = var_list2["e0_greek_Deltah_r1"].lower_bound
# var_list2["e0_greek_Deltah_r2"].fixed = False
# var_list2["e0_greek_Deltah_r2"].guess = var_list2["e0_greek_Deltah_r2"].lower_bound
# var_list2["e0_greek_Deltah_r3"].fixed = False
# var_list2["e0_greek_Deltah_r3"].guess = var_list2["e0_greek_Deltah_r3"].lower_bound

var_list2["e0_c_in_i1"].fixed = False
var_list2["e0_c_in_i1"].guess = var_list2["e0_c_in_i1"].lower_bound
var_list2["e0_c_in_i2"].fixed = False
var_list2["e0_c_in_i2"].guess = var_list2["e0_c_in_i2"].lower_bound
var_list2["e0_c_in_i3"].fixed = False
var_list2["e0_c_in_i3"].guess = var_list2["e0_c_in_i3"].lower_bound
var_list2["e0_c_in_i4"].fixed = False
var_list2["e0_c_in_i4"].guess = var_list2["e0_c_in_i4"].lower_bound
var_list2["e0_T_in"].fixed = False
var_list2["e0_T_in"].guess = var_list2["e0_T_in"].lower_bound
var_list2["e0_T_j"].fixed = False
var_list2["e0_T_j"].guess = var_list2["e0_T_j"].lower_bound
var_list2["e0_F"].fixed = False
var_list2["e0_F"].guess = var_list2["e0_F"].lower_bound

# start_time = datetime(2018, 1, 1, 1, 0, 0, 0) + timedelta(days=1)
# end_time = start_time + timedelta(seconds=1000)
var_list_oed = copy.deepcopy(var_list2)
# var_list3 = copy.deepcopy(var_list2)
# var_list_exp.write_data_opcua(start_time)
# var_list3.get_data_opcua(start_time, end_time)
pe = par_est.ParameterEstimation(m, [var_list2, var_list2])
# pe.optimize()

oed = par_est.OptimalExperimentalDesign(m, [var_list_oed], time_grid)
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
oed.optimize()
# pe.optimize(False)
