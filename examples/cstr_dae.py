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

    variable_list.add_variable(par_est.Algebraic_variable("e0_c_tot", 13.0))

    variable_list.add_variable(par_est.Parameter_variable("e0_E_r1", 9.6e4, 9.0e4, 10.0e4))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r2", 7.2e4, 6.8e4, 7.6e4))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r3", 6.9e4, 6.5e4, 7.3e4))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_pre_r1", 5.0e6, 4.5e6, 5.5e6))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_pre_r2", 1.0e7, 0.5e7, 1.5e7))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_pre_r3", 5.0e5, 4.5e5, 5.5e5))
    variable_list.add_variable(par_est.Parameter_variable("e0_U", 1.4, 1.0, 1.8))
    variable_list.add_variable(par_est.Parameter_variable("e0_c_p", 3.5, 3.0, 4.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_Deltah_r1", 4.5e-3, 4.0e-3, 5.0e-3))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_Deltah_r2", -5.5e-3, -6.0e-3, -5.0e-3))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_Deltah_r3", 4.5e-3, 4.0e-3, 5.0e-3))

    variable_list.add_variable(par_est.Control_variable("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.Control_variable("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.Control_variable("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.Control_variable("e0_c_in_i4", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.Control_variable("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.Control_variable("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.Control_variable("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))
    # fmt: on

    variable_list["e0_E_r1"].guess = variable_list["e0_E_r1"].lower_bound
    variable_list["e0_E_r2"].guess = variable_list["e0_E_r2"].lower_bound
    variable_list["e0_E_r3"].guess = variable_list["e0_E_r3"].lower_bound
    variable_list["e0_k_pre_r1"].guess = variable_list["e0_k_pre_r1"].lower_bound
    variable_list["e0_k_pre_r2"].guess = variable_list["e0_k_pre_r2"].lower_bound
    variable_list["e0_k_pre_r3"].guess = variable_list["e0_k_pre_r3"].lower_bound
    variable_list["e0_U"].guess = variable_list["e0_U"].lower_bound
    variable_list["e0_c_p"].guess = variable_list["e0_c_p"].lower_bound
    variable_list["e0_greek_Deltah_r1"].guess = variable_list[
        "e0_greek_Deltah_r1"
    ].lower_bound
    variable_list["e0_greek_Deltah_r2"].guess = variable_list[
        "e0_greek_Deltah_r2"
    ].lower_bound
    variable_list["e0_greek_Deltah_r3"].guess = variable_list[
        "e0_greek_Deltah_r3"
    ].lower_bound

    variable_list["e0_c_in_i1"].guess = variable_list["e0_c_in_i1"].lower_bound
    variable_list["e0_c_in_i2"].guess = variable_list["e0_c_in_i2"].lower_bound
    variable_list["e0_c_in_i3"].guess = variable_list["e0_c_in_i3"].lower_bound
    variable_list["e0_c_in_i4"].guess = variable_list["e0_c_in_i4"].lower_bound
    variable_list["e0_T_in"].guess = variable_list["e0_T_in"].lower_bound
    variable_list["e0_T_j"].guess = variable_list["e0_T_j"].lower_bound
    variable_list["e0_F"].guess = variable_list["e0_F"].lower_bound

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_T_in"].casadi_var - m._all_variables["e0_T"].casadi_var))) + (((m._all_variables["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (m._all_variables["e0_c_p"].casadi_var * e0_V))) * ((m._all_variables["e0_T_j"].casadi_var - m._all_variables["e0_T"].casadi_var)))) + (((-m._all_variables["e0_greek_Deltah_r1"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-m._all_variables["e0_greek_Deltah_r2"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-m._all_variables["e0_greek_Deltah_r3"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c1dot = ((((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i1"].casadi_var - m._all_variables["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c2dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i2"].casadi_var - m._all_variables["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c3dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i3"].casadi_var - m._all_variables["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c4dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i4"].casadi_var - m._all_variables["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))

    c_tot = m._all_variables["e0_c_tot"].casadi_var - m._all_variables["e0_c_i1"].casadi_var - m._all_variables["e0_c_i2"].casadi_var - m._all_variables["e0_c_i3"].casadi_var - m._all_variables["e0_c_i4"].casadi_var

    # fmt: on

    m.add_differential_equations([tdot, c1dot, c2dot, c3dot, c4dot])
    m.add_algebraic_equations([c_tot])

    return variable_list, m


if __name__ == "__main__":

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
    var_list_exp = sim_fixed.generate_exp_data()
    var_list_exp.plot_states()
