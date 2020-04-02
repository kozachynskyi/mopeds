# -*- coding: utf-8 -*-
"""
    Dummy conftest.py for par_est.

    If you don't know what this is for, just leave it empty.
    Read more about conftest.py under:
    https://pytest.org/latest/plugins.html
"""

import casadi as ca

import par_est


def pendulum_dae_1():
    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.VariableState("x", 3.0))
    variable_list.add_variable(par_est.VariableState("u", -1.0 / 3))

    variable_list.add_variable(par_est.VariableAlgebraic("y", 4.0))
    variable_list.add_variable(par_est.VariableAlgebraic("v", 1.0 / 4))
    variable_list.add_variable(par_est.VariableAlgebraic("lambda", 1147.0 / 720))

    variable_list.add_variable(par_est.VariableControl("L", 5.0))
    variable_list.add_variable(par_est.VariableParameter("g", 10.0))
    # fmt: on

    m = par_est.Model(variable_list)

    # fmt: off
    dydx1 = m._all_variables["u"].casadi_var
    dydx2 = m._all_variables["lambda"].casadi_var * m._all_variables["x"].casadi_var

    alg1 = m._all_variables["x"].casadi_var ** 2 + m._all_variables["y"].casadi_var ** 2 - m._all_variables["L"].casadi_var ** 2
    alg2 = m._all_variables["u"].casadi_var * m._all_variables["x"].casadi_var + m._all_variables["v"].casadi_var * m._all_variables["y"].casadi_var
    alg3 = m._all_variables["u"].casadi_var ** 2 - m._all_variables["g"].casadi_var * m._all_variables["y"].casadi_var + m._all_variables["v"].casadi_var ** 2 + m._all_variables["L"].casadi_var ** 2 * m._all_variables["lambda"].casadi_var

    m.add_differential_equations([dydx1, dydx2, ])
    m.add_algebraic_equations([alg1, alg2, alg3, ])
    # fmt: on

    return variable_list, m


def cstr_model_ode():
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
    variable_list.add_variable(par_est.VariableState("e0_T", 273.0, 10))
    variable_list.add_variable(par_est.VariableState("e0_c_i1", 3.0, 20))
    variable_list.add_variable(par_est.VariableState("e0_c_i2", 10.0, 30))
    variable_list.add_variable(par_est.VariableState("e0_c_i3", 0.0, 40))
    variable_list.add_variable(par_est.VariableState("e0_c_i4", 0.0, 50))

    variable_list.add_variable(par_est.VariableParameter("e0_E_r1", 9.6e4, 9.0e4, 10.0e4))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r2", 7.2e4, 6.8e4, 7.6e4))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r3", 6.9e4, 6.5e4, 7.3e4))
    variable_list.add_variable(par_est.VariableParameter("e0_k_pre_r1", 5.0e6, 4.5e6, 5.5e6))
    variable_list.add_variable(par_est.VariableParameter("e0_k_pre_r2", 1.0e7, 0.5e7, 1.5e7))
    variable_list.add_variable(par_est.VariableParameter("e0_k_pre_r3", 5.0e5, 4.5e5, 5.5e5))
    variable_list.add_variable(par_est.VariableParameter("e0_U", 1.4, 1.0, 1.8))
    variable_list.add_variable(par_est.VariableParameter("e0_c_p", 3.5, 3.0, 4.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_Deltah_r1", 4.5e-3, 4.0e-3, 5.0e-3))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_Deltah_r2", -5.5e-3, -6.0e-3, -5.0e-3))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_Deltah_r3", 4.5e-3, 4.0e-3, 5.0e-3))

    variable_list.add_variable(par_est.VariableControl("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i4", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.VariableControl("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))
    # fmt: on

    for var in variable_list.values():
        var.guess = var.lower_bound
    # variable_list["e0_E_r1"].guess = variable_list["e0_E_r1"].lower_bound
    # variable_list["e0_E_r2"].guess = variable_list["e0_E_r2"].lower_bound
    # variable_list["e0_E_r3"].guess = variable_list["e0_E_r3"].lower_bound
    # variable_list["e0_k_pre_r1"].guess = variable_list["e0_k_pre_r1"].lower_bound
    # variable_list["e0_k_pre_r2"].guess = variable_list["e0_k_pre_r2"].lower_bound
    # variable_list["e0_k_pre_r3"].guess = variable_list["e0_k_pre_r3"].lower_bound
    # variable_list["e0_U"].guess = variable_list["e0_U"].lower_bound
    # variable_list["e0_c_p"].guess = variable_list["e0_c_p"].lower_bound
    # variable_list["e0_greek_Deltah_r1"].guess = variable_list[
    #     "e0_greek_Deltah_r1"
    # ].lower_bound
    # variable_list["e0_greek_Deltah_r2"].guess = variable_list[
    #     "e0_greek_Deltah_r2"
    # ].lower_bound
    # variable_list["e0_greek_Deltah_r3"].guess = variable_list[
    #     "e0_greek_Deltah_r3"
    # ].lower_bound

    # variable_list["e0_c_in_i1"].guess = variable_list["e0_c_in_i1"].lower_bound
    # variable_list["e0_c_in_i2"].guess = variable_list["e0_c_in_i2"].lower_bound
    # variable_list["e0_c_in_i3"].guess = variable_list["e0_c_in_i3"].lower_bound
    # variable_list["e0_c_in_i4"].guess = variable_list["e0_c_in_i4"].lower_bound
    # variable_list["e0_T_in"].guess = variable_list["e0_T_in"].lower_bound
    # variable_list["e0_T_j"].guess = variable_list["e0_T_j"].lower_bound
    # variable_list["e0_F"].guess = variable_list["e0_F"].lower_bound

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_T_in"].casadi_var - m._all_variables["e0_T"].casadi_var))) + (((m._all_variables["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (m._all_variables["e0_c_p"].casadi_var * e0_V))) * ((m._all_variables["e0_T_j"].casadi_var - m._all_variables["e0_T"].casadi_var)))) + (((-m._all_variables["e0_greek_Deltah_r1"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-m._all_variables["e0_greek_Deltah_r2"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-m._all_variables["e0_greek_Deltah_r3"].casadi_var) / (e0_greek_rho * m._all_variables["e0_c_p"].casadi_var)) * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c1dot = ((((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i1"].casadi_var - m._all_variables["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c2dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i2"].casadi_var - m._all_variables["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-m._all_variables["e0_E_r2"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c3dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i3"].casadi_var - m._all_variables["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r1"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c4dot = ((m._all_variables["e0_F"].casadi_var / e0_V) * ((m._all_variables["e0_c_in_i4"].casadi_var - m._all_variables["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-m._all_variables["e0_E_r3"].casadi_var) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    # fmt: on

    m.add_differential_equations([tdot, c1dot, c2dot, c3dot, c4dot])

    return variable_list, m
