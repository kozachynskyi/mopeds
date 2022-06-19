from __future__ import annotations

import copy

import casadi as ca

import par_est


def empy_dae(piecewise_control: bool = False) -> tuple[par_est.VariableList, par_est.Model]:
    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.VariableState("X1", 0))

    variable_list.add_variable(par_est.VariableState("X2", 0))

    variable_list.add_variable(par_est.VariableAlgebraic("Z1", 0.0))

    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("C", 0.0, -1, 1))
    else:
        variable_list.add_variable(par_est.VariableControl("C", 0.0, -1, 1))
    variable_list.add_variable(par_est.VariableParameter("P", 0.0, -1, 1))
    # fmt: on

    m = par_est.Model(variable_list)

    # fmt: off
    dydx1 = m.varlist_all["C"].casadi_var * 0
    dydx2 = m.varlist_all["P"].casadi_var * 0
    alg1 = m.varlist_all["X1"].casadi_var + m.varlist_all["X2"].casadi_var + m.varlist_all["Z1"].casadi_var

    m.add_equations_differential([dydx1, dydx2, ])
    m.add_equations_algebraic([alg1, ])
    # fmt: on

    return variable_list, m


def pendulum_dae_1(piecewise_control: bool = False, variable_list: None | par_est.VariableList = None) -> tuple[par_est.VariableList, par_est.Model]:
    if variable_list is None:
        variable_list = par_est.VariableList()

        # fmt: off
        variable_list.add_variable(par_est.VariableState("x", 3.0))
        variable_list.add_variable(par_est.VariableState("u", -1.0 / 3))

        variable_list.add_variable(par_est.VariableAlgebraic("y", 4.0))
        variable_list.add_variable(par_est.VariableAlgebraic("v", 1.0 / 4))
        variable_list.add_variable(par_est.VariableAlgebraic("lambda", 1147.0 / 720))

        if piecewise_control:
            variable_list.add_variable(par_est.VariableControlPiecewiseConstant("L", 5.0))
        else:
            variable_list.add_variable(par_est.VariableControl("L", 5.0))
        variable_list.add_variable(par_est.VariableParameter("g", 10.0))
    # fmt: on

    m = par_est.Model(variable_list)

    # fmt: off
    dydx1 = m.varlist_all["u"].casadi_var
    dydx2 = m.varlist_all["lambda"].casadi_var * m.varlist_all["x"].casadi_var

    alg1 = m.varlist_all["x"].casadi_var ** 2 + m.varlist_all["y"].casadi_var ** 2 - m.varlist_all["L"].casadi_var ** 2
    alg2 = m.varlist_all["u"].casadi_var * m.varlist_all["x"].casadi_var + m.varlist_all["v"].casadi_var * m.varlist_all["y"].casadi_var
    alg3 = m.varlist_all["u"].casadi_var ** 2 - m.varlist_all["g"].casadi_var * m.varlist_all["y"].casadi_var + m.varlist_all["v"].casadi_var ** 2 + m.varlist_all["L"].casadi_var ** 2 * m.varlist_all["lambda"].casadi_var

    m.add_equations_differential([dydx1, dydx2, ])
    m.add_equations_algebraic([alg1, alg2, alg3, ])
    # fmt: on

    return variable_list, m


def cstr_ode(piecewise_control: bool = False) -> tuple[par_est.VariableList, par_est.Model]:
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
    variable_list.add_variable(par_est.VariableState("e0_T", 273.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i1", 3.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i2", 10.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i3", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i4", 0.0))

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

    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_c_in_i1", 5.0, 4.0, 6.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i4", 0.0, 0.0, 1.0))
    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_T_in", 373.0, 353.0, 393.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))
    # fmt: on

    for var in variable_list.values():
        if isinstance(var, (par_est.VariableParameter, par_est.VariableControl)):
            var.guess = var.lower_bound

    if piecewise_control:
        var = variable_list["e0_T_in"].variable_list.index(0)
        var.guess = var.lower_bound
        var = variable_list["e0_c_in_i1"].variable_list.index(0)
        var.guess = var.lower_bound

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_T_in"].casadi_var - m.varlist_all["e0_T"].casadi_var))) + (((m.varlist_all["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (m.varlist_all["e0_c_p"].casadi_var * e0_V))) * ((m.varlist_all["e0_T_j"].casadi_var - m.varlist_all["e0_T"].casadi_var)))) + (((-m.varlist_all["e0_greek_Deltah_r1"].casadi_var) / (e0_greek_rho * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r2"].casadi_var) / (e0_greek_rho * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r3"].casadi_var) / (e0_greek_rho * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c1dot = ((((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i1"].casadi_var - m.varlist_all["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c2dot = ((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i2"].casadi_var - m.varlist_all["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c3dot = ((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i3"].casadi_var - m.varlist_all["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c4dot = ((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i4"].casadi_var - m.varlist_all["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    # fmt: on

    m.add_equations_differential([tdot, c1dot, c2dot, c3dot, c4dot])

    return variable_list, m


def cstr_dae(piecewise_control: bool = False) -> tuple[par_est.VariableList, par_est.Model]:
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
    variable_list.add_variable(par_est.VariableState("e0_T", 273.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i1", 3.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i2", 10.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i3", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i4", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_tot", 13.0))

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

    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_c_in_i1", 5.0, 4.0, 6.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i4", 0.0, 0.0, 1.0))
    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_T_in", 373.0, 353.0, 393.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))
    # fmt: on

    for var in variable_list.values():
        if isinstance(var, (par_est.VariableParameter, par_est.VariableControl)):
            var.guess = var.lower_bound

    if piecewise_control:
        var = variable_list["e0_T_in"].variable_list.index(0)
        var.guess = var.lower_bound
        var = variable_list["e0_c_in_i1"].variable_list.index(0)
        var.guess = var.lower_bound

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_T_in"].casadi_var - m.varlist_all["e0_T"].casadi_var))) + (((m.varlist_all["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (m.varlist_all["e0_c_p"].casadi_var * e0_V))) * ((m.varlist_all["e0_T_j"].casadi_var - m.varlist_all["e0_T"].casadi_var)))) + (((-m.varlist_all["e0_greek_Deltah_r1"].casadi_var) / (e0_greek_rho * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r2"].casadi_var) / (e0_greek_rho * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r3"].casadi_var) / (e0_greek_rho * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c1dot = ((((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i1"].casadi_var - m.varlist_all["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c2dot = ((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i2"].casadi_var - m.varlist_all["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c3dot = ((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i3"].casadi_var - m.varlist_all["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))
    c4dot = ((m.varlist_all["e0_F"].casadi_var / e0_V) * ((m.varlist_all["e0_c_in_i4"].casadi_var - m.varlist_all["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (e0_R * m.varlist_all["e0_T"].casadi_var))))))

    ctot = m.varlist_all["e0_c_tot"].casadi_var - m.varlist_all["e0_c_i1"].casadi_var - m.varlist_all["e0_c_i2"].casadi_var - m.varlist_all["e0_c_i3"].casadi_var - m.varlist_all["e0_c_i4"].casadi_var
    # fmt: on

    m.add_equations_differential([tdot, c1dot, c2dot, c3dot, c4dot])
    m.add_equations_algebraic([ctot])

    return variable_list, m


def cstr_ode_constant(piecewise_control: bool = False) -> tuple[par_est.VariableList, par_est.Model]:

    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.VariableState("e0_T", 273.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i1", 3.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i2", 10.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i3", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i4", 0.0))

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

    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_c_in_i1", 5.0, 4.0, 6.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i4", 0.0, 0.0, 1.0))
    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_T_in", 373.0, 353.0, 393.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))

    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i1_r1", -1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i1_r2", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i2_r2", -1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i3_r1", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i1_r3", -1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i4_r3", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_rho", 800.0))
    variable_list.add_variable(par_est.VariableConstant("e0_A", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_R", 8.314))
    variable_list.add_variable(par_est.VariableConstant("e0_V", 1.0))
    # fmt: on

    for var in variable_list.values():
        var.guess = var.lower_bound

    if piecewise_control:
        var = variable_list["e0_T_in"].variable_list.index(0)
        var.guess = var.lower_bound
        var = variable_list["e0_c_in_i1"].variable_list.index(0)
        var.guess = var.lower_bound

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_T_in"].casadi_var - m.varlist_all["e0_T"].casadi_var))) + (((m.varlist_all["e0_U"].casadi_var * m.varlist_all["e0_A"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * (m.varlist_all["e0_c_p"].casadi_var * m.varlist_all["e0_V"].casadi_var))) * ((m.varlist_all["e0_T_j"].casadi_var - m.varlist_all["e0_T"].casadi_var)))) + (((-m.varlist_all["e0_greek_Deltah_r1"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r2"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r3"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c1dot = ((((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i1"].casadi_var - m.varlist_all["e0_c_i1"].casadi_var))) + (m.varlist_all["e0_greek_nu_i1_r1"].casadi_var * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (m.varlist_all["e0_greek_nu_i1_r2"].casadi_var * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (m.varlist_all["e0_greek_nu_i1_r3"].casadi_var * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c2dot = ((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i2"].casadi_var - m.varlist_all["e0_c_i2"].casadi_var))) + (m.varlist_all["e0_greek_nu_i2_r2"].casadi_var * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c3dot = ((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i3"].casadi_var - m.varlist_all["e0_c_i3"].casadi_var))) + (m.varlist_all["e0_greek_nu_i3_r1"].casadi_var * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c4dot = ((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i4"].casadi_var - m.varlist_all["e0_c_i4"].casadi_var))) + (m.varlist_all["e0_greek_nu_i4_r3"].casadi_var * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    # fmt: on

    m.add_equations_differential([tdot, c1dot, c2dot, c3dot, c4dot])

    return variable_list, m


def cstr_dae_constant(piecewise_control: bool = False) -> tuple[par_est.VariableList, par_est.Model]:
    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.VariableState("e0_T", 273.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i1", 3.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i2", 10.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i3", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_c_i4", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_tot", 13.0))

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

    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_c_in_i1", 5.0, 4.0, 6.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_c_in_i1", 5.0, 4.0, 6.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i2", 10.0, 9.0, 11.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i3", 0.0, 0.0, 1.0))
    variable_list.add_variable(par_est.VariableControl("e0_c_in_i4", 0.0, 0.0, 1.0))
    if piecewise_control:
        variable_list.add_variable(par_est.VariableControlPiecewiseConstant("e0_T_in", 373.0, 353.0, 393.0))
    else:
        variable_list.add_variable(par_est.VariableControl("e0_T_in", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_T_j", 373.0, 353.0, 393.0))
    variable_list.add_variable(par_est.VariableControl("e0_F", 6.5e-4, 6.0e-4, 7.0e-4))

    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i1_r1", -1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i1_r2", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i2_r2", -1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i3_r1", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i1_r3", -1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_nu_i4_r3", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_greek_rho", 800.0))
    variable_list.add_variable(par_est.VariableConstant("e0_A", 1.0))
    variable_list.add_variable(par_est.VariableConstant("e0_R", 8.314))
    variable_list.add_variable(par_est.VariableConstant("e0_V", 1.0))
    # fmt: on

    for var in variable_list.values():
        if not var.name == "e0_c_tot":
            var.guess = var.lower_bound

    if piecewise_control:
        var = variable_list["e0_T_in"].variable_list.index(0)
        var.guess = var.lower_bound
        var = variable_list["e0_c_in_i1"].variable_list.index(0)
        var.guess = var.lower_bound

    m = par_est.Model(variable_list)

    # fmt: off
    tdot = (((((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_T_in"].casadi_var - m.varlist_all["e0_T"].casadi_var))) + (((m.varlist_all["e0_U"].casadi_var * m.varlist_all["e0_A"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * (m.varlist_all["e0_c_p"].casadi_var * m.varlist_all["e0_V"].casadi_var))) * ((m.varlist_all["e0_T_j"].casadi_var - m.varlist_all["e0_T"].casadi_var)))) + (((-m.varlist_all["e0_greek_Deltah_r1"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r2"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (((-m.varlist_all["e0_greek_Deltah_r3"].casadi_var) / (m.varlist_all["e0_greek_rho"].casadi_var * m.varlist_all["e0_c_p"].casadi_var)) * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c1dot = ((((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i1"].casadi_var - m.varlist_all["e0_c_i1"].casadi_var))) + (m.varlist_all["e0_greek_nu_i1_r1"].casadi_var * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (m.varlist_all["e0_greek_nu_i1_r2"].casadi_var * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))) + (m.varlist_all["e0_greek_nu_i1_r3"].casadi_var * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c2dot = ((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i2"].casadi_var - m.varlist_all["e0_c_i2"].casadi_var))) + (m.varlist_all["e0_greek_nu_i2_r2"].casadi_var * (m.varlist_all["e0_k_pre_r2"].casadi_var * (m.varlist_all["e0_c_i2"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r2"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c3dot = ((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i3"].casadi_var - m.varlist_all["e0_c_i3"].casadi_var))) + (m.varlist_all["e0_greek_nu_i3_r1"].casadi_var * (m.varlist_all["e0_k_pre_r1"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r1"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))
    c4dot = ((m.varlist_all["e0_F"].casadi_var / m.varlist_all["e0_V"].casadi_var) * ((m.varlist_all["e0_c_in_i4"].casadi_var - m.varlist_all["e0_c_i4"].casadi_var))) + (m.varlist_all["e0_greek_nu_i4_r3"].casadi_var * (m.varlist_all["e0_k_pre_r3"].casadi_var * (m.varlist_all["e0_c_i1"].casadi_var * ca.exp(((-m.varlist_all["e0_E_r3"].casadi_var) / (m.varlist_all["e0_R"].casadi_var * m.varlist_all["e0_T"].casadi_var))))))

    ctot = m.varlist_all["e0_c_tot"].casadi_var - m.varlist_all["e0_c_i1"].casadi_var - m.varlist_all["e0_c_i2"].casadi_var - m.varlist_all["e0_c_i3"].casadi_var - m.varlist_all["e0_c_i4"].casadi_var
    # fmt: on

    m.add_equations_differential([tdot, c1dot, c2dot, c3dot, c4dot])
    m.add_equations_algebraic([ctot])

    return variable_list, m


def vle_nle_problem() -> tuple[par_est.VariableList, par_est.Model]:
    # Id. VLE of EtOH and Water

    # Variables
    variable_list = par_est.variables.VariableList()  # Preallocate variable_list

    # Define variables
    #     T in K
    #     x in 1
    #     P in Pa
    #     # EtOH = 1,      H2O = 2
    #     a = [5.24125,    5.19625] # a in 1
    #     b = [1592.864,   1730.630]# b in K
    #     c = [-46.9659,   -39.7239] # c in K

    variable_list.add_variable(par_est.VariableAlgebraic("T", 373))
    variable_list.add_variable(par_est.VariableControl("x", 0.5))
    variable_list.add_variable(par_est.VariableControl("P", 1e5))
    variable_list.add_variable(par_est.VariableParameter("a1", 5.24125))
    variable_list.add_variable(par_est.VariableParameter("a2", 5.19625))
    variable_list.add_variable(par_est.VariableParameter("b1", 1592.864))
    variable_list.add_variable(par_est.VariableParameter("b2", 1730.630))
    variable_list.add_variable(par_est.VariableParameter("c1", -46.9659))
    variable_list.add_variable(par_est.VariableParameter("c2", -39.7239))

    model = par_est.Model(variable_list)  # adding all variables to the model

    # Equations
    RES = model.varlist_all["P"].casadi_var - (
        model.varlist_all["x"].casadi_var
        * 10
        ** (
            model.varlist_all["a1"].casadi_var
            - model.varlist_all["b1"].casadi_var
            / (model.varlist_all["c1"].casadi_var + model.varlist_all["T"].casadi_var)
        )
        * 1e5
        + (1 - model.varlist_all["x"].casadi_var)
        * 10
        ** (
            model.varlist_all["a2"].casadi_var
            - model.varlist_all["b2"].casadi_var
            / (model.varlist_all["c2"].casadi_var + model.varlist_all["T"].casadi_var)
        )
        * 1e5
    )
    model.add_equations_algebraic([RES])  # adding the equations to model

    return variable_list, model


def bod_model() -> tuple[par_est.VariableList, par_est.Model, list[par_est.VariableList]]:
    # BOD data as used in Bates, Watts, Nonlinear regression analysis: Its applications
    variable_list = par_est.variables.VariableList()  # Preallocate variable_list

    variable_list.add_variable(par_est.VariableAlgebraic("f", 8.3))
    variable_list.add_variable(par_est.VariableControl("x", 1))
    variable_list.add_variable(par_est.VariableParameter("theta1", 20))
    variable_list.add_variable(par_est.VariableParameter("theta2", 0.24))

    m = par_est.Model(variable_list)  # adding all variables to the model

    f = m.varlist_all["f"].casadi_var  # noqa: E501
    x = m.varlist_all["x"].casadi_var  # noqa: E501
    theta1 = m.varlist_all["theta1"].casadi_var  # noqa: E501
    theta2 = m.varlist_all["theta2"].casadi_var  # noqa: E501

    equation = f - (theta1 * (1 - ca.exp(-theta2 * x)))
    # Equations
    m.add_equations_algebraic([equation])  # adding the equations to model

    data = [
        [1, 8.3],
        [2, 10.3],
        [3, 19.0],
        [4, 16.0],
        [5, 15.6],
        [7, 19.8],
    ]

    exp_data = []

    for x_i, f_i in data:
        var_list = copy.deepcopy(variable_list)
        var_list["f"].dataframe.iloc[0] = f_i
        var_list["x"].dataframe.iloc[0] = x_i
        var_list["theta1"].fixed = False
        var_list["theta1"].lower_bound = 0
        var_list["theta1"].upper_bound = 40
        var_list["theta2"].fixed = False
        var_list["theta2"].lower_bound = 0
        var_list["theta2"].upper_bound = 1
        exp_data.append(var_list)

    return variable_list, m, exp_data


def isomerization_model() -> tuple[par_est.VariableList, par_est.Model, list[par_est.VariableList]]:
    # Isomerization data as used in Bates, Watts, Nonlinear regression analysis: Its applications
    variable_list = par_est.variables.VariableList()  # Preallocate variable_list

    variable_list.add_variable(par_est.VariableAlgebraic("f", 8.3))
    variable_list.add_variable(par_est.VariableControl("x1", 1))
    variable_list.add_variable(par_est.VariableControl("x2", 1))
    variable_list.add_variable(par_est.VariableControl("x3", 1))
    variable_list.add_variable(par_est.VariableParameter("theta1", 35.92))
    variable_list.add_variable(par_est.VariableParameter("theta2", 0.0708))
    variable_list.add_variable(par_est.VariableParameter("theta3", 0.0377))
    variable_list.add_variable(par_est.VariableParameter("theta4", 0.167))

    m = par_est.Model(variable_list)  # adding all variables to the model

    f = m.varlist_all["f"].casadi_var  # noqa: E501
    x1 = m.varlist_all["x1"].casadi_var  # noqa: E501
    x2 = m.varlist_all["x2"].casadi_var  # noqa: E501
    x3 = m.varlist_all["x3"].casadi_var  # noqa: E501
    theta1 = m.varlist_all["theta1"].casadi_var  # noqa: E501
    theta2 = m.varlist_all["theta2"].casadi_var  # noqa: E501
    theta3 = m.varlist_all["theta3"].casadi_var  # noqa: E501
    theta4 = m.varlist_all["theta4"].casadi_var  # noqa: E501

    equation = f - (theta1 * theta3 * (x2 - x3 / 1.632)) / (
        1 + theta2 * x1 + theta3 * x2 + theta4 * x3
    )
    # Equations
    m.add_equations_algebraic([equation])  # adding the equations to model

    data = [
        [205.8, 90.9, 37.1, 3.541],
        [404.8, 92.9, 36.3, 2.397],
        [209.7, 174.9, 49.4, 6.694],
        [401.6, 187.2, 44.9, 4.722],
        [224.9, 92.7, 116.3, 0.593],
        [402.6, 102.2, 128.9, 0.268],
        [212.7, 186.9, 134.4, 2.797],
        [406.2, 192.6, 134.9, 2.451],
        [133.3, 140.8, 87.6, 3.196],
        [470.9, 144.2, 86.9, 2.021],
        [300.0, 68.3, 81.7, 0.896],
        [301.6, 214.6, 101.7, 5.084],
        [297.3, 142.2, 10.5, 5.686],
        [314.0, 146.7, 157.1, 1.193],
        [305.7, 142.0, 86.0, 2.648],
        [300.1, 143.7, 90.2, 3.303],
        [305.4, 141.1, 87.4, 3.054],
        [305.2, 141.5, 87.0, 3.302],
        [300.1, 83.0, 66.4, 1.271],
        [106.6, 209.6, 33.0, 11.648],
        [417.2, 83.9, 32.9, 2.002],
        [251.0, 294.4, 41.5, 9.604],
        [250.3, 148.0, 14.7, 7.754],
        [145.1, 291.0, 50.2, 11.590],
    ]

    exp_data = []

    for x1_i, x2_i, x3_i, f_i in data:
        var_list = copy.deepcopy(variable_list)
        var_list["f"].dataframe.iloc[0] = f_i
        var_list["x1"].dataframe.iloc[0] = x1_i
        var_list["x2"].dataframe.iloc[0] = x2_i
        var_list["x3"].dataframe.iloc[0] = x3_i
        var_list["theta1"].fixed = False
        var_list["theta2"].fixed = False
        var_list["theta3"].fixed = False
        var_list["theta4"].fixed = False
        var_list.set_bounds()
        exp_data.append(var_list)

    return variable_list, m, exp_data
