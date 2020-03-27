import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est


def initialize_problem_dae3():
    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.State_variable("x", 3.0))
    variable_list.add_variable(par_est.State_variable("u", -1.0 / 3))
    variable_list.add_variable(par_est.State_variable("y", 4.0))
    variable_list.add_variable(par_est.State_variable("v", 1.0 / 4))

    variable_list.add_variable(par_est.Algebraic_variable("lambda", 1147.0 / 720))

    variable_list.add_variable(par_est.Parameter_variable("L", 5.0))
    variable_list.add_variable(par_est.Parameter_variable("g", 10.0))
    # fmt: on

    m = par_est.Model(variable_list)

    # fmt: off
    dydx1 = m._all_variables["u"].casadi_var
    dydx2 = m._all_variables["lambda"].casadi_var * m._all_variables["x"].casadi_var
    dydx3 = m._all_variables["v"].casadi_var
    dydx4 = m._all_variables["lambda"].casadi_var * m._all_variables["y"].casadi_var + m._all_variables["g"].casadi_var

    alg1 = m._all_variables["x"].casadi_var ** 2 + m._all_variables["y"].casadi_var ** 2 - m._all_variables["L"].casadi_var ** 2

    m.add_differential_equations([dydx1, dydx2, dydx3, dydx4, ])
    m.add_algebraic_equations([alg1, ])
    # fmt: on

    return variable_list, m


def initialize_problem_dae1():
    variable_list = par_est.VariableList()

    # fmt: off
    variable_list.add_variable(par_est.State_variable("x", 3.0))
    variable_list.add_variable(par_est.State_variable("u", -1.0 / 3))

    variable_list.add_variable(par_est.Algebraic_variable("y", 4.0))
    variable_list.add_variable(par_est.Algebraic_variable("v", 1.0 / 4))
    variable_list.add_variable(par_est.Algebraic_variable("lambda", 1147.0 / 720))

    variable_list.add_variable(par_est.Parameter_variable("L", 5.0))
    variable_list.add_variable(par_est.Parameter_variable("g", 10.0))
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


if __name__ == "__main__":

    variable_list_dae1, m1 = initialize_problem_dae1()
    variable_list_dae3, m3 = initialize_problem_dae3()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(0.01, 1, 40)
    time_grid = np.insert(time_grid, 0, 0)

    var_list_fixed = copy.deepcopy(variable_list_dae1)
    for var in var_list_fixed.values():
        var.fixed = True

    sim_fixed = par_est.Simulator(m1, time_grid, var_list_fixed)
    res = sim_fixed.generate_exp_data()
    var_list_exp = sim_fixed.generate_exp_data()
    # var_list_exp.plot_states()
    r = sim_fixed.analyze()
    print(r)

    # # Generate experimental data
    # for [variable_list, m, sim_name] in zip(
    #     [variable_list_dae1, variable_list_dae3],
    #     [m1, m3],
    #     ["dae-index=1", "dae-index=3"],
    # ):
    #     var_list_fixed = copy.deepcopy(variable_list)
    #     for var in var_list_fixed.values():
    #         var.fixed = True

    #     try:
    #         sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
    #         res = sim_fixed.generate_exp_data()
    #         var_list_exp = sim_fixed.generate_exp_data()
    #         var_list_exp.plot_states()
    #     except Exception:
    #         print(f"Simulation failed for {sim_name}")
