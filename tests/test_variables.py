import par_est
from par_est.examples import cstr_ode
from par_est.variables import SameVariableNameError, PlottingError

import copy
import matplotlib.pyplot as plt
import numpy as np
import pytest


def test_variables():
    # Test VariableList() and subclasses
    variable_list = par_est.VariableList()
    counter = 0

    for variable_type in par_est.Variable.get_subclasses():
        var = variable_type(f"Name{counter}")
        variable_list.add_variable(var)
        counter += 1

    assert len(variable_list) == sum(1 for _ in par_est.Variable.get_subclasses())

    # Test SameVariableNameError
    variable_list = par_est.VariableList()

    var = par_est.VariableState("Name")
    variable_list.add_variable(var)
    assert var.casadi_var.name() == "Name"
    with pytest.raises(SameVariableNameError):
        variable_list.add_variable(var)

    assert len(variable_list) == 1

    # Test plot_states() and PlottingError
    variable_list, m = par_est.examples.cstr_ode()
    time_grid = np.linspace(0, 1, 2)
    time_grid = np.insert(time_grid, 0, 0)

    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
    res = sim_fixed.generate_exp_data()

    # PlottingError
    res_no_time = copy.deepcopy(res)
    res_no_var = copy.deepcopy(res)

    for var in res_no_time.values():
        var.value.time = np.array([])

    with pytest.raises(PlottingError):
        res_no_time.plot_states()

    for var in res_no_var.values():
        var.value.value = np.array([])

    with pytest.raises(PlottingError):
        res_no_var.plot_states()


if __name__ == "__main__":
    test_variables()
