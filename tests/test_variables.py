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

    # Test PlottingError
    variable_list = par_est.VariableList()

    var_1 = par_est.VariableState("Var1")
    var_2 = par_est.VariableState("Var2")

    variable_list.add_variable(var_1)
    variable_list.add_variable(var_2)

    variable_list["Var1"].value.value = [0, 1, 2]
    variable_list["Var2"].value.time = [0, 1, 2]

    with pytest.raises(PlottingError):
        variable_list.plot_states()

    with pytest.raises(PlottingError):
        variable_list["Var1"].value.time = [0, 1, 2]
        variable_list.plot_states()


if __name__ == "__main__":
    test_variables()
