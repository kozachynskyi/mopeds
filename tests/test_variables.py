import par_est
from par_est.variables import SameVariableNameError, PlottingError

import casadi as ca
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

    var_1 = par_est.VariableControlPiecewiseConstant("Var1", 20)
    var_2 = par_est.VariableControlPiecewiseConstant("Var1", 20)
    var_3 = par_est.VariableControlPiecewiseConstant("Var1", None)

    var_1 = par_est.VariableControlPiecewiseConstant("Var1", 20)
    assert len(var_1.time) == 1
    assert var_1.time == [0]
    assert var_1.var_at_time(0).value == 20
    assert var_1.var_at_time(0).fixed is True
    var_1.expand_horizon([11], [4])
    assert len(var_1.time) == 2
    assert var_1.time == [0, 11]
    assert var_1.var_at_time(10).value == 20
    assert var_1.var_at_time(11).value == 4
    assert var_1.var_at_time(12).value == 4
    with pytest.raises(NotImplementedError):
        var_1.expand_horizon([11], [4])
    with pytest.raises(NotImplementedError):
        var_1.set_horizon([11], [4])
    with pytest.raises(ValueError):
        var_1.expand_horizon([11, 12], [4])

    var_2.expand_horizon([11, 11.3], [4, 5])
    assert len(var_2.time) == 3
    assert var_2.time == [0, 11, 11.3]
    assert var_2.fixed is True
    assert var_2.var_at_time(10).value == 20
    assert var_2.var_at_time(11).value == 4
    assert var_2.var_at_time(11.299999).value == 4
    assert var_2.var_at_time(11.300001).value == 5
    assert var_2.var_at_time(12).value == 5

    assert var_3.var_at_time(0).fixed is True
    var_3.expand_horizon([11, 11.3], [4, None])
    assert len(var_3.time) == 3
    assert var_3.fixed is True
    assert var_3.var_at_time(10).get_value_based_on_fixed() is None
    assert var_3.var_at_time(11).get_value_based_on_fixed() == 4
    var_3.var_at_time(12).fixed = False
    assert isinstance(var_3.var_at_time(12).get_value_based_on_fixed(), ca.MX)

    var_3.fixed = True
    for var in var_3.variable_list.values():
        assert var.fixed is True

    var_3.fixed = False
    for var in var_3.to_dictionary().values():
        assert var.fixed is False
    assert list(var_3.to_dictionary().keys()) == [0, 11, 11.3]

    # Test constraints_idas prperty
    var = par_est.VariableAlgebraic("var")
    assert var.constraint_idas == 0
    var = par_est.VariableAlgebraic("var", None, -1, None)
    assert var.constraint_idas == 0
    var = par_est.VariableAlgebraic("var", None, 0, None)
    assert var.constraint_idas == 1
    var = par_est.VariableAlgebraic("var", None, None, 1)
    assert var.constraint_idas == 0
    var = par_est.VariableAlgebraic("var", None, None, -1)
    assert var.constraint_idas == -2
    var = par_est.VariableAlgebraic("var", None, None, 0)
    assert var.constraint_idas == -1


if __name__ == "__main__":
    test_variables()
