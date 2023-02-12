import par_est
import pandas as pd
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
    variable_list.add_variable(var_1)

    with pytest.raises(ValueError):
        variable_list["Var1"].set_dataframe_from_value_and_time([0, 1], [0])
    with pytest.raises(ValueError):
        variable_list["Var1"].set_dataframe_from_value_and_time([0], [0, 1])
    with pytest.raises(ValueError):
        variable_list["Var1"].set_dataframe_from_value_and_time([0, 1], [1, 2])

    variable_list["Var1"].set_dataframe_from_value_and_time([0, 1], [0, 2])
    a = variable_list["Var1"]

    var_1 = par_est.VariableControlPiecewiseConstant("Var1", 20)
    var_2 = par_est.VariableControlPiecewiseConstant("Var1", 20)

    var_2.variable_list.index(0).dataframe.rename(
        index={
            var_2.variable_list.index(0).origin_ts: pd.Timestamp(
                year=2021, month=1, day=1
            )
        },
        inplace=True,
    )

    var_3 = par_est.VariableControlPiecewiseConstant("Var1", None)

    var_1 = par_est.VariableControlPiecewiseConstant("Var1", 20)

    assert variable_list["Var1"].time_absolute.equals(
        pd.DatetimeIndex(
            ["1970-01-01 00:00:00", "1970-01-01 00:00:02"],
            dtype="datetime64[ns]",
            freq=None,
        )
    )
    assert variable_list["Var1"].time_relative == [0, 2]
    assert len(var_1.time_relative) == 1
    assert len(var_1.time_absolute) == 1
    assert var_1.time_relative == [0]
    var_1.get_variable_at_time_absolute("1970-01-01 00:00:00")
    assert var_1.get_variable_at_time_relative(
        0
    ) == var_1.get_variable_at_time_absolute("1970-01-01 00:00:00")
    assert var_1.get_variable_at_time_relative(0).fixed is True
    var_1.expand_horizon([11], [4])
    assert var_1.time_absolute.equals(
        pd.Series(
            ["1970-01-01 00:00:00", "1970-01-01 00:00:11"], dtype="datetime64[ns]"
        )
    )
    assert var_1.time_relative == [0, 11]
    assert len(var_1.time_relative) == 2
    assert len(var_1.time_absolute) == 2
    assert var_1.get_variable_at_time_relative(10).value == [20]
    assert var_1.get_variable_at_time_relative(11).value == [4]
    assert var_1.get_variable_at_time_relative(12).value == [4]
    with pytest.raises(NotImplementedError):
        var_1.expand_horizon([11], [4])
    with pytest.raises(NotImplementedError):
        var_1.set_horizon([11], [4])
    with pytest.raises(ValueError):
        var_1.expand_horizon([11, 12], [4])

    var_2.expand_horizon([11, 11.3], [4, 5])
    assert len(var_2.time_absolute) == 3
    assert len(var_2.time_relative) == 3
    assert var_2.time_relative == [0, 11, 11.3]
    assert var_2.time_absolute.equals(
        pd.Series(
            ["2021-01-01 00:00:00", "2021-01-01 00:00:11", "2021-01-01 00:00:11.300"],
            dtype="datetime64[ns]",
        )
    )
    assert var_2.fixed is True
    assert var_2.get_variable_at_time_relative(10).value == [20]
    var_2.value = 19
    assert var_2.get_variable_at_time_relative(10).value == [19]
    assert var_2.get_variable_at_time_relative(11).value == [4]
    assert var_2.get_variable_at_time_relative(11.299999).value == [4]
    assert var_2.get_variable_at_time_relative(11.300001).value == [5]
    assert var_2.get_variable_at_time_relative(12).value == [5]

    assert var_3.get_variable_at_time_relative(0).fixed is True
    var_3.expand_horizon([11, 11.3], [4, None])
    assert len(var_3.time_absolute) == 3
    assert len(var_3.time_relative) == 3
    assert var_3.fixed is True
    assert pd.isna(var_3.get_variable_at_time_relative(10).get_value_or_casadi())
    assert var_3.get_variable_at_time_relative(11).get_value_or_casadi() == 4
    var_3.get_variable_at_time_relative(12).fixed = False
    assert isinstance(
        var_3.get_variable_at_time_relative(12).get_value_or_casadi(), ca.MX
    )

    var_3.fixed = True
    for var in var_3.variable_list.values():
        assert var.fixed is True

    var_3.fixed = False
    for var in var_3.to_dictionary().values():
        assert var.fixed is False
    assert list(var_3.to_dictionary().keys()) == [0, 11, 11.3]

    # Test constraints_idas prperty
    var = par_est.VariableAlgebraic("var")
    assert var.get_constraint_idas == 0
    var = par_est.VariableAlgebraic("var", None, -1, None)
    assert var.get_constraint_idas == 0
    var = par_est.VariableAlgebraic("var", None, 0, None)
    assert var.get_constraint_idas == 1
    var = par_est.VariableAlgebraic("var", None, None, 1)
    assert var.get_constraint_idas == 0
    var = par_est.VariableAlgebraic("var", None, None, -1)
    assert var.get_constraint_idas == -2
    var = par_est.VariableAlgebraic("var", None, None, 0)
    assert var.get_constraint_idas == -1


if __name__ == "__main__":
    test_variables()
