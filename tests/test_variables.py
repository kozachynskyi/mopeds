import mopeds
import pandas as pd
from mopeds.variables import SameVariableNameError

import casadi as ca
import pytest


def test_options():
    mopeds.set_options(variable_scaling=True)
    assert mopeds.get_options()["variable_scaling"]
    with mopeds.options(variable_scaling=False):
        assert mopeds.get_options()["variable_scaling"] is False
    assert mopeds.get_options()["variable_scaling"]
    mopeds.set_options(variable_scaling=False)
    assert mopeds.get_options()["variable_scaling"] is False
    mopeds.set_options(variable_scaling=True)
    assert mopeds.get_options()["variable_scaling"]


def test_variables():
    # Test VariableList() and subclasses
    variable_list = mopeds.VariableList()
    counter = 0

    for variable_type in mopeds.Variable.get_subclasses():
        var = variable_type(f"Name{counter}")
        variable_list.add_variable(var)
        counter += 1

    assert len(variable_list) == sum(1 for _ in mopeds.Variable.get_subclasses())

    # Test SameVariableNameError
    variable_list = mopeds.VariableList()

    var = mopeds.VariableState("Name")
    variable_list.add_variable(var)
    assert var.casadi_var.name() == "Name"
    with pytest.raises(SameVariableNameError):
        variable_list.add_variable(var)

    assert len(variable_list) == 1

    # Test PlottingError
    variable_list = mopeds.VariableList()

    var_1 = mopeds.VariableState("Var1")
    variable_list.add_variable(var_1)

    with pytest.raises(ValueError):
        variable_list["Var1"].set_dataframe_from_value_and_time([0, 1], [0])
    with pytest.raises(ValueError):
        variable_list["Var1"].set_dataframe_from_value_and_time([0], [0, 1])
    with pytest.raises(ValueError):
        variable_list["Var1"].set_dataframe_from_value_and_time([0, 1], [1, 2])

    variable_list["Var1"].set_dataframe_from_value_and_time([0, 1], [0, 2])
    a = variable_list["Var1"]

    var_1 = mopeds.VariableControlPiecewiseConstant("Var1", 20)
    var_2 = mopeds.VariableControlPiecewiseConstant("Var1", 20)

    var_2.variable_list.index(0).dataframe.rename(
        index={
            var_2.variable_list.index(0).origin_ts: pd.Timestamp(
                year=2021, month=1, day=1
            )
        },
        inplace=True,
    )

    var_3 = mopeds.VariableControlPiecewiseConstant("Var1", None)

    var_1 = mopeds.VariableControlPiecewiseConstant("Var1", 20)

    assert list(variable_list["Var1"].time_absolute) == [
        pd.Timestamp("1970-01-01 00:00:00"),
        pd.Timestamp("1970-01-01 00:00:02"),
    ]
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
    assert list(var_1.time_absolute) == [
        pd.Timestamp("1970-01-01 00:00:00"),
        pd.Timestamp("1970-01-01 00:00:11"),
    ]
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
    assert list(var_2.time_absolute) == [
        pd.Timestamp("2021-01-01 00:00:00"),
        pd.Timestamp("2021-01-01 00:00:11"),
        pd.Timestamp("2021-01-01 00:00:11.300"),
    ]
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
    var = mopeds.VariableAlgebraic("var")
    assert var.get_constraint_idas == 0
    var = mopeds.VariableAlgebraic("var", None, -1, None)
    assert var.get_constraint_idas == 0
    var = mopeds.VariableAlgebraic("var", None, 0, None)
    assert var.get_constraint_idas == 1
    var = mopeds.VariableAlgebraic("var", None, None, 1)
    assert var.get_constraint_idas == 0
    var = mopeds.VariableAlgebraic("var", None, None, -1)
    assert var.get_constraint_idas == -2
    var = mopeds.VariableAlgebraic("var", None, None, 0)
    assert var.get_constraint_idas == -1


def test_varlist():
    var_1 = mopeds.VariableControlPiecewiseConstant("v1", 20)
    var_2 = mopeds.VariableControlPiecewiseConstant("v2", 20)
    var_3 = mopeds.VariableControlPiecewiseConstant("v3", 20)
    vl = mopeds.VariableList()
    vl_unordered = mopeds.VariableList()

    vl.add_variable(var_1)
    vl.add_variable(var_2)

    vl_unordered.add_variable(var_2)
    vl_unordered.add_variable(var_1)
    vl_unordered.add_variable(var_3)

    assert list(vl.keys()) == ["v1", "v2"]
    with pytest.raises(KeyError):
        vl._get_sorted_varlist(vl_unordered).keys()
    assert list(vl._get_sorted_varlist(vl_unordered, raise_error=False).keys()) == [
        "v2",
        "v1",
    ]


if __name__ == "__main__":
    test_options()
    # test_variables()
    test_varlist()
