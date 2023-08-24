import logging
import casadi as ca
import numpy as np

import copy
import par_est.examples
import par_est
import par_est.tools
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
def test_pe_intials_algebraic(piecewise):
    variable_list1, model = par_est.examples.empy_dae(piecewise)
    variable_list1["X1"].set_dataframe_from_value_and_time([1, 0], [0, 1])
    variable_list1["X2"].set_dataframe_from_value_and_time([0], [0])
    variable_list2 = copy.deepcopy(variable_list1)
    variable_list2["X1"].set_dataframe_from_value_and_time([1, 0], [0, 1])
    variable_list2["X2"].set_dataframe_from_value_and_time([1, 0], [0, 1])
    pe = par_est.ParameterEstimation(model, [variable_list1, variable_list2])
    pe.solver_settings = {
        "ipopt": {
            "max_iter": 0,
        },
    }
    assert pe.list_simulators[0]._initial_algebraic[0] == -1
    assert pe.list_simulators[1]._initial_algebraic[0] == -2


@pytest.mark.parametrize("piecewise", [True, False])
def test_pe_objective(piecewise):
    variable_list1, model = par_est.examples.empy_dae(piecewise)

    for var in variable_list1.values():
        if isinstance(
            var,
            (
                par_est.VariableControl,
                par_est.VariableParameter,
                par_est.VariableControlPiecewiseConstant,
            ),
        ):
            var.fixed = False

    variable_list1["X1"].set_dataframe_from_value_and_time([0, 1], [0, 1])
    variable_list1["X2"].set_dataframe_from_value_and_time([0], [0])
    variable_list1["X2"].variance = 10
    variable_list2 = copy.deepcopy(variable_list1)
    variable_list2["X1"].set_dataframe_from_value_and_time([0, 3], [0, 1])
    variable_list2["X1"].variance = 20
    variable_list2["X2"].set_dataframe_from_value_and_time([0, 4], [0, 2])
    variable_list2["X2"].variance = 40
    pe = par_est.ParameterEstimation(model, [variable_list1, variable_list2])
    pe.solver_settings = {
        "ipopt": {
            "max_iter": 0,
        },
    }

    weight = np.array([[1.5, 1.5], [1. , 1.], [1. , 1.]])
    var = np.array([[1., 0.1], [0.05 , 0.025], [0.05 , 0.025]])
    data = np.array([[1., 0], [3., 0], [0, 4.]])
    mask = np.array([[1, 0], [1, 0], [0, 1]])

    assert_numpy = np.testing.assert_array_equal
    assert_numpy(data, pe.array_data)
    assert_numpy(var, pe.array_inverted_variance)
    assert_numpy(weight, pe.experiments_weights)
    assert_numpy(mask, pe.array_data_mask)

    obj = np.sum(var * (data * mask) ** 2)
    obj_weight = np.sum(weight * var * (data * mask) ** 2)
    for switch in [True, False]:
        res = pe.optimize(switch)
        res_weight = pe.optimize(switch, scale_experiments=True)
        assert res["f"] == obj
        assert np.isclose(res_weight["f"], obj_weight)


@pytest.mark.parametrize("piecewise", [True, False])
def test_pe(piecewise):
    """Test that ParameterEstimation on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """

    """ ODE and DAE """
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        var_list, model = cstr_model(piecewise)
        time_grid = np.linspace(10, 10000, 3)
        time_grid = np.insert(time_grid, 0, 0)

        if piecewise:
            T_in = var_list["e0_T_in"]
            T_in.expand_horizon([2000, 4000], [373, 373])

        var_list_fixed = copy.deepcopy(var_list)
        for var in var_list_fixed.values():
            var.fixed = True
        var_list_exp = par_est.Simulator(
            model, time_grid, var_list_fixed
        ).generate_exp_data()

        for key, var in var_list_exp.items():
            var_list[key] = var

        for var in var_list.values():
            var.fixed = True

        var_list["e0_E_r1"].fixed = False

        var_list["e0_T"].dataframe = var_list["e0_T"].dataframe.iloc[:2]
        var_list["e0_c_i1"].lower_bound = 0
        var_list["e0_c_i1"].upper_bound = None

        pe = par_est.ParameterEstimation(
            model, [var_list]
        )

        if model.DAE:
            answer_scaled = 1.26485e-13
            answer = 1.02115e-11
        else:
            answer_scaled = 2.04875e-18
            answer = 6.48183e-12

        res = pe.optimize()
        logging.warning(
            f"Model.DAE: {model.DAE}, Result: {res['f']}, Expecting: {answer_scaled}"
        )
        assert np.isclose(res["f"], ca.DM(answer_scaled), rtol=0, atol=1.0e-9)

        res = pe.optimize(False)
        logging.warning(
            f"Model.DAE: {model.DAE}, Result: {res['f']}, Expecting: {answer}"
        )
        assert np.isclose(res["f"], ca.DM(answer), rtol=0, atol=1.0e-9)


if __name__ == "__main__":
    pass
    test_pe(True,True)
    test_pe_objective(False)
    test_pe_intials_algebraic()
