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

    obj = np.sum(var * (data * mask) ** 2) / 2
    obj_weight = np.sum(weight * var * (data * mask) ** 2) / 2
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


@pytest.mark.parametrize("piecewise", [True, False])
def test_oed(piecewise):
    """Test that OptimalExperimentalDesign on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        var_list, model = cstr_model(piecewise)
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)
        for var in var_list.values():
            var.fixed = True

        var_list["e0_E_r1"].fixed = False
        var_list["e0_T_in"].fixed = True
        var_list["e0_c_in_i1"].fixed = False

        if piecewise:
            c_in_i1 = var_list["e0_c_in_i1"]
            c_in_i1.expand_horizon([3340, 6670], [5, 5])
            c_in_i1.fixed = True
            c_in_i1.variable_list.index(0).fixed = False

        oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
        res = oed.optimize()

        logging.warning(f"{res['f']}")
        assert np.isclose(res["f"], ca.DM(45.1675), rtol=0, atol=1.0e-4)

        # For not functionality is turnded off
        # res = oed.optimize(True)

        # logging.warning(f"{res['f']}")
        # assert np.isclose(res["f"], ca.DM(45.1675), rtol=0, atol=1.0e-4)


def test_oed_piecewise():
    """Test if OED works consistently with different horizons"""
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        var_list, model = cstr_model(True)
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)
        for var in var_list.values():
            var.fixed = True

        var_list["e0_E_r1"].fixed = False
        var_list["e0_T_in"].fixed = True
        var_list["e0_c_in_i1"].fixed = False

        c_in_i1 = var_list["e0_T_in"]
        c_in_i1.expand_horizon([10, 3340, 6670], [373, 353, 393])
        c_in_i1.fixed = True

        if model.DAE:
            res_values = [
                31.16299587,
                7.9313022,
                10.47423272,
                7.91570605,
                10.69082675,
                9.69286661,
            ]
        else:
            res_values = [
                31.162953,
                7.93128491,
                10.47422549,
                7.91568793,
                10.69076264,
                9.69266026,
            ]

        list_orders = [
            [False, True, True],
            [False, False, True],
            [False, False, False],
            [True, False, True],
            [True, False, False],
            [True, True, False],
        ]
        for fix_order, result_expected in zip(list_orders, res_values):
            for fix, var in zip(fix_order, c_in_i1.variable_list.values()):
                var.fixed = fix

            oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
            # oed.solver_settings["ipopt"]["max_iter"] = 10
            oed.solver_settings = {
                "ipopt": {
                    "hessian_approximation": "limited-memory",
                    "max_iter": 5,
                    "print_level": 0,
                },
            }
            res = oed.optimize()

            logging.warning(f"{res['f'].full()}")
            logging.warning(f"{fix_order}")
            logging.warning(f"{cstr_model}")
            assert np.isclose(res["f"], ca.DM(result_expected))


@pytest.mark.parametrize("piecewise", [True, False])
def test_optimizer(piecewise):  # noqa: C901
    """Tests if optimizer can deal with variable list of fixed and unfixed parameters.
    Not well designed, and may yield false positives, but let it be.
    """
    variable_list, m = par_est.examples.cstr_ode(piecewise)

    for var in variable_list.values():
        if isinstance(
            var,
            (
                par_est.VariableControl,
                par_est.VariableParameter,
                par_est.VariableControlPiecewiseConstant,
            ),
        ):
            var.fixed = False

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 1000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    # Generate experimental data for parameter estimation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).generate_exp_data()

    # Replace empty state variables with results from simulation
    for key, var in var_list_exp.items():
        variable_list[key] = var

    for i in range(5):
        if i == 1:
            variable_list["e0_U"].fixed = True
        elif i == 2:
            variable_list["e0_T_in"].fixed = True
        elif i == 3:
            variable_list["e0_c_i2"].fixed = True
        elif i == 4:
            for var in variable_list.values():
                var.fixed = True
            variable_list["e0_U"].fixed = False
            variable_list["e0_T_in"].fixed = False

        for j in range(2):
            if j == 0:
                pe = par_est.ParameterEstimation(m, [variable_list])
            else:
                pe = par_est.ParameterEstimation(m, [variable_list, variable_list])

            pe.solver_settings = {
                "verbose": False,
                "ipopt": {
                    "hessian_approximation": "limited-memory",
                    "max_iter": 1,
                    "print_level": 0,
                },
            }

            oed = par_est.OptimalExperimentalDesign(m, [variable_list], time_grid)
            oed.solver_settings = {
                "verbose": False,
                "ipopt": {
                    "hessian_approximation": "limited-memory",
                    "max_iter": 1,
                    "print_level": 0,
                },
            }

            for ij in range(2):
                if ij == 0:
                    res_pe = pe.optimize()
                    if j == 0:
                        res_oed = oed.optimize()
                else:
                    res_pe = pe.optimize(False)
                    if j == 0:
                        res_oed = oed.optimize(False)

            if i == 0:
                assert res_pe["x"].size() == (11, 1)
                assert res_oed["x"].size() == (7, 1)
            elif i == 1:
                assert res_pe["x"].size() == (10, 1)
                assert res_oed["x"].size() == (7, 1)
            elif i == 2:
                assert res_pe["x"].size() == (10, 1)
                assert res_oed["x"].size() == (6, 1)
            elif i == 3:
                assert res_pe["x"].size() == (10, 1)
                assert res_oed["x"].size() == (6, 1)
            elif i == 4:
                assert res_pe["x"].size() == (1, 1)
                assert res_oed["x"].size() == (1, 1)


if __name__ == "__main__":
    pass
    # test_optimizer(True)
    # test_oed(True)
    # test_pe(True,True)
    test_pe_objective(False)
    # test_pe_intials_algebraic()
    # test_oed_piecewise()
