import logging
import casadi as ca
import numpy as np
import sys
import warnings

import copy
import mopeds.examples
import mopeds
import mopeds.tools
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
def test_scaling(piecewise, dae):
    sens_fim = []
    time_grid = np.linspace(0, 1000, 4)
    var_list, model = mopeds.examples.cstr(piecewise, dae)
    if piecewise:
        T_in = var_list["e0_T_in"]
        T_in.expand_horizon([2000, 2500, 4000], [363, 383, 393])

    varlist_i = copy.deepcopy(var_list)

    rng = np.random.default_rng(0)

    with mopeds.options(variable_scaling=False):
        var_list_exp = mopeds.Simulator(model, time_grid, var_list).simulate()[2]
        for var in var_list_exp.get_state().values():
            var.value = rng.normal(var.value[0], var.variance**0.5)


    for scaling in [True, False]:
        with mopeds.options(variable_scaling=scaling):
            for key, var in var_list_exp.items():
                varlist_i[key] = var

            varlist_i["e0_U"].fixed = False
            varlist_i["e0_E_r1"].fixed = False

            opts = {
                "expand": 1,
                "abstol": 1e-10,
                "reltol": 1e-8,
            }

            pe = mopeds.ParameterEstimation(model, [varlist_i], simulator_settings=opts)
            sens_fim.append(pe.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_E_r1": 9.6e4}))
    
    for key in sens_fim[0].keys():
        # print(key)
        v1 = sens_fim[0][key] 
        v2 = sens_fim[1][key] 
        if key not in ["jac_sorted", "jac_scaled_sorted", "jac_yao_sorted", "jac_wls", "hess_wls", "hess_tikh"]:
            # print(v2/v1)
            assert np.isclose(v1, v2, rtol=1e-4).all()


@pytest.mark.parametrize("piecewise", [True, False])
def test_pe_intials_algebraic(piecewise):
    variable_list1, model = mopeds.examples.empy_dae(piecewise)
    variable_list1["X1"].set_dataframe_from_value_and_time([1, 0], [0, 1])
    variable_list1["X2"].set_dataframe_from_value_and_time([0], [0])
    variable_list2 = copy.deepcopy(variable_list1)
    variable_list2["X1"].set_dataframe_from_value_and_time([1, 0], [0, 1])
    variable_list2["X2"].set_dataframe_from_value_and_time([1, 0], [0, 1])
    pe = mopeds.ParameterEstimation(model, [variable_list1, variable_list2])
    pe.solver_settings = {
        "ipopt": {
            "max_iter": 0,
        },
    }
    assert pe.list_simulators[0]._initial_algebraic[0] == -1
    assert pe.list_simulators[1]._initial_algebraic[0] == -2


@pytest.mark.parametrize("piecewise", [True, False])
def test_pe_objective(piecewise):
    variable_list1, model = mopeds.examples.empy_dae(piecewise)

    for var in variable_list1.values():
        if isinstance(
            var,
            (
                mopeds.VariableControl,
                mopeds.VariableParameter,
                mopeds.VariableControlPiecewiseConstant,
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
    pe = mopeds.ParameterEstimation(model, [variable_list1, variable_list2])
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
    assert_numpy(weight, pe.experiments_weights)
    assert_numpy(mask, pe.array_data_mask)

    obj = np.sum(var * (data * mask) ** 2)
    obj_weight = np.sum(weight * var * (data * mask) ** 2)
    res = pe.optimize()
    res_weight = pe.optimize(scale_experiments=True)
    assert res["f"] == obj
    assert np.isclose(res_weight["f"], obj_weight)

@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_pe(piecewise, dae, use_constant, scaling):
    """Test that ParameterEstimation on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """

    """ ODE and DAE """
    with mopeds.options(variable_scaling=scaling):
        var_list, model = mopeds.examples.cstr(piecewise, dae, use_constant)
        time_grid = np.linspace(10, 10000, 3)
        time_grid = np.insert(time_grid, 0, 0)

        if piecewise:
            T_in = var_list["e0_T_in"]
            T_in.expand_horizon([2000, 4000], [373, 373])

        var_list_fixed = copy.deepcopy(var_list)
        for var in var_list_fixed.values():
            var.fixed = True
        var_list_exp = mopeds.Simulator(
            model, time_grid, var_list_fixed
        ).simulate()[2]

        for key, var in var_list_exp.items():
            var_list[key] = var

        for var in var_list.values():
            var.fixed = True

        var_list["e0_E_r1"].fixed = False

        var_list["e0_T"].dataframe = var_list["e0_T"].dataframe.iloc[:2]
        var_list["e0_c_i1"].lower_bound = 0
        var_list["e0_c_i1"].upper_bound = None

        pe = mopeds.ParameterEstimation(
            model, [var_list]
        )
        pe.setup_regularization(0, np.array([95000]))

        if model.DAE:
            answer_scaled = 1.26485e-13
        else:
            answer_scaled = 2.04875e-18

        res = pe.optimize()
        assert np.isclose(res["f"], ca.DM(answer_scaled))
        assert np.isclose(pe.optimize(reuse_solver=True)["f"], ca.DM(answer_scaled))

        res_tikh = pe.optimize(objective_function="tikh")
        logging.warning(
            f"Model.DAE: {model.DAE}, Result: {res['f']}, Expecting: {answer_scaled}"
        )
        assert np.isclose(res_tikh["f"], ca.DM(answer_scaled))
        assert np.isclose(pe.optimize(reuse_solver=True)["f"], ca.DM(answer_scaled))


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_pe_regularization(piecewise, dae, scaling):
    """Test that ParameterEstimation on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    with mopeds.options(variable_scaling=True):
        var_list, model = mopeds.examples.cstr(piecewise, dae, True)
        true_parameters = {}
        for n, v in var_list.items():
            if isinstance(v, mopeds.VariableParameter):
                true_parameters[n] = v.value[0]
                v.fixed = False
            elif isinstance(v, mopeds.VariableState):
                v.variance = 0.001**2

        time_grid = np.linspace(10, 100, 3)
        time_grid = np.insert(time_grid, 0, 0)

        var_list_fixed = copy.deepcopy(var_list)
        for var in var_list_fixed.values():
            var.fixed = True
        var_list_exp = mopeds.Simulator(
            model, time_grid, var_list_fixed
        ).simulate()[2]

        for key, var in var_list_exp.items():
            if isinstance(var, mopeds.VariableControlPiecewiseConstant):
                var_list[key].variable_list = var.variable_list
            else:
                var_list[key].dataframe = var.dataframe

        pe = mopeds.ParameterEstimation(
            model, [var_list]
        )
        a = pe.parameter_identifiability_chu2012(true_parameters, true_parameters.keys())
        b = pe.parameter_identifiability_yao2003(true_parameters, true_parameters.keys())
        c = pe.parameter_identifiability_lopez2013(true_parameters, true_parameters.keys())
        d = pe.parameter_identifiability_quaiser2009(true_parameters, true_parameters.keys())

        with pytest.raises(NotImplementedError):
            e = pe.parameter_identifiability_brun2001(true_parameters, true_parameters.keys())

        with mopeds.options(variable_scaling=False):
            pe = mopeds.ParameterEstimation(model, [var_list])
            e = pe.parameter_identifiability_brun2001(true_parameters, true_parameters.keys())

        identifiable_a = ['e0_E_r2', 'e0_c_p']
        identifiable_b = ['e0_E_r1', 'e0_E_r3', 'e0_c_p']
        identifiable_d = ['e0_c_p', 'e0_E_r2', 'e0_E_r3']
        ranked_c = ['e0_c_p', 'e0_E_r2', 'e0_E_r3', 'e0_k_pre_r2', 'e0_k_pre_r3', 'e0_E_r1', 'e0_greek_Deltah_r2', 'e0_k_pre_r1', 'e0_U', 'e0_greek_Deltah_r3', 'e0_greek_Deltah_r1']
        ranked_d = ['e0_c_p', 'e0_E_r2', 'e0_E_r3', 'e0_k_pre_r2', 'e0_k_pre_r3', 'e0_E_r1', 'e0_U', 'e0_greek_Deltah_r2', 'e0_k_pre_r1', 'e0_greek_Deltah_r3', 'e0_greek_Deltah_r1']
        ranked_e = ['e0_greek_Deltah_r1', 'e0_greek_Deltah_r3', 'e0_greek_Deltah_r2', 'e0_k_pre_r1', 'e0_E_r1', 'e0_k_pre_r3', 'e0_k_pre_r2', 'e0_E_r3', 'e0_E_r2', 'e0_U', 'e0_c_p']
        identifiable_e = ['e0_E_r1', 'e0_E_r3', 'e0_k_pre_r2', 'e0_U', 'e0_greek_Deltah_r1']

        assert a["estimable"] == identifiable_a
        assert b["estimable"] == identifiable_b
        assert c["estimable"] == identifiable_a
        assert d["estimable"] == identifiable_d
        assert e["estimable"] == identifiable_e

        if not c["ranked"] == ranked_c:
            warnings.warn("Ranking is scaling dependent")
        if not d["ranked"] == ranked_d:
            warnings.warn("Ranking is scaling dependent")
        if not e["ranked"] == ranked_e:
            warnings.warn("Ranking is scaling dependent")


def test_pe_algebraic():
    """Test if dynamic pe always returns same results and sensitivity, 
    independently if used with or without algebraic variables.
    """
    var_list, model, _ = mopeds.examples.yeast_growth(piecewise=False, duplicate_states=True)
    sim = mopeds.Simulator(model, np.linspace(0, 20, 5), var_list)
    exp_data = sim.simulate(algebraic=True)[2]
    parameters = {"theta1": 0.3, "theta2": 0.17}

    variable_combinations = [
        ["x1_alg", "x2_alg"],
        ["x1", "x2"],
        ["x1", "x2_alg"]
    ]

    piecewise_options = [True, False]
    scaling_options = [False, True]

    benchmark_res_pe = None
    benchmark_sens = None

    for scaling in scaling_options:
        with mopeds.options(variable_scaling=scaling):
            for piecewise in piecewise_options:
                var_list, model, _ = mopeds.examples.yeast_growth(piecewise=piecewise, duplicate_states=True)
                for var_names in variable_combinations:
                    var_list_i = copy.deepcopy(var_list)
                    for var_name in var_names:
                        var_list_i[var_name] = exp_data[var_name]
                    var_list_i["theta1"].fixed = False
                    var_list_i["theta1"].guess = 1
                    var_list_i["theta2"].fixed = False
                    var_list_i["theta2"].guess = 1


                    pe = mopeds.ParameterEstimation(model, [var_list_i])

                    res_pe = pe.optimize()
                    if benchmark_res_pe is None:
                        benchmark_res_pe = res_pe
                    else:
                        for key in benchmark_res_pe.keys():
                            if key in ["lam_p", "lam_x", "x_dict", "x_dict_all", "x_unscaled"]:
                                continue
                            v1 = res_pe[key]
                            v2 = benchmark_res_pe[key]
                            # print(key)
                            # print(v1/v2)
                            assert np.isclose(v1, v2).all()

                    res_sens = pe.calculate_sensitivity_and_fim(parameters)
                    if benchmark_sens is None:
                        benchmark_sens = res_sens
                    else:
                        for key in benchmark_sens.keys():
                            if key in ["jac_sorted", "jac_scaled_sorted", "jac_yao_sorted", "jac_wls", "hess_wls", "hess_tikh"]:
                                continue
                            v1 = res_sens[key]
                            v2 = benchmark_sens[key]
                            # print(key)
                            # print(v1/v2)
                            assert np.isclose(v1, v2, rtol=1e-4).all()



if __name__ == "__main__":
    test_pe_algebraic()
    # test_pe(True,True)
    # test_pe_objective(False)
    # test_pe_intials_algebraic()
    # test_pe_regularization(True, True, True)
    # test_scaling(True, True)
