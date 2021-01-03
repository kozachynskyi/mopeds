import logging
import casadi as ca
import numpy as np

import copy
import par_est.examples
import par_est
import par_est.tools


def test_pe_intials_algebraic():
    variable_list1, model = par_est.examples.empy_dae()
    variable_list1["X1"].value.value = [1, 0]
    variable_list1["X2"].value.value = [0]
    variable_list1["X1"].value.time = [0, 1]
    variable_list1["X2"].value.time = [0]
    variable_list2 = copy.deepcopy(variable_list1)
    variable_list2["X1"].value.value = [1, 0]
    variable_list2["X2"].value.value = [1, 0]
    variable_list2["X1"].value.time = [0, 1]
    variable_list2["X2"].value.time = [0, 2]
    pe = par_est.ParameterEstimation(model, [variable_list1, variable_list2])
    pe.solver_settings = {
        "ipopt": {
            "max_iter": 0,
        },
    }
    assert pe.list_simulators[0]._initial_algebraic[0] == 0
    assert pe.list_simulators[1]._initial_algebraic[0] == 0
    pe._reinitialize_simulators_algebraic()
    assert pe.list_simulators[0]._initial_algebraic[0] == -1
    assert pe.list_simulators[1]._initial_algebraic[0] == -2


def test_pe_objective():
    variable_list1, model = par_est.examples.empy_dae()
    variable_list1["X1"].value.value = [0, 1]
    variable_list1["X2"].value.value = [0]
    variable_list1["X2"].variance = 10
    variable_list1["X1"].value.time = [0, 1]
    variable_list1["X2"].value.time = [0]
    variable_list2 = copy.deepcopy(variable_list1)
    variable_list2["X1"].value.value = [0, 3]
    variable_list2["X1"].variance = 20
    variable_list2["X2"].value.value = [0, 4]
    variable_list2["X2"].variance = 40
    variable_list2["X1"].value.time = [0, 1]
    variable_list2["X2"].value.time = [0, 2]
    pe = par_est.ParameterEstimation(model, [variable_list1, variable_list2])
    pe.solver_settings = {
        "ipopt": {
            "max_iter": 0,
        },
    }

    weight = np.array([1.5, 1.5, 1.0, 1.0, 1.0, 1.0])
    var = np.array([1.0, 0.1, 0.05, 0.025, 0.05, 0.025])
    data = np.array([1.0, 0.0, 3.0, 0.0, 0.0, 4.0])
    mask = np.array([1.0, 0.0, 1.0, 0.0, 0.0, 1.0])

    assert_numpy = np.testing.assert_array_equal
    assert_numpy(data, pe.array_data)
    assert_numpy(var, pe.inverted_variances)
    assert_numpy(weight, pe.experiments_weights)
    assert_numpy(mask, pe.array_data_mask)

    obj = np.sum(var * (data * mask) ** 2)
    obj_weight = np.sum(weight * var * (data * mask) ** 2)
    for switch in [True, False]:
        res = pe.optimize(switch)
        res_weight = pe.optimize(switch, scale_experiments=True)
        assert res["f"] == obj
        assert res_weight["f"] == obj_weight


def test_pe():
    """Test that ParameterEstimation on ODE and DAE always yields same result.
    Test that ParameterEstimationNLE on NLE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """

    """ ODE and DAE """
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        var_list, model = cstr_model()
        time_grid = np.linspace(10, 10000, 3)
        time_grid = np.insert(time_grid, 0, 0)

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

        var_list["e0_T"].value.value = var_list["e0_T"].value.value[0:2]
        var_list["e0_T"].value.time = var_list["e0_T"].value.time[0:2]

        pe = par_est.ParameterEstimation(model, [var_list])

        if model.DAE:
            answer_scaled = 1.26485e-13
            answer = 1.8412e-09
        else:
            answer_scaled = 2.04875e-18
            answer = 1.74658e-09

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

    """ NLE """
    variable_list, model = par_est.examples.vle_nle_problem()

    var_list_fixed = copy.deepcopy(variable_list)
    var_list_fixed.set_variable_list_fixed()

    var_list_fixed["x"].value = 0.5
    variable_list_optimizer = par_est.tools.generate_exp_data_list_NLE(
        model, var_list_fixed
    )
    variable_list_optimizer.set_variable_list_unfixed(["a2"])
    variable_list_optimizer.set_bounds(emerg_val=50)

    pe = par_est.optimization.ParameterEstimationNLE(
        model, [variable_list_optimizer, variable_list_optimizer]
    )
    res = pe.optimize(False)
    answer_f = 0
    answer_param = [5.19625]

    logging.warning(f"Model.NLE: {model}, Result: {res['f']}, Expecting: {answer_f}")
    assert np.isclose(res["f"], ca.DM(answer_f), rtol=0, atol=1.0e-9)

    logging.warning(
        f"Model.NLE: {model}, Result: {res['x']}, Expecting: {answer_param}"
    )
    assert np.all(np.isclose(res["x"], ca.DM(answer_param), rtol=0, atol=1.0e-9))


def test_oed():
    """Test that OptimalExperimentalDesign on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        var_list, model = cstr_model()
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)
        for var in var_list.values():
            var.fixed = True

        var_list["e0_E_r1"].fixed = False
        var_list["e0_c_in_i1"].fixed = False

        oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
        res = oed.optimize()

        logging.warning(f"{res['f']}")
        assert np.isclose(res["f"], ca.DM(45.1675), rtol=0, atol=1.0e-4)

        res = oed.optimize(True)

        logging.warning(f"{res['f']}")
        assert np.isclose(res["f"], ca.DM(45.1675), rtol=0, atol=1.0e-4)


def test_optimizer():
    """Tests if optimizer can deal with variable list of fixed and unfixed parameters.
    Not well designed, and may yield false positives, but let it be.
    """
    variable_list, m = par_est.examples.cstr_ode()

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
    # test_optimizer()
    # test_oed()
    # test_pe()
    # test_pe_objective()
    # test_pe_intials_algebraic()
