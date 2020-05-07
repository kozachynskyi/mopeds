import par_est
import par_est.tools
import logging
import casadi as ca
import numpy as np
import pytest

import copy
from conftest import cstr_model_ode, cstr_model_dae, pendulum_dae_1
import logging


# @pytest.mark.skip(reason="WIP")
def test_pe():
    for cstr_model in [cstr_model_ode, cstr_model_dae]:
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
        breakpoint()

        if model.DAE:
            answer_scaled = 1.26485e-13
            answer = 9.75963e-12
        else:
            answer_scaled = 8.39521e-15
            answer = 9.26777e-12

        res = pe.optimize()
        logging.warning(
            f"Model.DAE: {model.DAE}, Result: {res['f']}, Expecting: {answer_scaled}"
        )
        assert np.isclose(res["f"], ca.DM(answer_scaled), rtol=0, atol=1.0e-13)

        res = pe.optimize(False)
        logging.warning(
            f"Model.DAE: {model.DAE}, Result: {res['f']}, Expecting: {answer}"
        )
        assert np.isclose(res["f"], ca.DM(answer), rtol=0, atol=1.0e-12)


def test_covariance_manipulation():
    """ Covariance calculation produces square matrix with size (NumOfParam * NumOfTimepoints)
    Something in matrix multiplication makes covariance_difference_fail to be not excatly zero. """

    num_time = 10
    varlist, model = cstr_model_ode()

    time_grid = np.linspace(0, 1000, num_time + 1)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res = sim.simulate(True)

    num_param = 19
    num_state = 5
    covariance_measurement = np.eye(num_state)

    res_jacobian = res["jac_xf_p"]
    covariance_full = res_jacobian.T @ covariance_measurement @ res_jacobian

    covariance_measurement_reshape = np.kron(
        np.eye(num_time, dtype=int), covariance_measurement
    )

    split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param

    list_jacobian_at_timepoint = ca.horzsplit(res_jacobian, split_vector)

    jacobian_reshape = None
    covariance_full_sum = None

    for count_time_step in range(num_time):
        index_from = count_time_step * num_param
        index_till = num_param * count_time_step + (num_param)

        jac_at_timepoint = res_jacobian[:, index_from:index_till][:, 1:]
        cov_at_timepoint = (
            jac_at_timepoint.T @ covariance_measurement @ jac_at_timepoint
        )

        assert (
            list_jacobian_at_timepoint[count_time_step][:, 1:] - jac_at_timepoint
        ).is_zero()

        if jacobian_reshape is None:
            jacobian_reshape = jac_at_timepoint
        else:
            jacobian_reshape = ca.vertcat(jacobian_reshape, jac_at_timepoint)

        covariance_difference = (
            covariance_full[index_from + 1 : index_till, index_from + 1 : index_till]
            - cov_at_timepoint
        )
        assert covariance_difference.is_zero()

        if covariance_full_sum is None:
            covariance_full_sum = cov_at_timepoint
        else:
            covariance_full_sum = covariance_full_sum + cov_at_timepoint

    covariance_reshape = (
        jacobian_reshape.T @ covariance_measurement_reshape @ jacobian_reshape
    )

    covariance_difference_fail = covariance_full_sum - covariance_reshape
    assert not covariance_difference_fail.is_zero()


@pytest.mark.skip(reason="WIP")
def test_oed():
    var_list, model = cstr_model_ode()
    time_grid = np.linspace(10, 100, 3)
    time_grid = np.insert(time_grid, 0, 0)
    for var in var_list.values():
        var.fixed = True

    var_list["e0_E_r1"].fixed = False
    var_list["e0_c_in_i1"].fixed = False

    oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
    # fim = oed.get_fim_matrix()
    # print(oed.list_simulators[0].simulate())
    res = oed.optimize()
    res = oed.optimize()
    res = oed.optimize()
    res = oed.optimize()

    # print(oed.list_simulators[0].simulate())
    # breakpoint()
    # assert fim[0].size() == (20, 1)
    # assert fim[1].size() == (1, 1)

    # logging.warning(f"{res['f']}")
    # assert np.isclose(res["f"], ca.DM(45.16749745), rtol=0, atol=1.0e-8)

    # res = oed.optimize(True)
    # print(oed.list_simulators[0].simulate())
    # breakpoint()

    # logging.warning(f"{res['f']}")
    # assert np.isclose(res["f"], ca.DM(45.16749745), rtol=0, atol=1.0e-8)
    breakpoint()


def not_test_optimizer():
    variable_list, m = cstr_model_ode()

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
                "ipopt": {"max_iter": 1},
            }

            oed = par_est.OptimalExperimentalDesign(m, [variable_list], time_grid)
            oed.solver_settings = {
                "verbose": False,
                "ipopt": {"hessian_approximation": "limited-memory", "max_iter": 1},
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

            # oed.get_fim_matrix()

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


def test_scaling():
    var_list, m = cstr_model_dae()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    var_list_fixed = copy.deepcopy(var_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).generate_exp_data()

    for key, var in var_list_exp.items():
        var_list[key] = var

    for var in var_list.values():
        var.fixed = True

    var_list["e0_E_r1"].fixed = False

    pe = par_est.ParameterEstimation(m, [var_list])

    pe._setup_scaling(False)
    result_unscaled = pe.list_simulators[0].simulate(True)
    ev_unscaled = ca.Function(
        "aha",
        [pe.varlist_decision.get_casadi_var()],
        [result_unscaled["xf"], result_unscaled["jac_xf_p"]],
    )

    res1 = ev_unscaled(90000)

    pe._setup_scaling(True)
    result_scaled = pe.list_simulators[0].simulate(True)
    ev_scaled = ca.Function(
        "aha",
        [pe.varlist_decision.get_casadi_var()],
        [result_scaled["xf"], result_scaled["jac_xf_p"]],
    )
    res2 = ev_scaled(1)

    print(res2[0] - res1[0])

    difference_jacobian = res2[1] - res1[1]
    print(difference_jacobian[:, 0:19])
    print(difference_jacobian[:, 19:38])
    print(difference_jacobian[:, 38:57])
    print(difference_jacobian[:, 57:76])


def test_jacobian_manipulation():
    """ Covariance calculation produces square matrix with size (NumOfParam * NumOfTimepoints)
    Something in matrix multiplication makes covariance_difference_fail to be not excatly zero. """

    num_time = 3
    varlist, model = cstr_model_ode()

    time_grid = np.linspace(0, 1000, num_time + 1)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res = sim.simulate(True)

    num_param = 19
    num_state = 5
    covariance_measurement = np.eye(num_state)

    res_jacobian = res["jac_xf_p"]

    index_all_states = list(range(len(model.varlist_state)))


    split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param

    list_jacobian_at_timepoint = ca.horzsplit(res_jacobian, split_vector)

    for jacobian in list_jacobian_at_timepoint:
        selected = jacobian.get(False,index_all_states,[0,2])
        par_est.tools.plot_arrays([jacobian,selected])




if __name__ == "__main__":
    # test_ode_oed()
    # not_test_optimizer()
    # test_pe()
    # test_scaling()
    import timeit
    aa = timeit.timeit(test_oed,number=1)
    breakpoint()
    # test_covariance_manipulation()
    # test_jacobian_manipulation()
