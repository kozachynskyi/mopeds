import par_est
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
    # Covariance calculation produces square matrix with size (NumOfParam * NumOfTimepoints)

    # a, b, c, d, e, f, g, h = [1e-50,2e50,3e30,4e-30,5e40,6e-75,7e-100,8e300]

    # m1 = np.array([[a, e],[ b,f],[c,g],[d,h]])
    # c1 = np.eye(2)
    # r1 = m1 @ c1 @ m1.T

    # m2 = np.array([[a,e,c,g],[b,f,d,h]])
    # c2 = np.eye(4)
    # r2 = m2 @ c2 @ m2.T

    # r1s = r1[0:2,0:2] + r1[2:4,2:4]
    # print(r1s - r2)





    num_time = 6
    varlist, model = pendulum_dae_1()

    time_grid = np.linspace(0, 1, num_time+1)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res = sim.simulate(True)


    num_param = 1
    num_state = 2
    c1 = np.eye(num_state)

    from scipy.linalg import block_diag

    resj = res["jac_xf_p"].toarray()
    cov = resj.T @ c1 @ resj

    c_reshape = np.kron(np.eye(num_time,dtype=int),c1)
    # c_reshape = ca.repmat(c, num_time, num_time)

    # cov_reshape = resj_reshape.T @ c_reshape @ resj_reshape

    from matplotlib import pyplot as plt
    import matplotlib.cm as cm

    # resj_reshape_temp = np.zeros((15,19))
    resj_reshape = None
    cov_sum = None
    for test_time in range(num_time):
        index_from = test_time * num_param
        index_till = num_param * test_time + (num_param)
        print(index_from, index_till)
        jac_at_timepoint = resj.T[index_from:index_till,:]
        cov_at_timepoint = jac_at_timepoint @ c1 @ jac_at_timepoint.T

        if resj_reshape is None:
            resj_reshape = resj[:,index_from:index_till]
        else:
            resj_reshape = ca.vertcat(resj_reshape, resj[:,index_from:index_till])
        # resj_reshape_temp[test_time * num_state : (test_time *num_state +num_state), 0:19] = resj[:,index_from:index_till]


        ddd = cov[index_from:index_till, index_from:index_till] - cov_at_timepoint
        # print(ddd)
        print(ddd)
        # fig = plt.figure()
        # fig.add_subplot(121).imshow(ddd, cmap=cm.Greens_r)
        # fig.add_subplot(122).imshow(jac_at_timepoint, cmap=cm.Greens_r)
        # plt.show()

        if cov_sum is None:
            cov_sum = cov_at_timepoint
        else:
            cov_sum = cov_sum + cov_at_timepoint

        # plt.imshow(cov_sum, cmap=cm.Greens_r)
        # fig = plt.figure()
        # fig.add_subplot(131).imshow(jac_at_timepoint, cmap=cm.Greens_r)
        # fig.add_subplot(132).imshow(cov_at_timepoint, cmap=cm.Greens_r)
        # fig.add_subplot(133).imshow(cov_sum, cmap=cm.Greens_r)
        # plt.show()


    cov_reshape = resj_reshape.T @ c_reshape @ resj_reshape
    
    diff =cov_sum - cov_reshape
    print(diff)

    # fig = plt.figure()
    # # fig.add_subplot(121).imshow(diff, cmap=cm.prism)
    # # fig.add_subplot(122).imshow(resj_reshape, cmap=cm.bwr)
    plt.show()

    plt.imshow(diff, cmap=cm.Greens_r)
    # plt.show()


    breakpoint()



@pytest.mark.skip(reason="WIP")
def test_oed():
    var_list, model = cstr_model_ode()
    time_grid = np.linspace(10, 10000, 3)
    time_grid = np.insert(time_grid, 0, 0)
    for var in var_list.values():
        var.fixed = True
    var_list["e0_c_in_i1"].fixed = False
    var_list["e0_E_r1"].fixed = False

    oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
    fim = oed.get_fim_matrix()
    res = oed.optimize()

    assert fim[0].size() == (20, 1)
    assert fim[1].size() == (1, 1)

    logging.warning(f"{res['f']}")
    assert np.isclose(res["f"], ca.DM(45.16749745), rtol=0, atol=1.0e-8)

    res = oed.optimize(False)

    logging.warning(f"{res['f']}")
    assert np.isclose(res["f"], ca.DM(45.16749745), rtol=0, atol=1.0e-8)


def not_test_optimizer():
    variable_list, m = cstr_model_ode()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 1000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    # Generate experimental data for parameter estimation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).get_symbolic_result()

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

            oed.get_fim_matrix()

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


if __name__ == "__main__":
    # test_ode_oed()
    # test_pe()
    # test_scaling()
    # test_oed()
    test_covariance_manipulation()
