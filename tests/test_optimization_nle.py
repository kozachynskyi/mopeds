import copy
import logging

import casadi as ca
import numpy as np

import par_est
import par_est.examples
import par_est.tools


def test_pe():
    """Test that ParameterEstimationNLE on NLE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    variable_list, model = par_est.examples.vle_nle_problem()

    var_list_fixed = copy.deepcopy(variable_list)
    var_list_fixed.set_variable_list_fixed()

    control_bounds = {"x": [0.5, 0.5, 1]}

    (
        variable_list_optimizer,
        true_parameters,
    ) = par_est.tools.generate_varlist_with_data_NLE(
        model, variable_list, control_bounds, preturbate=False
    )
    variable_list_optimizer = variable_list_optimizer[0]
    variable_list_optimizer.set_variable_list_unfixed(["a2"])
    variable_list_optimizer.set_bounds(emerg_val=50)

    for i in range(3):
        if i == 0:
            simulator_name = "rootfinder"
        else:
            simulator_name = "ipopt"

        if i == 2:
            for var in variable_list_optimizer.values():
                if isinstance(var, par_est.VariableAlgebraic):
                    var.lower_bound = 350
                    var.upper_bound = 380

        pe = par_est.optimization.ParameterEstimationNLE(
            model,
            [variable_list_optimizer, variable_list_optimizer],
            simulator_name=simulator_name,
        )

        for switch in [True, False]:
            res = pe.optimize(switch)
            answer_f = 0
            answer_param = [5.19625]

            logging.warning(
                f"Model.NLE: {model}, Result: {res['f']}, Expecting: {answer_f}"
            )
            assert np.isclose(res["f"], ca.DM(answer_f), rtol=0, atol=1.0e-9)

            logging.warning(
                f"Model.NLE: {model}, Result: {res['x']}, Expecting: {answer_param}"
            )
            assert np.all(np.isclose(res["x"], ca.DM(answer_param), rtol=0, atol=1.0e-9))

        for sim in pe.list_simulators:
            if i == 1:
                assert sim._lower_bound[0] == -ca.inf
                assert sim._upper_bound[0] == ca.inf
                assert sim._rootfinder_bounds[0] == 0
            if i == 2:
                assert sim._lower_bound[0] == 350
                assert sim._upper_bound[0] == 380
                assert sim._rootfinder_bounds[0] == 2


def test_multivariate_pe():
    varlist, model = par_est.examples.simple_mixer()

    control_bounds = {"e0_F_s1": [0, 20, 3], "e0_F_s3": [10, 15, 3]}

    varlist["e0_F_s3"].fixed = False

    rng = np.random.default_rng(0)
    (
        variable_list_optimizer,
        true_parameters,
    ) = par_est.tools.generate_varlist_with_data_NLE(
        model, varlist, control_bounds, preturbate=True, rng=rng
    )

    mask_full = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    )

    data_full = np.array(
        [
            [0.12573022, -10.13210486, 0.64042265],
            [10.10490012, -0.53566937, 10.36159505],
            [21.30400005, 0, 19.29626476],
            [-1.26542147, -13.12327446, 0.04132598],
            [7.67496923, -2.71879166, 8.75408905],
            [19.26773265, 6.95574102, 19.68369984],
            [0.41163054, -13.95748663, -0.12853466],
            [11.36646347, -5.66519467, 0],
            [20.90347018, 5.0940123, 19.25650075],
        ]
    )

    var4 = variable_list_optimizer[2]["e0_F_s4"]
    var4.dataframe = var4._dataframe_from_value(None)

    var5 = variable_list_optimizer[7]["e0_F_s5"]
    var5.dataframe = var5._dataframe_from_value(None)

    parameters = {"e0_F_s3": 4}

    ols_f = 690.9195635362676
    ols_residuals = np.array(
        [
            [-0.12573022, 6.13210486, -0.64042265],
            [-0.10490012, 6.53566937, -0.36159505],
            [-1.30400005, 0.0, 0.70373524],
            [1.26542147, 9.12327446, -0.04132598],
            [2.32503077, 8.71879166, 1.24591095],
            [0.73226735, 9.04425898, 0.31630016],
            [-0.41163054, 9.95748663, 0.12853466],
            [-1.36646347, 11.66519467, 0.0],
            [-0.90347018, 10.9059877, 0.74349925],
        ]
    )

    pe = par_est.ParameterEstimationNLE(model, variable_list_optimizer)
    assert pe.names_of_measurements == ["e0_F_s2", "e0_F_s4", "e0_F_s5"]
    assert np.array_equal(pe.array_data_mask_new, mask_full)
    assert np.allclose(pe.array_data_new, data_full, equal_nan=True)
    assert pe.index_measurements_in_sim == [0, 1, 2]

    ols = pe.calculate_ols_value(parameters)
    assert np.isclose(ols_f, ols["f"])
    assert np.allclose(ols_residuals, ols["residuals"])

    # Remove one measurement variable
    for varlist in variable_list_optimizer:
        varlist["e0_F_s2"].dataframe = varlist["e0_F_s2"]._dataframe_from_value(None)

    pe = par_est.ParameterEstimationNLE(model, variable_list_optimizer)
    assert pe.names_of_measurements == ["e0_F_s4", "e0_F_s5"]
    assert np.array_equal(pe.array_data_mask_new, mask_full[:, 1:])
    assert np.allclose(pe.array_data_new, data_full[:, 1:], equal_nan=True)
    assert pe.index_measurements_in_sim == [1, 2]

    ols_f = 678.79613973
    ols_residuals = np.array(
        [
            [6.13210486, -0.64042265],
            [6.53566937, -0.36159505],
            [0.0, 0.70373524],
            [9.12327446, -0.04132598],
            [8.71879166, 1.24591095],
            [9.04425898, 0.31630016],
            [9.95748663, 0.12853466],
            [11.66519467, 0.0],
            [10.9059877, 0.74349925],
        ]
    )
    ols = pe.calculate_ols_value(parameters)
    assert np.isclose(ols_f, ols["f"])
    assert np.allclose(ols_residuals, ols["residuals"])


if __name__ == "__main__":
    pass
    test_pe()
    # test_multivariate_pe()
