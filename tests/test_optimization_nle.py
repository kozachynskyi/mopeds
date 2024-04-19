import copy
import logging
import pytest

import casadi as ca
import numpy as np

import mopeds
import mopeds.examples
import mopeds.tools


def test_parameter_jacobian():
    var_list, model = mopeds.examples.cstr_nle()

    list_var_list = []
    measurement_names = ["e0_k", "e0_T", "e0_c_c1"]
    for m_cat in [22, 18, 20]:
        var_list_i = copy.deepcopy(var_list)
        var_list_i["e0_m_Cat"].value = m_cat
        var_list_exp = mopeds.SimulatorNLE(model, var_list).simulate()[2]
        var_list_i["e0_cp"].fixed = False
        var_list_i["e0_E"].fixed = False

        for key, var in var_list_exp.items():
            if isinstance(var, mopeds.VariableAlgebraic):
                if key in measurement_names:
                    var_list_i[key].value = var.value[0]
        list_var_list.append(var_list_i)

    pe = mopeds.ParameterEstimationNLE(model, list_var_list)
    jac_pe = pe.calculate_sensitivity_and_fim({"e0_cp": 1.75, "e0_E": 135518.2})["jac_scaled_full_theory"]

    list_var_list[0]["e0_m_Cat"].fixed = False
    list_var_list[0]["e0_Q_Feed"].fixed = False

    oed = mopeds.OptimalExperimentalDesign_NLE(model, [list_var_list[0]], previous_measurements=[{"e0_m_Cat": 22, "e0_Q_Feed":30}, {"e0_m_Cat": 18, "e0_Q_Feed":30}], measurable_variables=measurement_names)

    jac_oed = oed.calculate_objective_and_jacobian({"e0_m_Cat": 20, "e0_Q_Feed":30})["jac"]
    
    rearange_index = [0, 3, 6,1,4,7,2,5,8]

    assert np.all(np.isclose(jac_pe, jac_oed[rearange_index, :]))

def test_pe():
    """Test that ParameterEstimationNLE on NLE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    variable_list, model = mopeds.examples.vle_nle_problem()

    variable_list["a2"].fixed = False
    variable_list.set_bounds(emerg_val=50)
    control_bounds = {"x": [0.5, 0.5, 1]}

    (
        variable_list_optimizer,
        true_parameters,
    ) = mopeds.tools.generate_artificial_data_from_grid_nle(
        model, variable_list, control_bounds, perturbate=False
    )
    variable_list_optimizer = variable_list_optimizer[0]

    for i in range(3):
        if i == 0:
            simulator_name = "rootfinder"
        else:
            simulator_name = "ipopt"

        if i == 2:
            for var in variable_list_optimizer.values():
                if isinstance(var, mopeds.VariableAlgebraic):
                    var.lower_bound = 350
                    var.upper_bound = 380

        pe = mopeds.optimization.ParameterEstimationNLE(
            model,
            [variable_list_optimizer, variable_list_optimizer],
            simulator_name=simulator_name,
        )

        for objective in ["ols", "wls", "fair", "tikh"]:
            for direct in [True, False]:
                res = pe.optimize(objective_function=objective, direct_optimization=direct)
                answer_f = 0
                answer_param = [5.19625]

                logging.warning(
                    f"Model.NLE: {model}, Result: {res['f']}, Expecting: {answer_f}"
                )
                assert np.isclose(res["f"], ca.DM(answer_f), rtol=0, atol=1.0e-9)
                assert np.isclose(pe.optimize(direct_optimization=direct, reuse_solver=True)["f"], ca.DM(answer_f), rtol=0, atol=1.0e-9)

                logging.warning(
                    f"Model.NLE: {model}, Result: {res['x']}, Expecting: {answer_param}"
                )
                assert np.all(
                    np.isclose(res["x"], ca.DM(answer_param), rtol=0, atol=1.0e-9)
                )
                assert np.all(
                    np.isclose(pe.optimize(direct_optimization=direct, reuse_solver=True)["x"], ca.DM(answer_param), rtol=0, atol=1.0e-9)
                )

        for sim in pe.list_simulators:
            if i == 1:
                assert sim._lower_bound[0] == -ca.inf
                assert sim._upper_bound[0] == ca.inf
            if i == 2:
                assert sim._lower_bound[0] == -1
                assert sim._upper_bound[0] == 1


def test_multivariate_pe():
    varlist, model = mopeds.examples.simple_mixer()

    control_bounds = {"e0_F_s1": [0, 20, 9]}

    varlist["e0_F_s3"].fixed = False

    varlist["e0_F_s2"].variance = 0.1
    varlist["e0_F_s4"].variance = 0.1

    rng = np.random.default_rng(0)
    (
        variable_list_optimizer,
        true_parameters,
    ) = mopeds.tools.generate_artificial_data_from_grid_nle(
        model, varlist, control_bounds, perturbate=True, rng=rng
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
            [0.03975939, -13.04177523, 0.64042265],
            [2.53317233, -10.66939353, 2.86159505],
            [5.41236102, 0.0, 4.29626476],
            [7.0998386, -5.69709669, 7.54132598],
            [9.26476071, -3.069188, 8.75408905],
            [12.26843673, -0.6721098, 12.18369984],
            [15.130169, 2.32967167, 14.87146534],
            [17.93211369, 4.28964697, 0.0],
            [20.28570236, 7.0297293, 19.25650075],
        ]
    )

    var4 = variable_list_optimizer[2]["e0_F_s4"]
    var4.dataframe = var4._dataframe_from_value(None)

    var5 = variable_list_optimizer[7]["e0_F_s5"]
    var5.dataframe = var5._dataframe_from_value(None)

    parameters = {"e0_F_s3": 4}

    ols_f = 661.7386240718675
    wls_f = 6588.050781719137
    ols_residuals = np.array(
        [
            [-0.03975939, 9.04177523, -0.64042265],
            [-0.03317233, 9.16939353, -0.36159505],
            [-0.41236102, 0.0, 0.70373524],
            [0.4001614, 9.19709669, -0.04132598],
            [0.73523929, 9.069188, 1.24591095],
            [0.23156327, 9.1721098, 0.31630016],
            [-0.130169, 8.67032833, 0.12853466],
            [-0.43211369, 9.21035303, 0.0],
            [-0.28570236, 8.9702707, 0.74349925],
        ]
    )
    jac_full = np.array(
        [
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [-1.0],
            [-1.0],
            [-0.0],
            [-1.0],
            [-1.0],
            [-1.0],
            [-1.0],
            [-1.0],
            [-1.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
        ]
    )

    jac_scaled_full = np.array(
        [
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [-9.78717763],
            [-9.78717763],
            [-0.0],
            [-9.78717763],
            [-9.78717763],
            [-9.78717763],
            [-9.78717763],
            [-9.78717763],
            [-9.78717763],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0],
        ]
    )
    ols_y = np.array(
        [
            [0.0, -4.0, 0.0],
            [2.5, -1.5, 2.5],
            [5.0, 1.0, 5.0],
            [7.5, 3.5, 7.5],
            [10.0, 6.0, 10.0],
            [12.5, 8.5, 12.5],
            [15.0, 11.0, 15.0],
            [17.5, 13.5, 17.5],
            [20.0, 16.0, 20.0],
        ]
    )

    fim = 8.0
    fim_scaled = 80
    cov_par = 0.0125
    jac_wls = -10.0103
    hess_wls = 160

    pe = mopeds.ParameterEstimationNLE(model, variable_list_optimizer)
    assert pe.names_of_measurements == ["e0_F_s2", "e0_F_s4", "e0_F_s5"]
    assert np.array_equal(pe.array_data_mask, mask_full)
    assert np.allclose(pe.array_data, data_full, equal_nan=True)
    assert pe.index_measurements_in_sim == [0, 1, 2]

    ols = pe.calculate_objective_and_residual(parameters, "ols")
    wls = pe.calculate_objective_and_residual(parameters, "wls")
    assert np.isclose(ols_f, ols["f"])
    assert np.allclose(ols_residuals, ols["residuals"])
    assert np.allclose(ols_y, ols["y"])
    assert np.isclose(wls_f, wls["f"])

    res_sens = pe.calculate_sensitivity_and_fim(true_parameters)
    assert np.array_equal(res_sens["jac_full"], jac_full)
    assert np.allclose(res_sens["jac_scaled_full"], jac_scaled_full)
    assert np.isclose(res_sens["cov_par"], cov_par)
    assert np.isclose(res_sens["jac_wls"], jac_wls)
    assert np.isclose(res_sens["fim"], fim)
    assert np.isclose(res_sens["fim_scaled"], fim_scaled)
    assert np.isclose(res_sens["hess_wls"], hess_wls)

    # Remove one measurement variable
    for varlist in variable_list_optimizer:
        varlist["e0_F_s2"].dataframe = varlist["e0_F_s2"]._dataframe_from_value(None)

    pe = mopeds.ParameterEstimationNLE(model, variable_list_optimizer)
    assert pe.names_of_measurements == ["e0_F_s4", "e0_F_s5"]
    assert np.array_equal(pe.array_data_mask, mask_full[:, 1:])
    assert np.allclose(pe.array_data, data_full[:, 1:], equal_nan=True)
    assert pe.index_measurements_in_sim == [1, 2]

    ols_f = 660.5262816912779
    wls_f = 6575.9273579132405
    ols_residuals = np.array(
        [
            [9.04177523, -0.64042265],
            [9.16939353, -0.36159505],
            [0.0, 0.70373524],
            [9.19709669, -0.04132598],
            [9.069188, 1.24591095],
            [9.1721098, 0.31630016],
            [8.67032833, 0.12853466],
            [9.21035303, 0.0],
            [8.9702707, 0.74349925],
        ]
    )
    jac_scaled_full = [
        [-7.75832712],
        [-7.75832712],
        [-0.0],
        [-7.75832712],
        [-7.75832712],
        [-7.75832712],
        [-7.75832712],
        [-7.75832712],
        [-7.75832712],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
        [0.0],
    ]
    ols = pe.calculate_objective_and_residual(parameters, "ols")
    wls = pe.calculate_objective_and_residual(parameters, "wls")
    assert np.isclose(ols_f, ols["f"])
    assert np.allclose(ols_residuals, ols["residuals"])
    assert np.allclose(ols_y[:, 1:], ols["y"])
    assert np.isclose(wls_f, wls["f"])

    cov_par = 0.0125

    res_sens = pe.calculate_sensitivity_and_fim(true_parameters)
    assert np.array_equal(res_sens["jac_full"], jac_full[9:])
    assert np.allclose(res_sens["jac_scaled_full"], jac_scaled_full)
    assert np.isclose(res_sens["cov_par"], cov_par)


def test_inference_bounds():
    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.bod_model()

    dict_of_params = {
        "theta1": 19.143,
        "theta2": 0.5311,
    }

    dict_of_controls = {
        "x": [0, 8, 5 + 1],
    }

    dict_of_responses = {
        "f": 1e1,
    }

    pe = mopeds.ParameterEstimationNLE(MODEL, EXP_DATA)
    exp_inference_results, exp_data, sim_data = pe.calculate_inference_bounds(
        dict_of_params,
        dict_of_responses,
        dict_of_controls,
    )

    R = np.array([[-1.95560506, -20.49667278], [0.0, -12.55188463]])
    simulation = np.array(
        [0.0, 10.95903204, 15.64421028, 17.64720608, 18.50352189, 18.8696119]
    )
    bound = np.array([0.0, 6.0464797, 4.44109375, 4.8561333, 6.55208547, 7.80504387])

    exp_data_expected = np.array([8.3, 10.3, 19.0, 16.0, 15.6, 19.8])
    sim_data_expected = np.array(
        [0.0, 10.95903204, 15.64421028, 17.64720608, 18.50352189, 18.8696119]
    )

    assert np.allclose(R, exp_inference_results["f"]["R"])
    assert np.allclose(simulation, exp_inference_results["f"]["simulation"])
    assert np.allclose(bound, exp_inference_results["f"]["bound"])
    assert np.allclose(exp_data["f"], exp_data_expected)
    assert np.allclose(sim_data["f"], sim_data_expected)

    exp_inference_results, exp_data, sim_data = pe.calculate_inference_bounds(
        dict_of_params,
        dict_of_responses,
        dict_of_controls,
        dict_of_controls,
        rng=np.random.default_rng(0),
    )
    exp_data_expected = np.array(
        [0.39759387, 10.54127978, 17.66940452, 17.97892938, 16.8095866, 20.01307587]
    )
    sim_data_expected = np.array(
        [0.0, 10.95903204, 15.64421028, 17.64720608, 18.50352189, 18.8696119]
    )
    R = np.array([[1.93684673, 15.16727704], [0.0, -11.82590011]])
    bound = np.array([0.0, 4.32448195, 3.22057785, 2.61916078, 3.18725966, 3.80380586])

    assert np.allclose(R, exp_inference_results["f"]["R"])
    assert np.allclose(sim_data_expected, exp_inference_results["f"]["simulation"])
    assert np.allclose(bound, exp_inference_results["f"]["bound"])
    assert np.allclose(exp_data["f"], exp_data_expected)
    assert np.allclose(sim_data["f"], sim_data_expected)


@pytest.mark.parametrize("meas_var", [None, ["y"]])
@pytest.mark.parametrize("criteria", ["A", "A_fd", "D"])
def test_oed_nle(meas_var, criteria):
    """Just compare the results of the optimization to some baseline.
    Variate scaling and direct_optimization, come to the same solution."""
    vl, m = mopeds.examples.linear_example()
    for var in vl.values():
        if isinstance(var, mopeds.VariableControl):
            var.fixed = False
        if isinstance(var, mopeds.VariableParameter):
            var.fixed = False
    vl["a"].fixed = True
    # vl["b"].fixed = True
    vl["y"].variance = 0.1
    vl["z"].variance = 0.3
    vl["q"].variance = 0.2

    results = []
    for scaling_flag in [True, False]:
        with mopeds.options(variable_scaling=scaling_flag):
            previous_meas = [
                {"v": 1.5, "u": 1.5},
                {"v": 2.5, "u": 4.5},
                {"v": 3.5, "u": 3.5},
                {"v": 4.5, "u": 3.5},
            ]
            oed = mopeds.OptimalExperimentalDesign_NLE(
                m,
                [vl],
                measurable_variables=meas_var,
                previous_measurements=previous_meas,
            )
            oed.solver_settings["ipopt"]["hessian_approximation"] = "exact"

            if criteria == "D":
                oed.objective_scaling = 1e10

            for direct_optimization_flag in [True, False]:
                results.append(oed.optimize(direct_optimization=direct_optimization_flag, objective_function=criteria))


    if meas_var is None:
        if "A" in criteria:
            expected_f = 0.001616787258285354
            expected_x = [2,4]
        elif "D" in criteria:
            expected_f = 0.02699111734101253
            expected_x = [1,4]
    else:
        if "A" in criteria:
            expected_f = 0.0471016387679191
            expected_x = [1,4]
        elif "D" in criteria:
            expected_f = 1.8867952943284358
            expected_x = [1,4]

    for res_i in results:
        res_f_i = res_i["f"]
        res_x_i = res_i["x"]
        # print("calc", oed.calculate_objective_and_jacobian(res_i["x_dict"])["f"])
        print("resf", float(res_f_i))
        print("resx", res_x_i)
        assert(np.isclose(expected_f, res_f_i))
        assert(np.isclose(expected_x, res_x_i.T, atol=1e-4).all())


if __name__ == "__main__":
    pass
    # test_pe()
    # test_multivariate_pe()
    # test_inference_bounds()
    test_parameter_jacobian()
    # test_oed_nle(None, "A")
    # test_oed_nle(["y"], "A")
    # test_oed_nle(None, "A_fd")
    # test_oed_nle(["y"], "A_fd")
    # test_oed_nle(["y"], "D")
    # test_oed_nle(None, "D")
