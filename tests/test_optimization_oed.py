import logging
import casadi as ca
import numpy as np

import copy
import mopeds.examples
import mopeds
import mopeds.tools
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_jacobian_weights(piecewise, dae, use_constant, scaling):
    """Test if jacobian weights correctly implemented"""
    with mopeds.options(variable_scaling=scaling):
        time_grid = np.linspace(0, 1000, 4)
        time_grid_modified = np.delete(time_grid, [2])
        time_grid_expanded = list(time_grid) + [2000, 4000]

        for weight_on in [True, False]:
            var_list, model = mopeds.examples.cstr(piecewise, dae, use_constant)
            if piecewise:
                T_in = var_list["e0_T_in"]
                T_in.expand_horizon([2000, 4000], [373, 373])

            if weight_on:
                var_list_exp = mopeds.Simulator(model, time_grid, var_list).simulate()[
                    2
                ]
            else:
                var_list_exp = mopeds.Simulator(
                    model, time_grid_modified, var_list
                ).simulate()[2]

            for key, var in var_list_exp.items():
                var_list[key] = var

            var_list["e0_U"].fixed = False
            var_list["e0_E_r1"].fixed = False

            # Needed to compensate for the fact that OED internally scaled parameters to -1 and 1
            var_list["e0_U"].lower_bound = -1.4
            var_list["e0_U"].upper_bound = 1.4
            var_list["e0_E_r1"].lower_bound = -9.6e4
            var_list["e0_E_r1"].upper_bound = 9.6e4

            pe = mopeds.ParameterEstimation(model, [var_list])
            jac_pe = pe.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_E_r1": 9.6e4})[
                "jac_scaled_full_theory"
            ].flatten()

            controls_dict = {
                "e0_T_in": 373,
                "weight_0": 1,
                "weight_1": 1,
                "weight_2": 1,
            }

            controls_dict_expanded = copy.deepcopy(controls_dict)
            controls_dict_expanded["weight_3"] = 1
            controls_dict_expanded["weight_4"] = 1

            oed_settings = mopeds.OptimalSampling(num_sampling_times=3)
            if not weight_on:
                controls_dict["weight_1"] = 0
                controls_dict_expanded["weight_1"] = 0

            oed_settings = mopeds.OptimalSampling(num_sampling_times=3)
            oed = mopeds.OptimalExperimentalDesign(
                model, [var_list], time_grid, oed_settings
            )
            oed_expanded = mopeds.OptimalExperimentalDesign(
                model, [var_list], time_grid_expanded, oed_settings
            )

            jac_oed = oed.calculate_objective_and_jacobian(controls_dict)["jac"]
            jac_oed_expanded = oed_expanded.calculate_objective_and_jacobian(
                controls_dict_expanded
            )["jac"]

            jac_oed = jac_oed[jac_oed != 0]
            jac_oed_expanded = jac_oed_expanded[jac_oed_expanded != 0]

            if piecewise:
                with pytest.raises(ValueError):
                    assert np.all(np.isclose(jac_pe, jac_oed))
                assert np.all(np.isclose(jac_pe, jac_oed_expanded))
            else:
                assert np.all(np.isclose(jac_pe, jac_oed))
                with pytest.raises(ValueError):
                    assert np.all(np.isclose(jac_pe, jac_oed_expanded))


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_oed(piecewise, dae, use_constant, scaling):
    """Test that OptimalExperimentalDesign on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    var_list, model = mopeds.examples.cstr(piecewise, dae, use_constant)
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)
    for var in var_list.values():
        var.fixed = True

    var_list["e0_E_r1"].fixed = False
    var_list["e0_T_in"].fixed = True
    var_list["e0_c_in_i1"].fixed = False

    with mopeds.options(variable_scaling=scaling):
        oed = mopeds.OptimalExperimentalDesign(model, [var_list], time_grid)
        res = oed.optimize()

    logging.warning(f"{res['f']}")
    expected = ca.DM(39.4989)

    assert np.isclose(res["f"], expected)
    assert np.isclose(res["x"], 6, rtol=1e-4)


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
def test_scaling(piecewise, dae, use_constant):
    """Test that OED methods yeilds same results scaled and unscaled"""
    var_list, model = mopeds.examples.cstr(piecewise, dae, use_constant)
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)
    for var in var_list.values():
        var.fixed = True

    var_list["e0_E_r1"].fixed = False
    var_list["e0_T_in"].fixed = True
    var_list["e0_c_in_i1"].fixed = False

    if piecewise:
        controls = {"e0_c_in_i1_t0": 5}
    else:
        controls = {"e0_c_in_i1": 5}

    res = [[], []]
    for scaling in [True, False]:
        with mopeds.options(variable_scaling=scaling):
            oed = mopeds.OptimalExperimentalDesign(model, [var_list], time_grid)
            res[0].append(oed.calculate_objective_and_jacobian(controls)["jac"])
            res[1].append(oed.generate_experimental_data(controls).dataframe)

    for v in res:
        assert np.isclose(v[0], v[1], equal_nan=True).all()


@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_oed_piecewise(dae, use_constant, scaling):
    """Test that OptimalExperimentalDesign on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    with mopeds.options(variable_scaling=scaling):
        var_list_peicewise, model_piecewise = mopeds.examples.cstr(
            True, dae, use_constant
        )
        var_list, model = mopeds.examples.cstr(False, dae, use_constant)
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)

        var_list["e0_E_r1"].fixed = False
        var_list["e0_T_in"].fixed = True
        var_list["e0_c_in_i1"].fixed = False

        var_list_peicewise["e0_E_r1"].fixed = False
        var_list_peicewise["e0_T_in"].fixed = True
        var_list_peicewise["e0_c_in_i1"].fixed = False

        oed_settings = None

        oed_piecewise = mopeds.OptimalExperimentalDesign(
            model_piecewise, [var_list_peicewise], time_grid, oed_settings
        )
        res_piecewise = oed_piecewise.optimize()

        oed = mopeds.OptimalExperimentalDesign(
            model, [var_list], time_grid, oed_settings
        )
        res = oed.optimize()

        exp_data = oed_piecewise.generate_experimental_data({"e0_c_in_i1_t0": 5})
        pe = mopeds.ParameterEstimation(model_piecewise, [exp_data])
        res_pe = pe.calculate_objective_and_residual({"e0_E_r1": 9.6e4})

        assert np.isclose(res_pe["f"], 0)
        assert np.isclose(res["f"], res_piecewise["f"])


if __name__ == "__main__":
    pass
    # test_jacobian_weights(True, True, True, False)
    # test_oed(True, False, True, False)
    # test_oed_piecewise(True, True, True)
    test_scaling(True, True, True)
