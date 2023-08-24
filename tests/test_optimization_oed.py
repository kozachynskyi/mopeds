import logging
import casadi as ca
import numpy as np

import copy
import par_est.examples
import par_est
import par_est.tools
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
def test_jacobian_weights(piecewise):
    """Test if jacobian weights correctly implemented"""
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        time_grid = np.linspace(0, 1000, 4)
        time_grid_modified = np.delete(time_grid, [2])
        time_grid_expanded = list(time_grid) + [2000, 4000]

        for weight_on in [True, False]:
            var_list, model = cstr_model(piecewise)
            if piecewise:
                T_in = var_list["e0_T_in"]
                T_in.expand_horizon([2000, 4000], [373, 373])

            if weight_on:
                var_list_exp = par_est.Simulator(model, time_grid, var_list).generate_exp_data()
            else:
                var_list_exp = par_est.Simulator(model, time_grid_modified, var_list).generate_exp_data()

            for key, var in var_list_exp.items():
                var_list[key] = var

            var_list["e0_U"].fixed = False
            var_list["e0_E_r1"].fixed = False

            pe = par_est.ParameterEstimation(
                model, [var_list]
            )
            jac_pe = pe.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_E_r1": 9.6e4})["jac_scaled_full_theory"].flatten()

            controls_dict = {"e0_T_in": 373}
            for i in range(len(time_grid_expanded)):
                controls_dict[f"weight_{i}"] = 1

            if not weight_on:
                controls_dict["weight_1"] = 0

            oed_settings = par_est.OEDsettings(measurement_weights=True)
            oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid, oed_settings)
            oed_expanded = par_est.OptimalExperimentalDesign(model, [var_list], time_grid_expanded, oed_settings)

            jac_oed = oed.calculate_objective_and_jacobian(controls_dict)["jac"]
            jac_oed_expanded = oed_expanded.calculate_objective_and_jacobian(controls_dict)["jac"]

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

        oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
        res = oed.optimize()

        logging.warning(f"{res['f']}")
        assert np.isclose(res["f"], ca.DM(39.499), rtol=0, atol=1.0e-4)

        # For not functionality is turnded off
        # res = oed.optimize(True)

        # logging.warning(f"{res['f']}")
        # assert np.isclose(res["f"], ca.DM(45.1675), rtol=0, atol=1.0e-4)

def test_oed_piecewise():
    """Test that OptimalExperimentalDesign on ODE and DAE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_ode_constant,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
    ]:
        var_list_peicewise, model_piecewise = cstr_model(piecewise_control=True)
        var_list, model = cstr_model(piecewise_control=False)
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)

        var_list["e0_E_r1"].fixed = False
        var_list["e0_T_in"].fixed = True
        var_list["e0_c_in_i1"].fixed = False

        var_list_peicewise["e0_E_r1"].fixed = False
        var_list_peicewise["e0_T_in"].fixed = True
        var_list_peicewise["e0_c_in_i1"].fixed = False

        oed_settings = par_est.OEDsettings(4, num_control_switches=0)

        oed_piecewise = par_est.OptimalExperimentalDesign(model_piecewise, [var_list_peicewise], time_grid, oed_settings)
        res_piecewise = oed_piecewise.optimize()

        oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid, oed_settings)
        res = oed.optimize()

        exp_data = oed_piecewise.generate_experimental_data({"e0_c_in_i1_t0": 5})
        pe = par_est.ParameterEstimation(model_piecewise, [exp_data])
        res_pe = pe.optimize()
        assert np.isclose(oed_piecewise.parameter_values, res_pe["x"])

        assert np.isclose(res["f"], res_piecewise["f"])
        # logging.warning(f"{res['f']}")
        # assert np.isclose(res["f"], ca.DM(39.499), rtol=0, atol=1.0e-4)


if __name__ == "__main__":
    pass
    test_jacobian_weights(False)
    # test_oed(True)
    # test_parameter_jacobian(True)
    # test_oed_piecewise()
    # test_oed_piecewise()
