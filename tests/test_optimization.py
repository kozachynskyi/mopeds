import logging
import casadi as ca
import numpy as np

import copy
import mopeds.examples
import mopeds
import mopeds.tools
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
def test_parameter_jacobian(piecewise):
    for cstr_model in [
        mopeds.examples.cstr_ode,
        mopeds.examples.cstr_ode_constant,
        mopeds.examples.cstr_dae,
        mopeds.examples.cstr_dae_constant,
    ]:
        var_list, model = cstr_model(piecewise)
        time_grid = np.linspace(0, 1000, 4)
        time_grid_expanded = list(time_grid) + [2000, 4000]

        if piecewise:
            T_in = var_list["e0_T_in"]
            T_in.expand_horizon([2000, 4000], [373, 373])

        var_list_exp = mopeds.Simulator(model, time_grid, var_list).generate_exp_data()

        for key, var in var_list_exp.items():
            var_list[key] = var

        var_list["e0_U"].fixed = False
        var_list["e0_E_r1"].fixed = False

        pe = mopeds.ParameterEstimation(
            model, [var_list]
        )
        jac_pe = pe.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_E_r1": 9.6e4})["jac_scaled_full_theory"]

        oed = mopeds.OptimalExperimentalDesign(model, [var_list], time_grid)
        oed_expanded = mopeds.OptimalExperimentalDesign(model, [var_list], time_grid_expanded)
        jac_oed = oed.calculate_objective_and_jacobian({"e0_T_in": 373})["jac"]
        jac_oed_expanded = oed_expanded.calculate_objective_and_jacobian({"e0_T_in": 373})["jac"]

        if piecewise:
            with pytest.raises(ValueError):
                assert np.all(np.isclose(jac_pe, jac_oed))
            assert np.all(np.isclose(jac_pe, jac_oed_expanded))
        else:
            assert np.all(np.isclose(jac_pe, jac_oed))
            with pytest.raises(ValueError):
                assert np.all(np.isclose(jac_pe, jac_oed_expanded))


@pytest.mark.parametrize("piecewise", [True, False])
def test_optimizer(piecewise):  # noqa: C901
    """Tests if optimizer can deal with variable list of fixed and unfixed parameters.
    Not well designed, and may yield false positives, but let it be.
    """
    variable_list, m = mopeds.examples.cstr_ode(piecewise)

    for var in variable_list.values():
        if isinstance(
            var,
            (
                mopeds.VariableControl,
                mopeds.VariableParameter,
                mopeds.VariableControlPiecewiseConstant,
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
    var_list_exp = mopeds.Simulator(m, time_grid, var_list_fixed).generate_exp_data()

    # Replace empty state variables with results from simulation
    for key, var in var_list_exp.items():
        if isinstance(var, mopeds.VariableState):
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
                pe = mopeds.ParameterEstimation(m, [variable_list])
            else:
                pe = mopeds.ParameterEstimation(m, [variable_list, variable_list])

            pe.solver_settings = {
                "verbose": False,
                "ipopt": {
                    "hessian_approximation": "limited-memory",
                    "max_iter": 1,
                    "print_level": 0,
                },
            }

            oed = mopeds.OptimalExperimentalDesign(m, [variable_list], time_grid)
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
                        res_oed = oed.optimize(1e-3)

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
    test_optimizer(True)
    test_parameter_jacobian(True)
