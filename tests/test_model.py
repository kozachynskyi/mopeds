import pytest
import casadi as ca
import mopeds.examples
import numpy as np


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_model(piecewise, dae, use_constant, scaling):
    var_list, model = mopeds.examples.cstr(use_constant=True)

    assert len(model.varlist_state(var_list)) == 5
    assert len(model.varlist_independent(var_list)) == 18
    assert len(model.varlist(var_list)) == 34
    assert len(model.varlist_all) == 34

    assert model.equations_differential.size() == (5, 1)

    with pytest.raises(NotImplementedError):
        model.add_equations_differential(model.equations_differential[1])

    var_list, model = mopeds.examples.pendulum_dae_1()

    assert len(model.varlist_state(var_list)) == 2
    assert len(model.varlist_independent(var_list)) == 2
    assert len(model.varlist_algebraic(var_list)) == 3
    assert len(model.varlist(var_list)) == 7
    assert len(model.varlist_all) == 7

    assert model.equations_differential.size() == (2, 1)
    assert model.equations_algebraic.size() == (3, 1)

    with pytest.raises(NotImplementedError):
        model.add_equations_algebraic(model.equations_algebraic[1])

    var_list = mopeds.VariableList()
    var_list.add_variable(mopeds.Variable("a_test"))

    with mopeds.options(variable_scaling=scaling):
        var_list, model = mopeds.examples.cstr(piecewise, dae, use_constant)
        ode_system = {
            "x": model.varlist_state(var_list).get_casadi_variables(),
            "p": ca.vertcat(model.varlist_independent(var_list).get_casadi_variables()),
            "ode": model.equations_differential,
        }

        function = ca.Function(
            "eq_sys",
            [ode_system["x"], ode_system["p"]],
            [ode_system["ode"]],
            ["x", "p"],
            ["ode"],
        )

        if model.DAE:
            ode_system["z"] = model.varlist_algebraic(var_list).get_casadi_variables()
            ode_system["alg"] = model.equations_algebraic

            function = ca.Function(
                "eq_sys",
                [ode_system["x"], ode_system["z"], ode_system["p"]],
                [ode_system["ode"], ode_system["alg"]],
                ["x", "z", "p"],
                ["ode", "alg"],
            )
        assert len(function.free_mx()) == 0


def test_varlist_model_reusability():
    """Test if model created from supplied varlist and
    without varlist provides same results."""
    variable_list, model = mopeds.examples.pendulum_dae_1(False)
    time_grid = np.linspace(0, 1, 3)
    variable_list["g"].value = 12.0
    simulation = mopeds.Simulator(model, time_grid, variable_list)
    res_before = simulation.simulate_sym()

    variable_list, model = mopeds.examples.pendulum_dae_1(False, variable_list)
    simulation = mopeds.Simulator(model, time_grid, variable_list)
    res_after_pickle = simulation.simulate_sym()

    assert np.isclose(
        ca.vertcat(res_before["xf"], res_before["zf"]),
        ca.vertcat(res_after_pickle["xf"], res_after_pickle["zf"]),
    ).all()


if __name__ == "__main__":
    test_model(True, True, True, True)
    test_varlist_model_reusability()
