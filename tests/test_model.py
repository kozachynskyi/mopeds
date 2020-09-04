import pytest
import casadi as ca
import par_est.examples


def test_model():

    var_list, model = par_est.examples.cstr_ode()

    assert len(model.varlist_state) == 5
    assert len(model.varlist_independent) == 18
    assert len(model.varlist_all) == 23

    assert model.equations_differential.size() == (5, 1)

    with pytest.raises(NotImplementedError):
        model.add_equations_differential(model.equations_differential[1])

    var_list, model = par_est.examples.pendulum_dae_1()

    assert len(model.varlist_state) == 2
    assert len(model.varlist_independent) == 2
    assert len(model.varlist_algebraic) == 3
    assert len(model.varlist_all) == 7

    assert model.equations_differential.size() == (2, 1)
    assert model.equations_algebraic.size() == (3, 1)

    with pytest.raises(NotImplementedError):
        model.add_equations_algebraic(model.equations_algebraic[1])

    for model in [par_est.examples.cstr_dae_constant, par_est.examples.cstr_ode_constant]:
        var_list, model = model()
        ode_system = {
            "x": model.varlist_state.get_casadi_var(),
            "p": ca.vertcat(model.varlist_independent.get_casadi_var()),
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
            ode_system["z"] = model.varlist_algebraic.get_casadi_var()
            ode_system["alg"] = model.equations_algebraic

            function = ca.Function(
                "eq_sys",
                [ode_system["x"], ode_system["z"], ode_system["p"]],
                [ode_system["ode"], ode_system["alg"]],
                ["x", "z", "p"],
                ["ode", "alg"],
            )
        assert len(function.free_mx()) == 0


if __name__ == "__main__":
    test_model()
