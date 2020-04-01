import pytest
from conftest import cstr_model_ode, pendulum_dae_1


def test_model():

    var_list, model = cstr_model_ode()

    assert len(model.states) == 5
    assert len(model.variables) == 18
    assert len(model._all_variables) == 23

    assert model.differential_equations.size() == (5, 1)

    with pytest.raises(NotImplementedError):
        model.add_differential_equations(model.differential_equations[1])

    var_list, model = pendulum_dae_1()

    assert len(model.states) == 2
    assert len(model.variables) == 2
    assert len(model.algebraic_variables) == 3
    assert len(model._all_variables) == 7

    assert model.differential_equations.size() == (2, 1)
    assert model.algebraic_equations.size() == (3, 1)

    with pytest.raises(NotImplementedError):
        model.add_algebraic_equations(model.algebraic_equations[1])


if __name__ == "__main__":
    test_model()
