import pytest
from conftest import cstr_model_ode, pendulum_dae_1


def test_model():

    var_list, model = cstr_model_ode()

    assert len(model.varlist_state) == 5
    assert len(model.varlist_independent) == 18
    assert len(model.varlist_all) == 23

    assert model.equations_differential.size() == (5, 1)

    with pytest.raises(NotImplementedError):
        model.add_equations_differential(model.equations_differential[1])

    var_list, model = pendulum_dae_1()

    assert len(model.varlist_state) == 2
    assert len(model.varlist_independent) == 2
    assert len(model.varlist_algebraic) == 3
    assert len(model.varlist_all) == 7

    assert model.equations_differential.size() == (2, 1)
    assert model.equations_algebraic.size() == (3, 1)

    with pytest.raises(NotImplementedError):
        model.add_equations_algebraic(model.equations_algebraic[1])


if __name__ == "__main__":
    test_model()
