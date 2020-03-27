import pytest
import casadi as ca
from conftest import VARIABLE_NAMES, generate_test_model


def test_model():

    model = generate_test_model()

    assert len(model.states) == 1
    assert len(model.variables) == len(VARIABLE_NAMES) - 1
    assert len(model._all_variables) == len(VARIABLE_NAMES)

    assert model.differential_equations.size() == (1, 1)

    equation = ca.MX.sym("test")
    with pytest.raises(NotImplementedError):
        model.add_differential_equations([equation])


test_model()
