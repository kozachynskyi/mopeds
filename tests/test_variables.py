from conftest import VARIABLE_NAMES, VARIABLE_TYPES, generate_test_variables


def test_variables():
    variable_list = generate_test_variables()

    assert len(variable_list) == len(VARIABLE_TYPES)

    for var_index, var_name in enumerate(VARIABLE_NAMES):
        assert variable_list.get_casadi_var()[var_index].name() == var_name
