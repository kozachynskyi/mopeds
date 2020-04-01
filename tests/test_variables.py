import par_est


def test_variables():
    variable_list = par_est.VariableList()

    for variable_type in par_est.Variable.get_subclasses():
        var = variable_type("Name")
        variable_list.add_variable(var)
        assert var.casadi_var.name() == "Name"

    assert len(variable_list) == 1

    variable_list = par_est.VariableList()
    counter = 0

    for variable_type in par_est.Variable.get_subclasses():
        var = variable_type(f"Name{counter}")
        variable_list.add_variable(var)
        counter += 1

    assert len(variable_list) == sum(1 for _ in par_est.Variable.get_subclasses())


if __name__ == "__main__":
    test_variables()
