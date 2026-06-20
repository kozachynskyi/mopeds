import copy

import mopeds


def initialize_problem():  # noqa: C901

    variable_list = mopeds.VariableList()
    # fmt:off

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_F_s2", 20.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_F_s4", 7.0, -1.0E9, 1.0E9))  # noqa: E501

    variable_list.add_variable(mopeds.VariableControl("e0_F_s1", 21.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_F_s3", 13.0, -1.0E9, 1.0E9))  # noqa: E501

    m = mopeds.Model(variable_list)

    e0_F_s1 = m.varlist_all["e0_F_s1"].casadi_var  # noqa: E501
    e0_F_s3 = m.varlist_all["e0_F_s3"].casadi_var  # noqa: E501
    e0_F_s2 = m.varlist_all["e0_F_s2"].casadi_var  # noqa: E501
    e0_F_s4 = m.varlist_all["e0_F_s4"].casadi_var  # noqa: E501

    EQ_alg1 = (0.0-((e0_F_s1-e0_F_s2)))  # noqa: E501,E226
    EQ_alg2 = (0.0-(((e0_F_s2-e0_F_s3)-e0_F_s4)))  # noqa: E501,E226

    list_algebraic_equations = [EQ_alg1, EQ_alg2, ]  # noqa: E501

    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


if __name__ == "__main__":

    variable_list, m = initialize_problem()

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = mopeds.SimulatorNLE(m, var_list_fixed)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate_fast()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.simulate()[2]
    print(res_simple)
    print(res.dataframe)
