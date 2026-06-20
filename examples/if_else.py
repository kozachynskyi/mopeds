import casadi as ca
import mopeds


def example() -> [mopeds.VariableList, mopeds.Model]:
    variable_list = mopeds.VariableList()
    variable_list.add_variable(mopeds.VariableAlgebraic("y", 1.0))
    variable_list.add_variable(mopeds.VariableControl("u", 1.0, 1, 2))
    variable_list.add_variable(mopeds.VariableParameter("a", 1.0))
    variable_list.add_variable(mopeds.VariableParameter("b", 2.0))

    m = mopeds.Model(variable_list)

    y = m.variables_all["y"]
    u = m.variables_all["u"]
    a = m.variables_all["a"]
    b = m.variables_all["b"]

    condition = a + b * u
    eq1 = y - ca.if_else(condition<0, 0, condition)

    m.add_equations_algebraic([eq1])

    return variable_list, m

vl, m = example()
sim = mopeds.SimulatorNLE(m, vl)
print(sim.simulate())

sim.change_independent_variables({"u": -10})
print(sim.simulate())

