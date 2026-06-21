import mopeds
import numpy as np
import matplotlib.pyplot as plt
import casadi as ca


def tank() -> [mopeds.VariableList, mopeds.Model]:
    variable_list = mopeds.VariableList()
    variable_list.add_variable(mopeds.VariableState("tau", 0, 0, 1200))
    variable_list.add_variable(mopeds.VariableState("T", 80, 0, 100))
    # variable_list.add_variable(mopeds.VariableState("Uint", 800.0, 0, 1200))
    variable_list.add_variable(mopeds.VariableState("I", 0, 0, 1200))
    variable_list.add_variable(mopeds.VariableAlgebraic("Qloss", 1.0, 0, 2))
    variable_list.add_variable(mopeds.VariableAlgebraic("Qheater", 16, 0, 2))
    variable_list.add_variable(mopeds.VariableAlgebraic("D", 0, 0, 100))
    variable_list.add_variable(
        mopeds.VariableControlPiecewiseConstant("Tsp", 70.0, 0, 100)
    )

    m = mopeds.Model(variable_list)

    tau = m.variables_all["tau"]
    Qloss = m.variables_all["Qloss"]
    # Uint = m.variables_all["Uint"]
    Qheater = m.variables_all["Qheater"]
    Tsp = m.variables_all["Tsp"]
    T = m.variables_all["T"]
    D = m.variables_all["D"]
    I = m.variables_all["I"]

    error = T - Tsp  # positive when no heating needed

    eq0 = 1
    dTdt = (-Qloss + Qheater) / 0.2
    eq1 = dTdt
    eq10 = error

    # eq2 = T - Uint/10
    eq3 = Qloss - T * 0.2

    eq4 = Qheater + 0.04 * (error * 1 + 10 * I + 10 * D)
    eq5 = D - dTdt

    m.add_equations_differential([eq0, eq1, eq10])
    m.add_equations_algebraic([eq3, eq4, eq5])

    variable_list["T"].ignore_plotting = False
    # variable_list["Qloss"].ignore_plotting = False
    variable_list["D"].ignore_plotting = False
    variable_list["Qheater"].ignore_plotting = False
    variable_list["Tsp"].ignore_plotting = False

    return variable_list, m


vl, m = tank()

time = np.linspace(0, 10, 2)
# vl["Tsp"].expand_horizon([5, 15], [80, 75])
# vl["Tsp"].expand_horizon([5, 15], [80, 75])
sim = mopeds.Simulator(m, time, vl)
res = sim.simulate(algebraic=True)[2]
res.plot(algebraic=True)
