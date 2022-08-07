import copy

import casadi as ca
import matplotlib.pyplot as plt
import numpy as np

import par_est


def get_model(linear=True):
    varlist = par_est.VariableList()
    varlist.add_variable(par_est.VariableAlgebraic("y", 1))
    varlist.add_variable(par_est.VariableControl("x", 1))
    varlist.add_variable(par_est.VariableParameter("theta1", 1))
    varlist.add_variable(par_est.VariableParameter("theta2", 2))

    model = par_est.Model(varlist)

    y = model.varlist_all["y"].casadi_var
    x = model.varlist_all["x"].casadi_var
    theta1 = model.varlist_all["theta1"].casadi_var
    theta2 = model.varlist_all["theta2"].casadi_var

    if linear:
        eq1 = y - (theta1 + theta2 * x)
    else:
        eq1 = y - (theta1 * ca.log(theta2 / x))

    model.add_equations_algebraic([eq1])

    return model, varlist


def generate_data(
    linear,
    x_bounds: list[float],
    preturbate=True,
    num_x: int = 20,
    parameters: dict = None,
):
    model, varlist = get_model(linear)
    if parameters is not None:
        for par, val in parameters.items():
            varlist[par].value = val
    simulator = par_est.SimulatorNLE(model, varlist)
    x_range = np.linspace(x_bounds[0], x_bounds[1], num_x)
    y_range = []
    rng = np.random.default_rng(12345)
    # rng = np.random.default_rng()

    for x in x_range:
        simulator._independent_variables[0] = x
        res = simulator.simulate_sym()
        y = float(res["x"])

        if preturbate:
            y = rng.normal(y, 0.05)
        y_range.append(y)

    y_range = np.array(y_range)
    return x_range, y_range


def generate_varlist(linear, x_range, y_range):
    model, varlist = get_model(linear)
    varlist["theta1"].fixed = False
    varlist["theta2"].fixed = False
    varlist_list = []

    for x, y in zip(x_range, y_range):
        varlist_i = copy.deepcopy(varlist)
        varlist_i["x"].value = x
        varlist_i["y"].value = y
        varlist_list.append(varlist_i)

    return varlist_list


if __name__ == "__main__":
    LINEAR = True
    TRUE_PARAMETERS = {"theta1": 1, "theta2": 2}
    PARAMETERS_2 = {"theta1": 1, "theta2": -2}
    model, varlist = get_model(LINEAR)
    BOUNDS = [0.5, 1.5]
    BOUNDS2 = [3.5, 4.5]

    for parameters in [TRUE_PARAMETERS]:  # , PARAMETERS_2]:
        for bounds in [BOUNDS, BOUNDS2]:
            x, y = generate_data(LINEAR, bounds, False, 30, parameters)
            x_data, y_data = generate_data(LINEAR, bounds, True, 5, parameters)
            plt.errorbar(x_data, y_data, 0.3, linestyle="None", marker=".")
            plt.plot(x, y)
            varlist_list = generate_varlist(LINEAR, x_data, y_data)

    for bounds in [BOUNDS, BOUNDS2]:
        x_data, y_data = generate_data(LINEAR, bounds, True, 20, parameters)
        varlist_list = generate_varlist(LINEAR, x_data, y_data)
        pe = par_est.ParameterEstimationNLE(model, varlist_list)
        pe.parameter_analysis(TRUE_PARAMETERS)
        plt.title(str(bounds))

    plt.show()
    for parameters in [TRUE_PARAMETERS, PARAMETERS_2]:
        theta1 = []
        theta2 = []
        for i in range(1):
            x_data, y_data = generate_data(LINEAR, BOUNDS, True, 20, parameters)
            varlist_list = generate_varlist(LINEAR, x_data, y_data)
            pe = par_est.ParameterEstimationNLE(model, varlist_list)
            pe.guess = np.array(list(parameters.values()))
            if i == 0:
                pe.parameter_analysis(parameters)
            res = pe.optimize()
            theta1.append(float(res["x"][0]))
            theta2.append(float(res["x"][1]))

        plt.scatter(theta1, theta2)

    plt.show()
    breakpoint()
