import par_est

y = par_est.VariableAlgebraic("y", 8.3)
x = par_est.VariableControl("x", 1)
C = par_est.VariableConstant("C", 1)
theta1 = par_est.VariableParameter("theta1", 20)
theta2 = par_est.VariableParameter("theta2", 0.24)

variable_list = par_est.VariableList()
variable_list.add_variable(y)
variable_list.add_variable(x)
variable_list.add_variable(C)
variable_list.add_variable(theta1)
variable_list.add_variable(theta2)

model = par_est.Model(variable_list)

y = model.varlist_all["y"].casadi_var
x = model.varlist_all["x"].casadi_var
C = model.varlist_all["C"].casadi_var
theta1 = model.varlist_all["theta1"].casadi_var
theta2 = model.varlist_all["theta2"].casadi_var

import casadi as ca

equation = y - (theta1 * (C - ca.exp(-theta2 * x)))

model.add_equations_algebraic([equation])

simulator = par_est.SimulatorNLE(model, variable_list)
result = simulator.generate_exp_data()

print(result["y"].value)
# >>> [4.267442778668931]
print(result["y"].dataframe)
# >>>                   y
# >>>1970-01-01  4.267443

variable_list["x"].value = 20
simulator = par_est.SimulatorNLE(model, variable_list)
result = simulator.generate_exp_data()

print(result["y"].value)
# >>> [19.8354050590196]

import copy

var_list_1 = copy.deepcopy(variable_list)
var_list_1["x"].value = 1
var_list_1["y"].value = 8.3

var_list_2 = copy.deepcopy(variable_list)
var_list_2["x"].value = 7
var_list_2["y"].value = 19.8

experimental_data = [var_list_1, var_list_2]

for var_list in experimental_data:
    var_list["theta1"].fixed = False
    var_list["theta1"].guess = 40
    var_list["theta1"].lower_bound = 0
    var_list["theta1"].upper_bound = 40

pe = par_est.ParameterEstimationNLE(model, experimental_data)
result = pe.optimize()
print(result["x"])
# >>> 25.2727

simulator = par_est.SimulatorNLE(model, variable_list)

solver_settings = {
    "nlpsol": "ipopt",
    "verbose": False,
    "print_in": False,
    "print_out": False,
    "expand": True,
    "nlpsol_options": {
        "ipopt.hessian_approximation": "limited-memory",
        "ipopt.max_iter": 300,
        "ipopt.print_level": 0,
        "print_time": False,
    },
}

simulator = par_est.SimulatorNLE(
    model,
    variable_list,
    solver_settings=solver_settings,
    solver_name="rootfinder",
    use_bounds=True,
)

solver_settings = {
    "nlpsol": "ipopt",
    "verbose": False,
    "print_in": False,
    "print_out": False,
    "expand": True,
    "nlpsol_options": {
        "ipopt.hessian_approximation": "exact",
        "ipopt.linear_solver": "ma57",
        "ipopt.ma57_automatic_scaling": "yes",
        "ipopt.max_iter": 300,
        "ipopt.print_level": 0,
        "print_time": False,
    },
}

simulator = par_est.SimulatorNLE(
    model,
    variable_list,
    solver_settings=solver_settings,
    solver_name="rootfinder",
    use_bounds=True,
)


solver_settings = {
    "nlpsol": "qrsqp",
    "verbose": False,
    "print_in": False,
    "print_out": False,
    "expand": True,
    "nlpsol_options": {
        "print_iteration": False,
    },
}

simulator = par_est.SimulatorNLE(
    model,
    variable_list,
    solver_settings=solver_settings,
    solver_name="rootfinder",
    use_bounds=True,
)


solver_settings = {
    "nlpsol": "ipopt",
    "verbose": False,
    "print_in": False,
    "print_out": False,
    "expand": True,
    "nlpsol_options": {
        "ipopt.hessian_approximation": "limited-memory",
        "ipopt.max_iter": 300,
        "ipopt.print_level": 0,
        "print_time": False,
    },
}

pe = par_est.ParameterEstimationNLE(
    model,
    experimental_data,
    simulator_settings=solver_settings,
    simulator_name="rootfinder",
    use_simulator_bounds=True,
)

pe.solver_name = "ipopt"
pe.solver_settings = {
    "verbose": False,
    "ipopt": {"max_iter": 300},
}

pe = par_est.ParameterEstimationNLE(model, experimental_data)
pe.solver_name = "qrsqp"
pe.solver_settings = {
    "verbose": False,
    "max_iter": 300,
}

pe.optimize()
