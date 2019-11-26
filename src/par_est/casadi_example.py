import matplotlib.pyplot as plt
import numpy as np
import casadi as ca

# Model Generation:

# Create variables
e0_T = ca.MX.sym("e0_T")
e0_c_i1 = ca.MX.sym("e0_c_i1")
e0_c_i2 = ca.MX.sym("e0_c_i2")
e0_c_i3 = ca.MX.sym("e0_c_i3")
e0_c_i4 = ca.MX.sym("e0_c_i4")

# Create controls
e0_c_in_i1 = ca.MX.sym("e0_c_in_i1")
e0_c_in_i2 = ca.MX.sym("e0_c_in_i2")
# e0_c_in_i3 = ca.MX.sym("e0_c_in_i3")
# e0_c_in_i4 = ca.MX.sym("e0_c_in_i4")
# e0_T_in = ca.MX.sym("e0_T_in")
# e0_T_j = ca.MX.sym("e0_T_j")

# Create parameters
e0_k_pre_r1 = ca.MX.sym("e0_k_pre_r1")  # 5000000.0
e0_k_pre_r2 = ca.MX.sym("e0_k_pre_r2")  # 1.0e7
e0_k_pre_r3 = ca.MX.sym("e0_k_pre_r3")  # 500000.0
e0_U = ca.MX.sym("e0_U")

# Parameter needed to have a flexible iteration step for Integrator
# If you use ca.Integrator than you have to specify a start_time and end_time for each instance.
# Meaning, that you when you evaluate an integrator you always have a fixed time step.
# The trick here is to create one instance of the integrator which integrate from zero to one. And if you
# multiply the right-han-side by tau parameter, you'll get a variable time-step. More info: https://groups.google.com/d/msg/casadi-users/iX66dHLUAqU/sK_STHgHBgAJ
tau = ca.MX.sym("tau")

# Create constants
e0_greek_nu_i1_r1 = -1.0
e0_greek_nu_i1_r2 = 1.0
e0_greek_nu_i2_r2 = -1.0
e0_greek_nu_i3_r1 = 1.0
e0_greek_nu_i1_r3 = -1.0
e0_greek_nu_i4_r3 = 1.0
e0_greek_Deltah_r1 = 0.0045
e0_greek_Deltah_r2 = -0.0055
e0_greek_Deltah_r3 = 0.0045
e0_greek_rho = 800.0
e0_A = 1.0
e0_E_r1 = 96000.0
e0_c_p = 3.5
e0_E_r2 = 72000.0
e0_E_r3 = 69000.0
e0_F = 6.5e-4
e0_R = 8.314
e0_V = 1.0
# e0_c_in_i1 = 5.0
# e0_c_in_i2 = 10.0
e0_c_in_i3 = 0.0
e0_c_in_i4 = 0.0
e0_T_in = 373.0
e0_T_j = 373.0

# Create equations
# fmt: off
tdot = (((((e0_F / e0_V) * ((e0_T_in - e0_T))) + (((e0_U * e0_A) / (e0_greek_rho * (e0_c_p * e0_V))) * ((e0_T_j - e0_T)))) + (((-e0_greek_Deltah_r1) / (e0_greek_rho * e0_c_p)) * (e0_k_pre_r1 * (e0_c_i1 * ca.exp(((-e0_E_r1) / (e0_R * e0_T))))))) + (((-e0_greek_Deltah_r2) / (e0_greek_rho * e0_c_p)) * (e0_k_pre_r2 * (e0_c_i2 * ca.exp(((-e0_E_r2) / (e0_R * e0_T))))))) + (((-e0_greek_Deltah_r3) / (e0_greek_rho * e0_c_p)) * (e0_k_pre_r3 * (e0_c_i1 * ca.exp(((-e0_E_r3) / (e0_R * e0_T))))))
c1dot = ((((e0_F / e0_V) * ((e0_c_in_i1 - e0_c_i1))) + (e0_greek_nu_i1_r1 * (e0_k_pre_r1 * (e0_c_i1 * ca.exp(((-e0_E_r1) / (e0_R * e0_T))))))) + (e0_greek_nu_i1_r2 * (e0_k_pre_r2 * (e0_c_i2 * ca.exp(((-e0_E_r2) / (e0_R * e0_T))))))) + (e0_greek_nu_i1_r3 * (e0_k_pre_r3 * (e0_c_i1 * ca.exp(((-e0_E_r3) / (e0_R * e0_T))))))
c2dot = ((e0_F / e0_V) * ((e0_c_in_i2 - e0_c_i2))) + (e0_greek_nu_i2_r2 * (e0_k_pre_r2 * (e0_c_i2 * ca.exp(((-e0_E_r2) / (e0_R * e0_T))))))
c3dot = ((e0_F / e0_V) * ((e0_c_in_i3 - e0_c_i3))) + (e0_greek_nu_i3_r1 * (e0_k_pre_r1 * (e0_c_i1 * ca.exp(((-e0_E_r1) / (e0_R * e0_T))))))
c4dot = ((e0_F / e0_V) * ((e0_c_in_i4 - e0_c_i4))) + (e0_greek_nu_i4_r3 * (e0_k_pre_r3 * (e0_c_i1 * ca.exp(((-e0_E_r3) / (e0_R * e0_T))))))
# fmt: on

rhs = ca.vertcat(tdot, c1dot, c2dot, c3dot, c4dot)
rhs_tau = ca.vertcat(tdot, c1dot, c2dot, c3dot, c4dot) * tau

controls = ca.vertcat(e0_c_in_i1, e0_c_in_i2)
states = ca.vertcat(e0_T, e0_c_i1, e0_c_i2, e0_c_i3, e0_c_i4)
parameters = ca.vertcat(e0_k_pre_r1, e0_k_pre_r2, e0_k_pre_r3, e0_U)

parameters_MAN = ca.vertcat(parameters, controls)
parameters_TAU = ca.vertcat(parameters_MAN, tau)

# Time grid is an array with timepoints at which I want to have variable values and their derivatives
time_grid = [0]
time_grid.extend(np.linspace(time_grid[-1] + 10, 1000, 2).tolist())
num_steps = len(time_grid) - 1

# Here a few iterators for testing and usage are intialized
ode_system_MAN = {"x": states, "p": parameters_MAN, "ode": rhs}
ode_system_TAU = {"x": states, "p": parameters_TAU, "ode": rhs_tau}

# This integrator returns values at each time_point, but you can get derivatives only at last time stamp
integrator_MAN = ca.integrator(
    "integrator",
    "cvodes",
    ode_system_MAN,
    {"grid": time_grid, "output_t0": False, "print_stats": False},
)

# This integrator is just for testing purposes
integrator_MAN_single = ca.integrator(
    "integrator",
    "cvodes",
    ode_system_MAN,
    {"tf": time_grid[-1], "output_t0": False, "print_stats": False},
)

# This integrator uses tau variable and is used in PE and OED
integrator_MAN_tau = ca.integrator(
    "integrator",
    "cvodes",
    ode_system_TAU,
    {
        "tf": 1,
        "output_t0": False,
        "print_stats": False,
        "linear_multistep_method": "adams",
    },
)

# Here initial values for integrators and optimization problems, as well
# as bound of optimization variables are defined. Scaling is needed for better behaviour of optimizer
states_init = ca.DM([273.0, 3.0, 10.0, 0.0, 0.0])
parameters_values = ca.DM([5000000.0, 1.0e7, 500000.0, 1.4])
controls_values = ca.DM([5.0, 10.0])

parameters_MAN_values = ca.vertcat(parameters_values, controls_values)

parameters_scale = parameters_values[0:4] * 0.1
controls_scale = controls_values * 0.1

parameters_guess = ca.DM([4000000.0, 1.0e6, 400000.0, 2.4]) / parameters_scale
parameters_lb = ca.DM([4000000.0, 1.0e6, 400000.0, 0.4]) / parameters_scale
parameters_ub = ca.DM([6000000.0, 1.0e8, 600000.0, 3.4]) / parameters_scale

controls_guess = controls_values
controls_lb = ca.DM([2.0, 8.0])
controls_ub = ca.DM([8.0, 12.0])


# Call integrator to generate experimental data for PE
res_integration = integrator_MAN(x0=states_init, p=parameters_MAN_values)
exp_data = res_integration["xf"]


# Formulation of PE Problem.
res_states_compare = []
x_init = states_init
prev_time_step = 0
# This loop runs iterator for given time_points
for time_step in time_grid[1:]:
    res_integration_tau_accum = integrator_MAN_tau(
        x0=x_init,
        p=ca.vertcat(
            parameters * parameters_scale, controls_values, time_step - prev_time_step
        ),
    )
    prev_time_step = time_step
    x_init = res_integration_tau_accum["xf"]
    if time_step == time_grid[1]:
        res_states_compare = res_integration_tau_accum["xf"]
    else:
        res_states_compare = ca.horzcat(
            res_states_compare, res_integration_tau_accum["xf"]
        )

# Here the objective function is calculated
error = exp_data - res_states_compare
objective = 0.5 * ca.dot(error, error)

# eval_error = ca.Function("eval_objective", [parameters_MAN], [error])
# eval_objective = ca.Function("eval_objective", [parameters_MAN], [objective])

# Settings for optimizer are set here
nlp = {"x": parameters, "f": objective}
nlp_solver = ca.nlpsol(
    "solver",
    "ipopt",
    nlp,
    {"verbose": False, "ipopt": {"derivative_test": "first-order"}},
)

res_solver = nlp_solver(x0=parameters_guess, lbx=parameters_lb, ubx=parameters_ub)
print(res_solver["x"])
print(res_solver["x"] * parameters_scale)

"""
ODE Routine is currently not yeilding great results
"""

# # ODE Routine
# res_states_ode = []
# x_init = states_init
# prev_time_step = 0
# for time_step in time_grid[1:]:
#     res_integration_tau_ode = integrator_MAN_tau(
#         x0=x_init, p=ca.vertcat(parameters, controls, time_step - prev_time_step),
#     )
#
#     integrator_jac = integrator_MAN_tau.factory("I_fwd", ["x0", "p"], ["jac:xf:p"])
#
#     res_integration_jac = integrator_jac(
#         x0=x_init, p=ca.vertcat(parameters, controls, time_step - prev_time_step),
#     )
#
#     prev_time_step = time_step
#     x_init = res_integration_tau_ode["xf"]
#     if time_step == time_grid[1]:
#         res_jac_ode = res_integration_jac["jac_xf_p"][:, 0:4]
#     else:
#         res_jac_ode = ca.vertcat(res_jac_ode, res_integration_jac["jac_xf_p"][:, 0:4])
#
# # res_jacobian = ca.jacobian(res_states_ode, ca.vertcat(parameters, controls))
# # res_jacobian = ca.jacobian(res_states_ode, parameters)
#
# # Check objective
# # eval_jacobian = ca.Function("eval_jacobian", [parameters, controls], [res_jac_ode])
# # sensitivity_matrix = eval_jacobian(parameters_values, controls)[:, 0:4]
#
# # # Calculate OBJ TRACE[FIM]
# sigma_diag = ca.DM([1, 1, 1, 1, 1]) * 1e20
# sigma_full = ca.diag(sigma_diag)
#
# measurement_matrix = ca.repmat(sigma_full, num_steps, num_steps)
# sensitivity_matrix = res_jac_ode
# fim_matrix = sensitivity_matrix.T @ measurement_matrix @ sensitivity_matrix
#
# eval_fim_matrix = ca.Function("eval_fim", [controls], [fim_matrix])
# fim_matrix_inv = ca.inv(fim_matrix)
# trace = ca.trace(fim_matrix_inv)
#
# eval_trace = ca.Function("eval_trace", [controls], [trace])
#
# # fim = sensitivity_matrix.T @ measurement_matrix @ sensitivity_matrix
# # trace = ca.trace(ca.inv(fim))
#
# nlp_ode = {"x": ca.vertcat(parameters, controls), "f": trace}
# nlp_solver_ode = ca.nlpsol(
#     "solver",
#     "ipopt",
#     nlp_ode,
#     #     {"verbose": False, "ipopt": {"hessian_approximation": "exact", "max_iter": 200,"derivative_test": 'first-order'}},
#     {
#         "verbose": True,
#         "ipopt": {
#             "hessian_approximation": "limited-memory",
#             "max_iter": 200,
#             "derivative_test": "first-order",
#         },
#     },
# )
#
# res_solver_ode = nlp_solver_ode(
#     x0=ca.vertcat(parameters_values, controls_guess),
#     lbx=ca.vertcat(parameters_values, controls_lb),
#     ubx=ca.vertcat(parameters_values, controls_ub),
# )
# # res_solver_ode = nlp_solver_ode(x0=controls_guess, lbx=controls_lb, ubx=controls_ub)
# print(res_solver_ode["x"])
