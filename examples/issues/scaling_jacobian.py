import copy

import casadi as ca
import numpy as np

import par_est
import par_est.examples
import par_est.tools

var_list, m = par_est.examples.cstr_dae()
# Create time-grid. Zero should be first
time_grid = np.linspace(10, 10000, 4)
time_grid = np.insert(time_grid, 0, 0)

var_list_fixed = copy.deepcopy(var_list)
for var in var_list_fixed.values():
    var.fixed = True
var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).generate_exp_data()

for key, var in var_list_exp.items():
    var_list[key] = var

for var in var_list.values():
    var.fixed = True

var_list["e0_E_r1"].fixed = False

pe = par_est.ParameterEstimation(m, [var_list])

pe._setup_scaling(False)
sim_unscaled = pe.list_simulators[0].simulate_jac()
ev_unscaled = ca.Function(
    "unscaled",
    [pe.varlist_decision.get_casadi_variables()],
    [sim_unscaled["xf"], sim_unscaled["jac_xf_p"]],
)

result_unscaled = ev_unscaled(90000)

pe._setup_scaling(True)
sim_scaled = pe.list_simulators[0].simulate_jac()
ev_scaled = ca.Function(
    "scaled",
    [pe.varlist_decision.get_casadi_variables()],
    [sim_scaled["xf"], sim_scaled["jac_xf_p"]],
)
results_scaled = ev_scaled(1)


# Simulation results are the same
assert (results_scaled[0] - result_unscaled[0]).is_zero()

difference_jacobian = results_scaled[1] - result_unscaled[1]
# Jacobian at all timepoints > 1 is the same
assert difference_jacobian[:, 19:38].is_zero()
assert difference_jacobian[:, 38:57].is_zero()
assert difference_jacobian[:, 57:76].is_zero()

# Jacobian at first time point is off by small values.
print(difference_jacobian[:, 0:19])
