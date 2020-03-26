import copy

import numpy as np

import par_est
from conftest import cstr_model_ode


# def test_simulation():
#     variable_list, m = cstr_model_ode()

#     # Create time-grid. Zero should be first
#     time_grid = np.linspace(10, 10000, 40)
#     time_grid = np.insert(time_grid, 0, 0)

#     # Generate experimental data
#     var_list_fixed = copy.deepcopy(variable_list)
#     for var in var_list_fixed.values():
#         var.fixed = True

#     sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
#     res = sim_fixed.generate_exp_data()
#     var_list_exp = sim_fixed.generate_exp_data()
