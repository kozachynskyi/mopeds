import numpy as np

import mopeds
import mopeds.examples

"""
Not sure if jacobian of a function that ipopt uses for solver direction
is actually derived correctly for OED function.

1. nlp_grad_f["f"] is NaN, not sure if it's supposed to be so
2. If you turn on verbose output of a simulator:
               "print_in": True,
               "print_out": True,
   you see that some functions, like asens1_asens5_integrator_tau, do not use correct vector for input ["p"] - sone of the values are zero. It can be connected to mapping of integrator in Simulator clas, which can be tested while using simple loop instead of mapaccum for "integrator".
"""

var_list, model = mopeds.examples.cstr_ode()
time_grid = np.linspace(10, 10000, 4)
time_grid = np.insert(time_grid, 0, 0)
for var in var_list.values():
    var.fixed = True

var_list["e0_E_r1"].fixed = False
# var_list["e0_E_r2"].fixed = False
# var_list["e0_E_r3"].fixed = False
# var_list["e0_k_pre_r1"].fixed = False
# var_list["e0_k_pre_r2"].fixed = False
# var_list["e0_k_pre_r3"].fixed = False

var_list["e0_c_in_i1"].fixed = False
# var_list["e0_c_in_i2"].fixed = False
# var_list["e0_c_in_i3"].fixed = False
# var_list["e0_c_in_i4"].fixed = False
# var_list["e0_T_in"].fixed = False
# var_list["e0_T_j"].fixed = False
# var_list["e0_F"].fixed = False

oed = mopeds.OptimalExperimentalDesign(model, [var_list], time_grid)
oed.solver_settings = {
    "verbose": True,
    "monitor": ["nlp_grad_f", "nlp_f"],
    "ipopt": {
        "hessian_approximation": "limited-memory",
        "max_iter": 1,
        # "print_level": 6,
    },
}

res = oed.optimize()
