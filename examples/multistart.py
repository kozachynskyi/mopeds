import copy
import numpy as np

import par_est
import par_est.examples

if __name__ == "__main__":

    variable_list, m = par_est.examples.cstr_ode()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 1000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    # Generate experimental data for parameter estimation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).generate_exp_data()

    # Replace empty state variables with results from simulation
    variable_list_optimizer = copy.deepcopy(variable_list)
    for key, var in var_list_exp.items():
        variable_list_optimizer[key] = var

    # variable_list_optimizer["e0_E_r1"].fixed = True
    # variable_list_optimizer["e0_E_r2"].fixed = True
    # variable_list_optimizer["e0_E_r3"].fixed = True
    # variable_list_optimizer["e0_k_pre_r1"].fixed = True
    # variable_list_optimizer["e0_k_pre_r2"].fixed = True
    # variable_list_optimizer["e0_k_pre_r3"].fixed = True
    # variable_list_optimizer["e0_U"].fixed = True
    # variable_list_optimizer["e0_c_p"].fixed = True
    # variable_list_optimizer["e0_greek_Deltah_r1"].fixed = True
    # variable_list_optimizer["e0_greek_Deltah_r2"].fixed = True
    # variable_list_optimizer["e0_greek_Deltah_r3"].fixed = True

    # variable_list_optimizer["e0_c_in_i1"].fixed = True
    # variable_list_optimizer["e0_c_in_i2"].fixed = True
    # variable_list_optimizer["e0_c_in_i3"].fixed = True
    # variable_list_optimizer["e0_c_in_i4"].fixed = True
    variable_list_optimizer["e0_T_in"].fixed = True
    variable_list_optimizer["e0_T_j"].fixed = True
    # variable_list_optimizer["e0_F"].fixed = True

    pe = par_est.ParameterEstimation(
        m, [variable_list_optimizer, variable_list_optimizer]
    )

    res1 = pe.optimize_multistart(3, True, 10)

    oed = par_est.OptimalExperimentalDesign(m, [variable_list_optimizer], time_grid)
    print(oed.identifiability_analysis(True))
    res2 = oed.optimize_multistart(3, False, 2)
