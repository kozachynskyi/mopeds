import copy

import dae_ode.hyfo_dae  # type: ignore
import numpy as np

import mopeds
import mopeds.examples

if __name__ == "__main__":

    variable_list, m = dae_ode.hyfo_dae.initialize_problem()

    variable_list.set_bounds(0.01)

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 1000, 8)
    time_grid = np.insert(time_grid, 0, 0)

    variable_list_optimizer = copy.deepcopy(variable_list)

    variable_list_optimizer["e0_K_cat_e1"].fixed = False
    variable_list_optimizer["e0_K_cat_e2"].fixed = False
    variable_list_optimizer["e0_greek_DeltaG_r1"].fixed = False
    variable_list_optimizer["e0_E_r1"].fixed = False
    variable_list_optimizer["e0_K_r1_e1"].fixed = False
    variable_list_optimizer["e0_K_r1_e2"].fixed = False
    variable_list_optimizer["e0_K_LM"].fixed = False
    # variable_list_optimizer["e0_P_Surfactant"].fixed = False
    variable_list_optimizer["e0_P_trig_r1"].fixed = False
    variable_list_optimizer["e0_k_LM_r1"].fixed = False
    variable_list_optimizer["e0_k_ref_r1"].fixed = False
    variable_list_optimizer["e0_E_r2"].fixed = False
    variable_list_optimizer["e0_E_r3"].fixed = False
    variable_list_optimizer["e0_K_r3_e1"].fixed = False
    variable_list_optimizer["e0_K_r3_e2"].fixed = False
    variable_list_optimizer["e0_K_r3_e3"].fixed = False
    variable_list_optimizer["e0_E_r4"].fixed = False
    variable_list_optimizer["e0_P_trig_Hyfo"].fixed = False
    variable_list_optimizer["e0_k_ref_r4"].fixed = False
    variable_list_optimizer["e0_k_LM_Hyfo"].fixed = False
    variable_list_optimizer["e0_E_r5"].fixed = False
    variable_list_optimizer["e0_K_r5_e1"].fixed = False
    variable_list_optimizer["e0_K_r5_e2"].fixed = False
    variable_list_optimizer["e0_K_r5_e3"].fixed = False
    variable_list_optimizer["e0_E_r6"].fixed = False
    variable_list_optimizer["e0_k_ref_r6"].fixed = False

    variable_list_optimizer["e0_T"].fixed = False
    # variable_list_optimizer["e0_p_Reactor"].fixed = False

    oed = mopeds.OptimalExperimentalDesign(m, [variable_list_optimizer], time_grid)

    unfix_parameters = oed.identifiability_analysis()
    print(unfix_parameters)

    oed.optimize()
