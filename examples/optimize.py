import numpy as np

import par_est
import par_est.examples

if __name__ == "__main__":

    variable_list, m = par_est.examples.cstr_ode()
    for var in variable_list.values():
        var.fixed = True

    variable_list["e0_U"].fixed = False
    variable_list["e0_T_j"].fixed = False
    variable_list["e0_T"].variance = 0.1

    # Create time-grid. Zero should be first
    time_grid1 = np.linspace(0, 1000, 4)
    time_grid2 = np.linspace(0, 1000, 8)

    data1 = par_est.tools.generate_exp_data(variable_list, m, time_grid1)
    data2 = par_est.tools.generate_exp_data(variable_list, m, time_grid2)

    # If data is not available for all simulated points, PE works
    data2["e0_T"].value.value = np.delete(data2["e0_T"].value.value, 2)
    data2["e0_T"].value.time = np.delete(data2["e0_T"].value.time, 2)

    data2["e0_c_i1"].value.value = [data2["e0_c_i1"].value.value[0]]
    data2["e0_c_i1"].value.time = [data2["e0_c_i1"].value.time[0]]

    pe = par_est.ParameterEstimation(m, [data1, data2])
    pe.optimize(scale_experiments=True)
    # pe.optimize(False)

    oed = par_est.OptimalExperimentalDesign(m, [data1], time_grid1)
    # oed.optimize()
    # oed.optimize(False)
