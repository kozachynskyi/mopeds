import numpy as np
from matplotlib import pyplot as plt

import par_est
import par_est.examples

plt.ion()

if __name__ == "__main__":

    piecewiseswitch = False
    variable_list, m = par_est.examples.cstr_dae(piecewiseswitch)
    res_dict = {}
    for var in variable_list.values():
        var.fixed = True
        if isinstance(var, par_est.VariableParameter):
            par_dict[var.name] = var.value[0]

    variable_list["e0_U"].fixed = False
    variable_list["e0_E_r1"].fixed = False
    variable_list["e0_E_r2"].fixed = False
    variable_list["e0_E_r3"].fixed = False
    variable_list["e0_T_in"].fixed = False
    variable_list["e0_T"].variance = 0.1

    # Create time-grid. Zero should be first
    time_grid1 = np.linspace(0, 1000, 2)
    time_grid2 = np.linspace(0, 1000, 4)

    e0_T_in = variable_list["e0_T_in"]
    if isinstance(e0_T_in, par_est.VariableControlPiecewiseConstant):
        e0_T_in.expand_horizon([10, 723], [363, 453])
        e0_T_in.variable_list.index(0).fixed = False
        e0_T_in.variable_list.index(1).fixed = True
        e0_T_in.variable_list.index(2).fixed = False

    data1 = par_est.tools.generate_varlist_with_data(variable_list, m, time_grid1, True)
    data2 = par_est.tools.generate_varlist_with_data(variable_list, m, time_grid2, True)
    # data1.show()

    # If data is not available for all simulated points, PE works
    for d in [data1, data2]:
        e0_T = d["e0_T"]
        e0_T.dataframe = e0_T._dataframe_from_value(e0_T.value[0])

    for i in range(1,2):
        e0_c_i1_df = data2[f"e0_c_i{i}"].dataframe
        e0_c_i1_df.drop(e0_c_i1_df.index[2:], inplace=True)

    # Perturbate alg variable from "ideal solution" to see its affect on PE
    a = data2["e0_c_tot"].dataframe
    data2["e0_c_tot"].dataframe = data2["e0_c_tot"].dataframe * 1.05

    pe = par_est.ParameterEstimation(m, [data1, data2])
    a = pe.calculate_objective_and_residual(res_dict, "wls")

    res = pe.optimize()

    breakpoint()
    # pe_state = par_est.ParameterEstimation(m, [data2])
    # print(pe_state.optimize(True))

    # pe_alg = par_est.ParameterEstimation(m, [data1, data2], use_algebraic_vars=True)
    # print(pe_alg.optimize(True))

    # oed = par_est.OptimalExperimentalDesign(m, [data1], time_grid1)
    # oed.optimize()
