import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples

plt.ion()

if __name__ == "__main__":
    piecewiseswitch = False
    variable_list, m = mopeds.examples.cstr()
    for var in variable_list.values():
        var.fixed = True

    variable_list["e0_U"].fixed = False
    variable_list["e0_c_p"].fixed = True
    # variable_list["e0_E_r1"].fixed = False
    variable_list["e0_T_in"].fixed = False
    variable_list["e0_T"].variance = 1

    # Create time-grid. Zero should be first
    time_grid1 = np.linspace(0, 1000, 4)
    time_grid2 = np.linspace(0, 1000, 8)

    e0_T_in = variable_list["e0_T_in"]
    variable_list["e0_F"].fixed = False
    if isinstance(e0_T_in, mopeds.VariableControlPiecewiseConstant):
        e0_T_in.expand_horizon([10, 723], [363, 453])
        e0_T_in.variable_list.index(0).fixed = False
        e0_T_in.variable_list.index(1).fixed = True
        e0_T_in.variable_list.index(2).fixed = False

    data1 = mopeds.tools.generate_varlist_with_data(variable_list, m, time_grid1, True)
    data2 = mopeds.tools.generate_varlist_with_data(variable_list, m, time_grid2, True)
    # data1.show()

    # If data is not available for all simulated points, PE works
    e0_T = data2["e0_T"]
    e0_T.dataframe = e0_T._dataframe_from_value(e0_T.value[0])

    e0_c_i1_df = data2["e0_c_i1"].dataframe
    e0_c_i1_df.drop(e0_c_i1_df.index[2:], inplace=True)

    # Perturbate alg variable from "ideal solution" to see its affect on PE
    a = data2["e0_c_tot"].dataframe
    data2["e0_c_tot"].dataframe = data2["e0_c_tot"].dataframe * 1.05

    pe_state = mopeds.ParameterEstimation(m, [data1, data2])
    # pe_state = mopeds.ParameterEstimation(m, [data1])
    v = pe_state.calculate_sensitivity_and_fim(
        {"e0_U": 1.4, "e0_c_p": 3.5, "e0_E_r1": 9.6e4}
    )
    print(a)
    # print(pe_state.optimize(True))

    # pe_alg = mopeds.ParameterEstimation(m, [data1, data2], use_algebraic_vars=True)
    # print(pe_alg.optimize(True))

    data1["e0_T"].fixed = False
    data1["e0_U"].fixed = False
    data1["e0_c_p"].fixed = True
    data1["e0_E_r1"].fixed = True
    data1["e0_F"].fixed = True
    data1["e0_c_i1"].fixed = False
    # mes_names = ["e0_T", "e0_c_i1"]
    mes_names = ["e0_c_i1", "e0_T"]
    # mes_names = ["e0_T"]
    # mes_names = None

    oed = mopeds.OptimalExperimentalDesign(
        m,
        [data1],
        time_grid1,
        time_grid1,
        measurable_variables=mes_names,
        simulator_name="idas",
    )
    oed.guess[0] = 373
    # print(oed.calculate_objective_and_jacobian({"e0_T_in": 373})["jac"])

    print(oed.optimize(objective_function="A_fd"))
    print(oed.optimize())
    breakpoint()
    # oed.optimize()
