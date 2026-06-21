import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples

plt.ion()

if __name__ == "__main__":
    piecewiseswitch = False
    variable_list, m = mopeds.examples.cstr_dae(piecewiseswitch)
    for var in variable_list.values():
        var.fixed = True

    variable_list["e0_U"].fixed = False
    variable_list["e0_c_p"].fixed = False
    variable_list["e0_E_r1"].fixed = False
    variable_list["e0_T_in"].fixed = True
    variable_list["e0_T"].variance = 0.1**2
    variable_list["e0_c_i1"].variance = 0.01**2
    variable_list["e0_c_i2"].variance = 0.01**2
    variable_list["e0_c_i3"].variance = 0.01**2
    variable_list["e0_c_i4"].variance = 0.01**2

    # Create time-grid. Zero should be first
    time_grid1 = np.linspace(0, 3000, 4)
    time_grid2 = np.linspace(0, 3000, 6)

    e0_T_in = variable_list["e0_T_in"]
    variable_list["e0_F"].fixed = False
    if isinstance(e0_T_in, mopeds.VariableControlPiecewiseConstant):
        e0_T_in.expand_horizon([10, 723], [363, 453])
        e0_T_in.variable_list.index(0).fixed = False
        e0_T_in.variable_list.index(1).fixed = True
        e0_T_in.variable_list.index(2).fixed = False

    data1 = mopeds.tools.generate_varlist_with_data(
        variable_list, m, time_grid1, True, perturbate=True
    )
    data2 = mopeds.tools.generate_varlist_with_data(
        variable_list, m, time_grid2, True, perturbate=True
    )

    pe = mopeds.ParameterEstimation(m, [data1, data2])
    # pe_state = mopeds.ParameterEstimation(m, [data1])
    # a = pe_state.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_c_p": 3.5, "e0_E_r1": 9.6e4})
    gcv_all = []
    # gcv_all = [2.3303451708286274e-11, 1.4993406699321207e-10, 2.817053414754338e-10, 4.0305248454108245e-10, 1.7860756306094184e-08, 4.399712463716451e-06, 0.0010992587399780297, 0.2553359281879561, 0.5549853507348389, np.nan, np.nan]

    for gamma in np.geomspace(1e-3, 1e3, 11):
        pe.setup_regularization(
            gamma, reference_parameters=np.zeros((len(pe.varlist_decision), 1))
        )
        pe.solver_settings["ipopt"]["max_iter"] = 50
        a = pe.optimize(True, objective_function="tikh")
        rrr = pe.calculate_sensitivity_and_fim(a["x_dict"])
        b = rrr["hess_wls"]
        g = rrr["hess_tikh"]
        neff = np.trace(b @ np.linalg.inv(g) @ b @ np.linalg.inv(g))
        # print(np.trace(rrr["hess_wls"]))
        # print(np.trace(rrr["hess_tikh"]))
        # neff = np.trace(pe.calculate_sensitivity_and_fim(a["x_dict"])["hess_tikh"])
        obj = pe.calculate_objective_and_residual(
            a["x_dict"], objective_function="wls"
        )["f"]
        gcv = obj / (pe.dof - neff)
        gcv_all.append(gcv)
        # print(neff)
        print(gcv)

    plt.scatter(np.geomspace(1e-3, 1e3, 11), gcv_all)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

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

    # oed = mopeds.OptimalExperimentalDesign(m, [data1], time_grid1, time_grid1, measurable_variables=mes_names, simulator_name="idas")
    oed.guess[0] = 373
    # print(oed.calculate_objective_and_jacobian({"e0_T_in": 373})["jac"])

    print(oed.optimize(objective_function="A_fd"))
    print(oed.optimize())
    breakpoint()
    # oed.optimize()
