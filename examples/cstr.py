import copy

import casadi as ca
import matplotlib.pyplot as plt

import mopeds
import numpy as np


def initialize_problem(mode="old"):  # noqa: C901

    variable_list = mopeds.VariableList()

    # fmt:off
    def fun_175223__arrhenius(std_T,std_E_A,std_R,std_k_pre):  # noqa: E501,E231,E306
        std_k = (std_k_pre*ca.exp((-(std_E_A/(std_R*std_T)))))  # noqa: E501,E226
        return std_k
    def fun_175242__activity_proton(std_T,std_x_A,std_x_B,std_x_cat,std_Param_Hplus_A,std_Param_Hplus_B,std_Param_Hplus_C):  # noqa: E501,E231,E306
        std_a_Hplus = (((1.0-std_x_A))*(((std_x_cat+(std_Param_Hplus_A*std_x_B))))**(1.0*(std_Param_Hplus_B*ca.exp((std_Param_Hplus_C/std_T)))))  # noqa: E501,E226
        return std_a_Hplus

    if mode == "gc":
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i1_gc", 0.186, -10000.0, 1.0E9))  # noqa: E501
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i2_gc", 0.616, -1.0E9, 1.0E9))  # noqa: E501
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i3_gc", 0.004, -1.0E9, 1.0E9))  # noqa: E501
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i4_gc", 0.004, -1.0E9, 1.0E9))  # noqa: E501

    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_i1", -1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_i2", -1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_i3", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_i4", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_R", 8.314))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_K_reac", 30.313046554204817, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_a_Hplus_l1", 0.035817754492527894, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_for", 2.588504783617644, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_i4_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i1_l1", 0.186, -10000.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i2_l1", 0.616, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i3_l1", 0.004, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i4_l1", 0.194, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_i1_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_i2_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_i3_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_F_i2", 0.8, -1.0E9, 1.0E9))  # noqa: E501



    variable_list.add_variable(mopeds.VariableControl("e0_greek_gamma_i1_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_gamma_i2_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    # variable_list.add_variable(mopeds.VariableControl("e0_F_F", 1, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_F_F", 2.8e-2, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_gamma_i3_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    # variable_list.add_variable(mopeds.VariableControl("e0_F_I", 1, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_F_I", 2.8e-2, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_F_i1", 0.2, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_F_i3", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_F_i4", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_T", 353.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_gamma_i4_l1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_cat", 0.11, -1.0E9, 1.0E9))  # noqa: E501

    variable_list.add_variable(mopeds.VariableParameter("e0_E_for_A", 48576.2078, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_reac_A", -6948.0098, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_Param_Hplus_A", 8.49E-4, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_Param_Hplus_B", 1.016, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_Param_Hplus_C", 116.6, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_k_for_pre", 3.993E7, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_k_reac_pre", 2.841, -1.0E9, 1.0E9))  # noqa: E501

    m = mopeds.Model(variable_list)

    e0_greek_gamma_i1_l1 = m.varlist_all["e0_greek_gamma_i1_l1"].casadi_var  # noqa: E501
    e0_greek_gamma_i2_l1 = m.varlist_all["e0_greek_gamma_i2_l1"].casadi_var  # noqa: E501
    e0_F_F = m.varlist_all["e0_F_F"].casadi_var  # noqa: E501
    e0_greek_gamma_i3_l1 = m.varlist_all["e0_greek_gamma_i3_l1"].casadi_var  # noqa: E501
    e0_F_I = m.varlist_all["e0_F_I"].casadi_var  # noqa: E501
    e0_x_F_i1 = m.varlist_all["e0_x_F_i1"].casadi_var  # noqa: E501
    e0_x_F_i2 = m.varlist_all["e0_x_F_i2"].casadi_var  # noqa: E501
    e0_x_F_i3 = m.varlist_all["e0_x_F_i3"].casadi_var  # noqa: E501
    e0_x_F_i4 = m.varlist_all["e0_x_F_i4"].casadi_var  # noqa: E501
    e0_T = m.varlist_all["e0_T"].casadi_var  # noqa: E501
    e0_greek_gamma_i4_l1 = m.varlist_all["e0_greek_gamma_i4_l1"].casadi_var  # noqa: E501
    e0_x_cat = m.varlist_all["e0_x_cat"].casadi_var  # noqa: E501
    e0_K_reac = m.varlist_all["e0_K_reac"].casadi_var  # noqa: E501
    e0_a_Hplus_l1 = m.varlist_all["e0_a_Hplus_l1"].casadi_var  # noqa: E501
    e0_k_for = m.varlist_all["e0_k_for"].casadi_var  # noqa: E501
    e0_greek_nu_i1 = m.varlist_all["e0_greek_nu_i1"].casadi_var  # noqa: E501
    e0_greek_nu_i2 = m.varlist_all["e0_greek_nu_i2"].casadi_var  # noqa: E501
    e0_greek_nu_i3 = m.varlist_all["e0_greek_nu_i3"].casadi_var  # noqa: E501
    e0_greek_nu_i4 = m.varlist_all["e0_greek_nu_i4"].casadi_var  # noqa: E501
    e0_E_for_A = m.varlist_all["e0_E_for_A"].casadi_var  # noqa: E501
    e0_E_reac_A = m.varlist_all["e0_E_reac_A"].casadi_var  # noqa: E501
    e0_Param_Hplus_A = m.varlist_all["e0_Param_Hplus_A"].casadi_var  # noqa: E501
    e0_Param_Hplus_B = m.varlist_all["e0_Param_Hplus_B"].casadi_var  # noqa: E501
    e0_Param_Hplus_C = m.varlist_all["e0_Param_Hplus_C"].casadi_var  # noqa: E501
    e0_R = m.varlist_all["e0_R"].casadi_var  # noqa: E501
    e0_k_for_pre = m.varlist_all["e0_k_for_pre"].casadi_var  # noqa: E501
    e0_k_reac_pre = m.varlist_all["e0_k_reac_pre"].casadi_var  # noqa: E501
    e0_r_i4_l1 = m.varlist_all["e0_r_i4_l1"].casadi_var  # noqa: E501
    e0_x_i1_l1 = m.varlist_all["e0_x_i1_l1"].casadi_var  # noqa: E501
    e0_x_i2_l1 = m.varlist_all["e0_x_i2_l1"].casadi_var  # noqa: E501
    e0_x_i3_l1 = m.varlist_all["e0_x_i3_l1"].casadi_var  # noqa: E501
    e0_x_i4_l1 = m.varlist_all["e0_x_i4_l1"].casadi_var  # noqa: E501
    e0_r_i1_l1 = m.varlist_all["e0_r_i1_l1"].casadi_var  # noqa: E501
    e0_r_i2_l1 = m.varlist_all["e0_r_i2_l1"].casadi_var  # noqa: E501
    e0_r_i3_l1 = m.varlist_all["e0_r_i3_l1"].casadi_var  # noqa: E501

    if mode == "gc":
        e0_x_i1_gc = m.varlist_all["e0_x_i1_gc"].casadi_var  # noqa: E501
        e0_x_i2_gc = m.varlist_all["e0_x_i2_gc"].casadi_var  # noqa: E501
        e0_x_i3_gc = m.varlist_all["e0_x_i3_gc"].casadi_var  # noqa: E501
        e0_x_i4_gc = m.varlist_all["e0_x_i4_gc"].casadi_var  # noqa: E501

    EQ_alg1 = (e0_r_i1_l1-((e0_a_Hplus_l1*(e0_greek_nu_i1*(e0_k_for*(((e0_x_i1_l1*(e0_greek_gamma_i1_l1*(e0_x_i2_l1*e0_greek_gamma_i2_l1)))-((e0_x_i3_l1*(e0_greek_gamma_i3_l1*(e0_x_i4_l1*e0_greek_gamma_i4_l1)))/e0_K_reac))))))))  # noqa: E501,E226
    EQ_alg2 = (e0_r_i2_l1-((e0_a_Hplus_l1*(e0_greek_nu_i2*(e0_k_for*(((e0_x_i1_l1*(e0_greek_gamma_i1_l1*(e0_x_i2_l1*e0_greek_gamma_i2_l1)))-((e0_x_i3_l1*(e0_greek_gamma_i3_l1*(e0_x_i4_l1*e0_greek_gamma_i4_l1)))/e0_K_reac))))))))  # noqa: E501,E226
    EQ_alg3 = (e0_r_i3_l1-((e0_a_Hplus_l1*(e0_greek_nu_i3*(e0_k_for*(((e0_x_i1_l1*(e0_greek_gamma_i1_l1*(e0_x_i2_l1*e0_greek_gamma_i2_l1)))-((e0_x_i3_l1*(e0_greek_gamma_i3_l1*(e0_x_i4_l1*e0_greek_gamma_i4_l1)))/e0_K_reac))))))))  # noqa: E501,E226
    EQ_alg4 = (e0_r_i4_l1-((e0_a_Hplus_l1*(e0_greek_nu_i4*(e0_k_for*(((e0_x_i1_l1*(e0_greek_gamma_i1_l1*(e0_x_i2_l1*e0_greek_gamma_i2_l1)))-((e0_x_i3_l1*(e0_greek_gamma_i3_l1*(e0_x_i4_l1*e0_greek_gamma_i4_l1)))/e0_K_reac))))))))  # noqa: E501,E226
    EQ_alg5 = (0.0-((((e0_F_F*e0_x_F_i1)-(e0_F_I*e0_x_i1_l1))+e0_r_i1_l1)))  # noqa: E501,E226
    EQ_alg6 = (0.0-((((e0_F_F*e0_x_F_i2)-(e0_F_I*e0_x_i2_l1))+e0_r_i2_l1)))  # noqa: E501,E226
    EQ_alg7 = (0.0-((((e0_F_F*e0_x_F_i3)-(e0_F_I*e0_x_i3_l1))+e0_r_i3_l1)))  # noqa: E501,E226
    EQ_alg8 = (0.0-((((e0_F_F*e0_x_F_i4)-(e0_F_I*e0_x_i4_l1))+e0_r_i4_l1)))  # noqa: E501,E226
    EQ_alg00 = e0_x_F_i2 + e0_x_F_i1 - 1

    if mode == "gc":
        # sum = e0_x_i1_l1 + e0_x_i2_l1 + e0_x_i3_l1
        sum = e0_x_i1_l1 + e0_x_i2_l1 + e0_x_i3_l1
        EQ_alg01 = e0_x_i1_gc - (e0_x_i1_l1 / sum)
        EQ_alg02 = e0_x_i2_gc - (e0_x_i2_l1 / sum)
        EQ_alg03 = e0_x_i3_gc - (e0_x_i3_l1 / sum)
        EQ_alg04 = e0_x_i4_gc - (e0_x_i4_l1 / sum)
        # EQ_alg01 = e0_x_i1_gc - e0_x_i1_l1
        # EQ_alg02 = e0_x_i2_gc - e0_x_i2_l1
        # EQ_alg03 = e0_x_i3_gc - e0_x_i3_l1
        # EQ_alg04 = e0_x_i4_gc - e0_x_i4_l1


    list_algebraic_equations = []
    # list_algebraic_equations.extend([EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6, EQ_alg7, EQ_alg8, ])  # noqa: E501
    if mode == "gc":
        # list_algebraic_equations.extend([EQ_alg01, EQ_alg02, EQ_alg03])
        list_algebraic_equations.extend([EQ_alg01, EQ_alg02, EQ_alg03, EQ_alg04])
    list_algebraic_equations.extend([EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6, EQ_alg7, EQ_alg8, ])  # noqa: E501
    list_algebraic_equations.append(EQ_alg00)

    # list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6, EQ_alg7, EQ_alg8, ]  # noqa: E501
    try:
        Eq_fun_e0_K_reac = m.varlist_all["e0_K_reac"].casadi_var - fun_175223__arrhenius(e0_T,e0_E_reac_A,e0_R,e0_k_reac_pre)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_K_reac)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_a_Hplus_l1 = m.varlist_all["e0_a_Hplus_l1"].casadi_var - fun_175242__activity_proton(e0_T,e0_x_i4_l1,e0_x_i2_l1,e0_x_cat,e0_Param_Hplus_A,e0_Param_Hplus_B,e0_Param_Hplus_C)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_a_Hplus_l1)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_for = m.varlist_all["e0_k_for"].casadi_var - fun_175223__arrhenius(e0_T,e0_E_for_A,e0_R,e0_k_for_pre)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_for)  # noqa: E501
    except KeyError:
        pass

    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


if __name__ == "__main__":
    mode = "old"
    # mode = "gc"
    variable_list, m = initialize_problem(mode)

    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    sim_fixed = mopeds.SimulatorNLE(m, var_list_fixed)
    res_simple = sim_fixed.simulate_fast()
    res = sim_fixed.simulate()[2]

    unfix_names = [
        "e0_E_for_A",
        # "e0_E_reac_A",
        # "e0_Param_Hplus_A",
        # "e0_Param_Hplus_B",
        # "e0_Param_Hplus_C",
        "e0_k_for_pre",
        # "e0_k_reac_pre",
    ]

    if mode == "old":
        meas_names = [
            "e0_x_i1_l1",
            # "e0_x_i2_l1",
            # "e0_x_i3_l1",
            # "e0_x_i4_l1",
        ]
    elif mode == "gc":
        meas_names = [
            "e0_x_i1_gc",
            "e0_x_i2_gc",
            "e0_x_i3_gc",
        ]

    # meas_names = m.varlist_algebraic.keys()

    controls = {"e0_x_F_i1": (0.2, 0.8, 5), "e0_T": (350, 370, 5)}
    variable_list.set_variable_list_unfixed(unfix_names)
    variable_list.set_bounds(0.001)
    variable_list["e0_k_for_pre"].lower_bound = 1e7
    variable_list["e0_k_for_pre"].upper_bound = 5e7
    variable_list["e0_k_for_pre"].guess = 1e7

    for var_name in meas_names:
        variable_list[var_name].variance = 0.0005**2

    grid = mopeds.tools.create_grid(list(controls.values()))

    rng = np.random.default_rng(0)
    list_varlist, true_par = mopeds.tools.generate_varlist_with_data_NLE(
        m,
        variable_list,
        controls,
        perturbate=True,
        measurement_names=meas_names,
        rng=rng,
    )

    pe = mopeds.ParameterEstimationNLE(m, list_varlist)

    if False:
        residuals = pe.calculate_objective_and_residual(true_par)["residuals"]
        data = np.column_stack((np.array(grid), residuals[:, 0]))
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.scatter(*np.hsplit(data, 3))
        plt.show()
        breakpoint()

    res = pe.optimize(objective_function="wls")
    x_dict = res["x_dict"]

    analysis = pe.parameter_analysis(x_dict)
    for var_name in unfix_names:
        new = x_dict[var_name]
        true = true_par[var_name]
        print(var_name, new, true, f"{(new - true) * 100 / true}%")

    breakpoint()
