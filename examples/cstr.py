# import copy
# from collections import OrderedDict
# from datetime import datetime, timedelta

# import casadi as ca
# import matplotlib.cm as cm
# import numpy as np
# from matplotlib import pyplot as plt
# from opcua import ua
# from opcua.ua import NumericNodeId
# from optipal.client import OptiPALClient

def main():
    e0_greek_nu_i1_r1 = -1.0
    e0_greek_nu_i1_r2 = 1.0
    e0_greek_nu_i2_r2 = -1.0
    e0_greek_nu_i3_r1 = 1.0
    e0_greek_nu_i1_r3 = -1.0
    e0_greek_nu_i4_r3 = 1.0
    e0_greek_Deltah_r1 = 0.0045
    e0_greek_Deltah_r2 = -0.0055
    e0_greek_Deltah_r3 = 0.0045
    e0_greek_rho = 800.0
    e0_A = 1.0
    e0_E_r1 = 96000.0
    e0_c_p = 3.5
    e0_E_r2 = 72000.0
    e0_E_r3 = 69000.0
    e0_F = 6.5e-4
    e0_R = 8.314
    e0_V = 1.0

    variable_list = VariableList()

    variable_list.add_variable(State_variable("e0_T", 273.0, 10))
    variable_list.add_variable(State_variable("e0_c_i1", 3.0, 20))
    variable_list.add_variable(State_variable("e0_c_i2", 10.0, 30))
    variable_list.add_variable(State_variable("e0_c_i3", 0.0, 40))
    variable_list.add_variable(State_variable("e0_c_i4", 0.0, 50))

    variable_list.add_variable(Parameter_variable("e0_k_pre_r1", 5000000.0))
    variable_list.add_variable(Parameter_variable("e0_k_pre_r2", 1.0e7))
    variable_list.add_variable(Parameter_variable("e0_k_pre_r3", 500000.0))
    variable_list.add_variable(Parameter_variable("e0_U", 1.4))

    variable_list.add_variable(Control_variable("e0_c_in_i1", 5.0))
    variable_list.add_variable(Control_variable("e0_c_in_i2", 10.0))
    variable_list.add_variable(Control_variable("e0_c_in_i3", 0.0))
    variable_list.add_variable(Control_variable("e0_c_in_i4", 0.0))
    variable_list.add_variable(Control_variable("e0_T_in", 373.0))
    variable_list.add_variable(Control_variable("e0_T_j", 373.0))

    m = Model(variable_list)

    # fmt: off
    tdot = (((((e0_F / e0_V) * ((m._all_variables["e0_T_in"].casadi_var - m._all_variables["e0_T"].casadi_var))) + (((m._all_variables["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (e0_c_p * e0_V))) * ((m._all_variables["e0_T_j"].casadi_var - m._all_variables["e0_T"].casadi_var)))) + (((-e0_greek_Deltah_r1) / (e0_greek_rho * e0_c_p)) * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-e0_greek_Deltah_r2) / (e0_greek_rho * e0_c_p)) * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-e0_greek_Deltah_r3) / (e0_greek_rho * e0_c_p)) * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c1dot = ((((e0_F / e0_V) * ((m._all_variables["e0_c_in_i1"].casadi_var - m._all_variables["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c2dot = ((e0_F / e0_V) * ((m._all_variables["e0_c_in_i2"].casadi_var - m._all_variables["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c3dot = ((e0_F / e0_V) * ((m._all_variables["e0_c_in_i3"].casadi_var - m._all_variables["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c4dot = ((e0_F / e0_V) * ((m._all_variables["e0_c_in_i4"].casadi_var - m._all_variables["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    # fmt: on

    m.add_equations([tdot, c1dot, c2dot, c3dot, c4dot])

    time_grid = np.linspace(10, 1000, 20)
    time_grid = np.insert(time_grid, 0, 0)

    var_list1 = copy.deepcopy(variable_list)
    for var in var_list1.values():
        var.fixed = True
    s = Simulator(m, time_grid, var_list1)
    res = s.simulate()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="\t")

    var_list_exp = s.generate_exp_data()

    var_list2 = copy.deepcopy(var_list1)
    for key, var in var_list_exp.items():
        var_list2[key] = var

    # var_list2["e0_T"].value = Experimental_Data()
    # var_list2["e0_c_i4"].value = Experimental_Data()

    var_list2["e0_k_pre_r1"].fixed = True
    var_list2["e0_k_pre_r1"].guess = 4000000.0
    var_list2["e0_k_pre_r1"].lower_bound = 4000000.0
    var_list2["e0_k_pre_r1"].upper_bound = 6000000.0

    var_list2["e0_k_pre_r2"].fixed = False
    var_list2["e0_k_pre_r2"].guess = 1.0e6
    var_list2["e0_k_pre_r2"].lower_bound = 1.0e6
    var_list2["e0_k_pre_r2"].upper_bound = 1.0e8

    var_list2["e0_k_pre_r3"].fixed = True
    var_list2["e0_k_pre_r3"].guess = 400000.0
    var_list2["e0_k_pre_r3"].lower_bound = 400000.0
    var_list2["e0_k_pre_r3"].upper_bound = 600000.0

    var_list2["e0_U"].fixed = False
    var_list2["e0_U"].guess = 1.1
    var_list2["e0_U"].lower_bound = 1.0
    var_list2["e0_U"].upper_bound = 3.0

    var_list2["e0_c_in_i1"].fixed = False
    var_list2["e0_c_in_i1"].guess = 5.0
    var_list2["e0_c_in_i1"].lower_bound = 4.0
    var_list2["e0_c_in_i1"].upper_bound = 6.0

    var_list2["e0_c_in_i2"].fixed = False
    var_list2["e0_c_in_i2"].guess = 10.0
    var_list2["e0_c_in_i2"].lower_bound = 9.0
    var_list2["e0_c_in_i2"].upper_bound = 11.0

    var_list2["e0_c_in_i3"].fixed = False
    var_list2["e0_c_in_i3"].guess = 0.0
    var_list2["e0_c_in_i3"].lower_bound = 0.0
    var_list2["e0_c_in_i3"].upper_bound = 0.0

    var_list2["e0_c_in_i4"].fixed = False
    var_list2["e0_c_in_i4"].guess = 0.0
    var_list2["e0_c_in_i4"].lower_bound = 0.0
    var_list2["e0_c_in_i4"].upper_bound = 0.0

    var_list2["e0_T_in"].fixed = False
    var_list2["e0_T_in"].guess = 373.0
    var_list2["e0_T_in"].lower_bound = 353.0
    var_list2["e0_T_in"].upper_bound = 393.0

    var_list2["e0_T_j"].fixed = False
    var_list2["e0_T_j"].guess = 373.0
    var_list2["e0_T_j"].lower_bound = 353.0
    var_list2["e0_T_j"].upper_bound = 393.0

    start_time = datetime(2018, 1, 1, 1, 0, 0, 0) + timedelta(days=1)
    end_time = start_time + timedelta(seconds=1000)
    var_list_oed = copy.deepcopy(var_list2)
    # var_list3 = copy.deepcopy(var_list2)
    # var_list_exp.write_data_opcua(start_time)
    # var_list3.get_data_opcua(start_time, end_time)
    pe = ParameterEstimation(m, var_list2)
    # pe.optimize()

    oed = OptimalExperimentalDesign(m, var_list_oed, time_grid)
    a = oed.get_fim_matrix()
    b = a[0].toarray()
    # b = ca.fabs(b)
    c = a[1].toarray()
    # turn_off_states = np.array([1, 1, 1, 1, 1])
    # sc_states = [1, 0.01, 0.01, 0.01, 0.01]
    # sc = np.diagflat(np.tile(turn_off_states, len(time_grid) - 1))
    # sc_full = np.diagflat(np.tile(turn_off_states / sc_states, len(time_grid) - 1))
    # # num_states = 5
    # # sc_full_params = np.tile(sc_params, ((len(time_grid) - 1) * num_states, 1)).T
    # sc_params = [5000000.0, 10000000.0, 500000.0, 1.4]
    # sc_full_params = np.diagflat(sc_params)
    # b_scaled = sc @ b
    # b_scaled_full = sc_full @ (b @ sc_full_params)
    # # sc = np.tile(sc, 2)
    fig = plt.figure()
    fig.add_subplot(151).imshow(b, cmap=cm.Greens_r)
    fig.add_subplot(152).imshow(ca.inv(c), cmap=cm.Greens_r)
    plt.show()
    oed.optimize()
    # pe.optimize(False)

if __name__ == "__main__":
    main()

