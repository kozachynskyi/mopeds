import copy
import casadi as ca
import numpy as np

import mopeds


def initialize_problem():  # noqa: C901

    variable_list = mopeds.VariableList()
    # fmt:off
    def fun_205666__aux_enthalpy_component_vapor(std_T,std_greek_Deltah_V,std_P_V_i1,std_P_V_i2,std_P_V_i3,std_P_V_i4,std_P_V_i5,std_P_V_i6):  # noqa: E501,E231,E306
        std_h = (((((((std_P_V_i1*((std_T/1000.0)))+((std_P_V_i2/2.0)*(((std_T/1000.0)))**(1.0*2.0)))+((std_P_V_i3/3.0)*(((std_T/1000.0)))**(1.0*3.0)))+((std_P_V_i4/4.0)*(((std_T/1000.0)))**(1.0*4.0)))-(std_P_V_i5*(((std_T/1000.0)))**(1.0*(-1.0))))+std_P_V_i6)+std_greek_Deltah_V)  # noqa: E501,E226
        return std_h
    def fun_205665__aux_enthalpy_component_liquid(std_T,std_greek_Deltah_L,std_P_L_i1,std_P_L_i2,std_P_L_i3,std_P_L_i4,std_P_L_i5,std_P_L_i6):  # noqa: E501,E231,E306
        std_h = (((((((std_P_L_i1*((std_T/1000.0)))+((std_P_L_i2/2.0)*(((std_T/1000.0)))**(1.0*2.0)))+((std_P_L_i3/3.0)*(((std_T/1000.0)))**(1.0*3.0)))+((std_P_L_i4/4.0)*(((std_T/1000.0)))**(1.0*4.0)))-(std_P_L_i5*(((std_T/1000.0)))**(1.0*(-1.0))))+std_P_L_i6)+std_greek_Deltah_L)  # noqa: E501,E226
        return std_h

    variable_list.add_variable(mopeds.VariableConstant("e0_greek_Deltah_L_c1", -276000.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_Deltah_V_c1", -234000.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_Deltah_L_c2", -285830.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_Deltah_V_c2", -241830.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_lambda_i1", 1411.64))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_lambda_i2", 4078.0508))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_pi", 3.14))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_tau_LC", 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_tau_PC", 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_tau_TC", 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_c1", 5.24677))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_c2", 5.0768))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_c1", 1598.673))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_c2", 1659.793))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_c1", -46.424))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_c2", -45.854))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_K_LC", 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_K_PC", 10.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_K_TC", -10000.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i1_c1", 102538.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i1_c1", 5385.58))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i1_c2", -203.606))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i1_c2", 30.092))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i2_c1", -138.44))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i2_c1", 236.1088))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i2_c2", 1523.29))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i2_c2", 6.832514))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i3_c1", -0.03469))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i3_c1", 0.1237))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i3_c2", -3196.413))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i3_c2", 6.793435))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i4_c1", 20.4367))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i4_c1", 2.3E-5))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i4_c2", 2474.455))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i4_c2", -2.53448))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i5_c1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i5_c1", 3.7E-5))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i5_c2", 3.855326))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i5_c2", 0.082139))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i6_c1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i6_c1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_L_i6_c2", -256.5478))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_V_i6_c2", -250.881))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_R", 0.0083145))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_v_L_c1", 1.66E-4))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_v_L_c2", 5.58E-5))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_F_c1", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_F_c2", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_L_c1", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_L_c2", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_V_c1", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_V_c2", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_H_L", 2.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_Q", 395124.7, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p", 0.8, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_V", 11.8988, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_Lambda_i1", 0.2078426, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_Lambda_i2", 0.7418247, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_gamma_c1", 3.407189, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_gamma_c2", 1.029053, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_L", 23590.7, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_V", 2.138771E-4, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_V_L", 1.57, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_V_Tank", 2.355, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_V_V", 0.785, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_F", -279985.1, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_L", 68.1012, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_L", -281593.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_V", -237574.9, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_LV_o_c1", 1.083214, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_LV_o_c2", 0.4737136, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_c1", 0.09756492, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_c2", 0.9024351, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_F_c2", 0.85, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c1", 0.4501054, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c2", 0.5498946, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_T", 353.15, -1.0E9, 1.0E9))  # noqa: E501

    variable_list.add_variable(mopeds.VariableState("e0_I_PC", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableState("e0_I_LC", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableState("e0_HU_c1", 2301.625))  # noqa: E501
    variable_list.add_variable(mopeds.VariableState("e0_HU_c2", 21289.07))  # noqa: E501
    variable_list.add_variable(mopeds.VariableState("e0_h", -6.642976E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableState("e0_I_TC", 0.0))  # noqa: E501

    variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_H_L_set", 2.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_T_set", 353.15, -9.9999972685E8, 1.00000027315E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_Q_set", 400000.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_p_set", 0.8, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_V_set", 11.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_D", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_F", 80.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_H", 3.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_T_F", 300.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_L_set", 68.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_x_F_c1", 0.15, -1.0E9, 1.0E9))  # noqa: E501


    m = mopeds.Model(variable_list)

    e0_H_L_set = m.varlist_all["e0_H_L_set"].casadi_var  # noqa: E501
    e0_T_set = m.varlist_all["e0_T_set"].casadi_var  # noqa: E501
    e0_Q_set = m.varlist_all["e0_Q_set"].casadi_var  # noqa: E501
    e0_p_set = m.varlist_all["e0_p_set"].casadi_var  # noqa: E501
    e0_V_set = m.varlist_all["e0_V_set"].casadi_var  # noqa: E501
    e0_D = m.varlist_all["e0_D"].casadi_var  # noqa: E501
    e0_F = m.varlist_all["e0_F"].casadi_var  # noqa: E501
    e0_H = m.varlist_all["e0_H"].casadi_var  # noqa: E501
    e0_T_F = m.varlist_all["e0_T_F"].casadi_var  # noqa: E501
    e0_L_set = m.varlist_all["e0_L_set"].casadi_var  # noqa: E501
    e0_x_F_c1 = m.varlist_all["e0_x_F_c1"].casadi_var  # noqa: E501
    e0_h_F_c1 = m.varlist_all["e0_h_F_c1"].casadi_var  # noqa: E501
    e0_h_F_c2 = m.varlist_all["e0_h_F_c2"].casadi_var  # noqa: E501
    e0_h_L_c1 = m.varlist_all["e0_h_L_c1"].casadi_var  # noqa: E501
    e0_h_L_c2 = m.varlist_all["e0_h_L_c2"].casadi_var  # noqa: E501
    e0_h_V_c1 = m.varlist_all["e0_h_V_c1"].casadi_var  # noqa: E501
    e0_h_V_c2 = m.varlist_all["e0_h_V_c2"].casadi_var  # noqa: E501
    e0_greek_Deltah_L_c1 = m.varlist_all["e0_greek_Deltah_L_c1"].casadi_var  # noqa: E501
    e0_greek_Deltah_V_c1 = m.varlist_all["e0_greek_Deltah_V_c1"].casadi_var  # noqa: E501
    e0_greek_Deltah_L_c2 = m.varlist_all["e0_greek_Deltah_L_c2"].casadi_var  # noqa: E501
    e0_greek_Deltah_V_c2 = m.varlist_all["e0_greek_Deltah_V_c2"].casadi_var  # noqa: E501
    e0_greek_lambda_i1 = m.varlist_all["e0_greek_lambda_i1"].casadi_var  # noqa: E501
    e0_greek_lambda_i2 = m.varlist_all["e0_greek_lambda_i2"].casadi_var  # noqa: E501
    e0_greek_pi = m.varlist_all["e0_greek_pi"].casadi_var  # noqa: E501
    e0_greek_tau_LC = m.varlist_all["e0_greek_tau_LC"].casadi_var  # noqa: E501
    e0_greek_tau_PC = m.varlist_all["e0_greek_tau_PC"].casadi_var  # noqa: E501
    e0_greek_tau_TC = m.varlist_all["e0_greek_tau_TC"].casadi_var  # noqa: E501
    e0_A_c1 = m.varlist_all["e0_A_c1"].casadi_var  # noqa: E501
    e0_A_c2 = m.varlist_all["e0_A_c2"].casadi_var  # noqa: E501
    e0_B_c1 = m.varlist_all["e0_B_c1"].casadi_var  # noqa: E501
    e0_B_c2 = m.varlist_all["e0_B_c2"].casadi_var  # noqa: E501
    e0_C_c1 = m.varlist_all["e0_C_c1"].casadi_var  # noqa: E501
    e0_C_c2 = m.varlist_all["e0_C_c2"].casadi_var  # noqa: E501
    e0_K_LC = m.varlist_all["e0_K_LC"].casadi_var  # noqa: E501
    e0_K_PC = m.varlist_all["e0_K_PC"].casadi_var  # noqa: E501
    e0_K_TC = m.varlist_all["e0_K_TC"].casadi_var  # noqa: E501
    e0_P_L_i1_c1 = m.varlist_all["e0_P_L_i1_c1"].casadi_var  # noqa: E501
    e0_P_V_i1_c1 = m.varlist_all["e0_P_V_i1_c1"].casadi_var  # noqa: E501
    e0_P_L_i1_c2 = m.varlist_all["e0_P_L_i1_c2"].casadi_var  # noqa: E501
    e0_P_V_i1_c2 = m.varlist_all["e0_P_V_i1_c2"].casadi_var  # noqa: E501
    e0_P_L_i2_c1 = m.varlist_all["e0_P_L_i2_c1"].casadi_var  # noqa: E501
    e0_P_V_i2_c1 = m.varlist_all["e0_P_V_i2_c1"].casadi_var  # noqa: E501
    e0_P_L_i2_c2 = m.varlist_all["e0_P_L_i2_c2"].casadi_var  # noqa: E501
    e0_P_V_i2_c2 = m.varlist_all["e0_P_V_i2_c2"].casadi_var  # noqa: E501
    e0_P_L_i3_c1 = m.varlist_all["e0_P_L_i3_c1"].casadi_var  # noqa: E501
    e0_P_V_i3_c1 = m.varlist_all["e0_P_V_i3_c1"].casadi_var  # noqa: E501
    e0_P_L_i3_c2 = m.varlist_all["e0_P_L_i3_c2"].casadi_var  # noqa: E501
    e0_P_V_i3_c2 = m.varlist_all["e0_P_V_i3_c2"].casadi_var  # noqa: E501
    e0_P_L_i4_c1 = m.varlist_all["e0_P_L_i4_c1"].casadi_var  # noqa: E501
    e0_P_V_i4_c1 = m.varlist_all["e0_P_V_i4_c1"].casadi_var  # noqa: E501
    e0_P_L_i4_c2 = m.varlist_all["e0_P_L_i4_c2"].casadi_var  # noqa: E501
    e0_P_V_i4_c2 = m.varlist_all["e0_P_V_i4_c2"].casadi_var  # noqa: E501
    e0_P_L_i5_c1 = m.varlist_all["e0_P_L_i5_c1"].casadi_var  # noqa: E501
    e0_P_V_i5_c1 = m.varlist_all["e0_P_V_i5_c1"].casadi_var  # noqa: E501
    e0_P_L_i5_c2 = m.varlist_all["e0_P_L_i5_c2"].casadi_var  # noqa: E501
    e0_P_V_i5_c2 = m.varlist_all["e0_P_V_i5_c2"].casadi_var  # noqa: E501
    e0_P_L_i6_c1 = m.varlist_all["e0_P_L_i6_c1"].casadi_var  # noqa: E501
    e0_P_V_i6_c1 = m.varlist_all["e0_P_V_i6_c1"].casadi_var  # noqa: E501
    e0_P_L_i6_c2 = m.varlist_all["e0_P_L_i6_c2"].casadi_var  # noqa: E501
    e0_P_V_i6_c2 = m.varlist_all["e0_P_V_i6_c2"].casadi_var  # noqa: E501
    e0_R = m.varlist_all["e0_R"].casadi_var  # noqa: E501
    e0_v_L_c1 = m.varlist_all["e0_v_L_c1"].casadi_var  # noqa: E501
    e0_v_L_c2 = m.varlist_all["e0_v_L_c2"].casadi_var  # noqa: E501
    e0_H_L = m.varlist_all["e0_H_L"].casadi_var  # noqa: E501
    e0_Q = m.varlist_all["e0_Q"].casadi_var  # noqa: E501
    e0_I_PC = m.varlist_all["e0_I_PC"].casadi_var  # noqa: E501
    e0_p = m.varlist_all["e0_p"].casadi_var  # noqa: E501
    e0_V = m.varlist_all["e0_V"].casadi_var  # noqa: E501
    e0_I_LC = m.varlist_all["e0_I_LC"].casadi_var  # noqa: E501
    e0_greek_Lambda_i1 = m.varlist_all["e0_greek_Lambda_i1"].casadi_var  # noqa: E501
    e0_greek_Lambda_i2 = m.varlist_all["e0_greek_Lambda_i2"].casadi_var  # noqa: E501
    e0_greek_gamma_c1 = m.varlist_all["e0_greek_gamma_c1"].casadi_var  # noqa: E501
    e0_greek_gamma_c2 = m.varlist_all["e0_greek_gamma_c2"].casadi_var  # noqa: E501
    e0_HU_c1 = m.varlist_all["e0_HU_c1"].casadi_var  # noqa: E501
    e0_HU_c2 = m.varlist_all["e0_HU_c2"].casadi_var  # noqa: E501
    e0_HU_L = m.varlist_all["e0_HU_L"].casadi_var  # noqa: E501
    e0_HU_V = m.varlist_all["e0_HU_V"].casadi_var  # noqa: E501
    e0_V_L = m.varlist_all["e0_V_L"].casadi_var  # noqa: E501
    e0_V_Tank = m.varlist_all["e0_V_Tank"].casadi_var  # noqa: E501
    e0_V_V = m.varlist_all["e0_V_V"].casadi_var  # noqa: E501
    e0_h = m.varlist_all["e0_h"].casadi_var  # noqa: E501
    e0_h_F = m.varlist_all["e0_h_F"].casadi_var  # noqa: E501
    e0_L = m.varlist_all["e0_L"].casadi_var  # noqa: E501
    e0_h_L = m.varlist_all["e0_h_L"].casadi_var  # noqa: E501
    e0_h_V = m.varlist_all["e0_h_V"].casadi_var  # noqa: E501
    e0_p_LV_o_c1 = m.varlist_all["e0_p_LV_o_c1"].casadi_var  # noqa: E501
    e0_p_LV_o_c2 = m.varlist_all["e0_p_LV_o_c2"].casadi_var  # noqa: E501
    e0_x_c1 = m.varlist_all["e0_x_c1"].casadi_var  # noqa: E501
    e0_x_c2 = m.varlist_all["e0_x_c2"].casadi_var  # noqa: E501
    e0_x_F_c2 = m.varlist_all["e0_x_F_c2"].casadi_var  # noqa: E501
    e0_y_c1 = m.varlist_all["e0_y_c1"].casadi_var  # noqa: E501
    e0_y_c2 = m.varlist_all["e0_y_c2"].casadi_var  # noqa: E501
    e0_I_TC = m.varlist_all["e0_I_TC"].casadi_var  # noqa: E501
    e0_T = m.varlist_all["e0_T"].casadi_var  # noqa: E501

    EQ_diff1 = (e0_H_L_set-e0_H_L)  # noqa: E501,E226
    EQ_diff2 = (e0_T_set-e0_T)  # noqa: E501,E226
    EQ_diff3 = (e0_p_set-e0_p)  # noqa: E501,E226
    EQ_diff4 = (((e0_F*e0_x_F_c1)-(e0_V*e0_y_c1))-(e0_L*e0_x_c1))  # noqa: E501,E226
    EQ_diff5 = (((e0_F*e0_x_F_c2)-(e0_V*e0_y_c2))-(e0_L*e0_x_c2))  # noqa: E501,E226
    EQ_diff6 = ((((e0_F*e0_h_F)-(e0_V*e0_h_V))-(e0_L*e0_h_L))+e0_Q)  # noqa: E501,E226

    EQ_alg7 = (e0_L-((e0_L_set-(e0_K_LC*((((e0_H_L_set-e0_H_L))+(e0_I_LC/e0_greek_tau_LC)))))))  # noqa: E501,E226
    EQ_alg8 = (e0_Q-((e0_Q_set-(e0_K_TC*((((e0_T_set-e0_T))+(e0_I_TC/e0_greek_tau_TC)))))))  # noqa: E501,E226
    EQ_alg9 = (e0_V-((e0_V_set-(e0_K_PC*((((e0_p_set-e0_p))+(e0_I_PC/e0_greek_tau_PC)))))))  # noqa: E501,E226
    EQ_alg10 = (e0_HU_c1-(((e0_x_c1*e0_HU_L)+(e0_y_c1*e0_HU_V))))  # noqa: E501,E226
    EQ_alg11 = (e0_HU_c2-(((e0_x_c2*e0_HU_L)+(e0_y_c2*e0_HU_V))))  # noqa: E501,E226
    EQ_alg12 = (e0_h-(((e0_HU_L*e0_h_L)+(e0_HU_V*e0_h_V))))  # noqa: E501,E226
    EQ_alg13 = (e0_H_L-((4.0*(e0_V_L/(e0_greek_pi*((e0_D))**(1.0*2.0))))))  # noqa: E501,E226
    EQ_alg14 = (e0_V_Tank-((e0_V_L+e0_V_V)))  # noqa: E501,E226
    EQ_alg15 = (e0_V_L-((((((e0_v_L_c1*e0_x_c1)+(e0_v_L_c2*e0_x_c2))))*e0_HU_L)))  # noqa: E501,E226
    EQ_alg16 = (e0_V_V-(((e0_HU_V*(e0_R*(1000.0*e0_T)))/e0_p)))  # noqa: E501,E226
    EQ_alg17 = (e0_V_Tank-((e0_greek_pi*((((e0_D))**(1.0*2.0)/4.0)*e0_H))))  # noqa: E501,E226
    EQ_alg18 = ((e0_y_c1*e0_p)-((e0_x_c1*(e0_greek_gamma_c1*e0_p_LV_o_c1))))  # noqa: E501,E226
    EQ_alg19 = ((e0_y_c2*e0_p)-((e0_x_c2*(e0_greek_gamma_c2*e0_p_LV_o_c2))))  # noqa: E501,E226
    EQ_alg20 = (((e0_x_F_c1+e0_x_F_c2))-(1.0))  # noqa: E501,E226
    EQ_alg21 = (((e0_x_c1+e0_x_c2))-(1.0))  # noqa: E501,E226
    EQ_alg22 = (((e0_y_c1+e0_y_c2))-(1.0))  # noqa: E501,E226
    EQ_alg23 = (e0_p_LV_o_c1-(((10.0))**(1.0*(e0_A_c1-(e0_B_c1/(e0_T+e0_C_c1))))))  # noqa: E501,E226
    EQ_alg24 = (e0_p_LV_o_c2-(((10.0))**(1.0*(e0_A_c2-(e0_B_c2/(e0_T+e0_C_c2))))))  # noqa: E501,E226
    EQ_alg25 = (e0_h_F-((((e0_x_F_c1*e0_h_F_c1)+(e0_x_F_c2*e0_h_F_c2)))))  # noqa: E501,E226
    EQ_alg26 = (e0_h_L-((((e0_x_c1*e0_h_L_c1)+(e0_x_c2*e0_h_L_c2)))))  # noqa: E501,E226
    EQ_alg27 = (e0_h_V-((((e0_y_c1*e0_h_V_c1)+(e0_y_c2*e0_h_V_c2)))))  # noqa: E501,E226
    EQ_alg28 = (e0_greek_gamma_c1-(ca.exp(((-ca.log((e0_x_c1+(e0_greek_Lambda_i1*e0_x_c2))))+(e0_x_c2*(((e0_greek_Lambda_i1/(e0_x_c1+(e0_greek_Lambda_i1*e0_x_c2)))-(e0_greek_Lambda_i2/(e0_x_c2+(e0_greek_Lambda_i2*e0_x_c1))))))))))  # noqa: E501,E226
    EQ_alg29 = (e0_greek_gamma_c2-(ca.exp(((-ca.log((e0_x_c2+(e0_greek_Lambda_i2*e0_x_c1))))-(e0_x_c1*(((e0_greek_Lambda_i1/(e0_x_c1+(e0_greek_Lambda_i1*e0_x_c2)))-(e0_greek_Lambda_i2/(e0_x_c2+(e0_greek_Lambda_i2*e0_x_c1))))))))))  # noqa: E501,E226
    EQ_alg30 = (e0_greek_Lambda_i1-(((e0_v_L_c2/e0_v_L_c1)*ca.exp((-(e0_greek_lambda_i1/(e0_R*(e0_T*1000.0))))))))  # noqa: E501,E226
    EQ_alg31 = (e0_greek_Lambda_i2-(((e0_v_L_c1/e0_v_L_c2)*ca.exp((-(e0_greek_lambda_i2/(e0_R*(e0_T*1000.0))))))))  # noqa: E501,E226

    order_state_var = ["e0_I_PC", "e0_I_LC", "e0_HU_c1", "e0_HU_c2", "e0_h", "e0_I_TC", ]  # noqa: E501
    order_eqs_diff = {"e0_I_LC": EQ_diff1, "e0_I_TC": EQ_diff2, "e0_I_PC": EQ_diff3, "e0_HU_c1": EQ_diff4, "e0_HU_c2": EQ_diff5, "e0_h": EQ_diff6, }  # noqa: E501

    order_eqs_corrected = list(order_eqs_diff[i] for i in order_state_var)
    m.add_equations_differential(order_eqs_corrected)
    list_algebraic_equations = [EQ_alg7, EQ_alg8, EQ_alg9, EQ_alg10, EQ_alg11, EQ_alg12, EQ_alg13, EQ_alg14, EQ_alg15, EQ_alg16, EQ_alg17, EQ_alg18, EQ_alg19, EQ_alg20, EQ_alg21, EQ_alg22, EQ_alg23, EQ_alg24, EQ_alg25, EQ_alg26, EQ_alg27, EQ_alg28, EQ_alg29, EQ_alg30, EQ_alg31, ]  # noqa: E501

    try:
        Eq_fun_e0_h_L_c2 = m.varlist_all["e0_h_L_c2"].casadi_var - fun_205665__aux_enthalpy_component_liquid(e0_T,e0_greek_Deltah_L_c2,e0_P_L_i1_c2,e0_P_L_i2_c2,e0_P_L_i3_c2,e0_P_L_i4_c2,e0_P_L_i5_c2,e0_P_L_i6_c2)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_L_c2)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_F_c1 = m.varlist_all["e0_h_F_c1"].casadi_var - fun_205665__aux_enthalpy_component_liquid(e0_T_F,e0_greek_Deltah_L_c1,e0_P_L_i1_c1,e0_P_L_i2_c1,e0_P_L_i3_c1,e0_P_L_i4_c1,e0_P_L_i5_c1,e0_P_L_i6_c1)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_F_c1)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_V_c1 = m.varlist_all["e0_h_V_c1"].casadi_var - fun_205666__aux_enthalpy_component_vapor(e0_T,e0_greek_Deltah_V_c1,e0_P_V_i1_c1,e0_P_V_i2_c1,e0_P_V_i3_c1,e0_P_V_i4_c1,e0_P_V_i5_c1,e0_P_V_i6_c1)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_V_c1)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_V_c2 = m.varlist_all["e0_h_V_c2"].casadi_var - fun_205666__aux_enthalpy_component_vapor(e0_T,e0_greek_Deltah_V_c2,e0_P_V_i1_c2,e0_P_V_i2_c2,e0_P_V_i3_c2,e0_P_V_i4_c2,e0_P_V_i5_c2,e0_P_V_i6_c2)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_V_c2)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_F_c2 = m.varlist_all["e0_h_F_c2"].casadi_var - fun_205665__aux_enthalpy_component_liquid(e0_T_F,e0_greek_Deltah_L_c2,e0_P_L_i1_c2,e0_P_L_i2_c2,e0_P_L_i3_c2,e0_P_L_i4_c2,e0_P_L_i5_c2,e0_P_L_i6_c2)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_F_c2)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_L_c1 = m.varlist_all["e0_h_L_c1"].casadi_var - fun_205665__aux_enthalpy_component_liquid(e0_T,e0_greek_Deltah_L_c1,e0_P_L_i1_c1,e0_P_L_i2_c1,e0_P_L_i3_c1,e0_P_L_i4_c1,e0_P_L_i5_c1,e0_P_L_i6_c1)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_L_c1)  # noqa: E501
    except KeyError:
        pass

    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


if __name__ == "__main__":
    var_list, m = initialize_problem()

    # Create time-grid. Zero should be first
    for var in var_list.values():
        if isinstance(var, mopeds.VariableState):
            var.ignore_plotting = True
        elif isinstance(var, mopeds.VariableControlPiecewiseConstant):
            var.ignore_plotting = False

    var_list["e0_H_L"].ignore_plotting = False
    var_list["e0_L"].ignore_plotting = False
    var_list["e0_T"].ignore_plotting = False
    var_list["e0_Q"].ignore_plotting = False
    var_list["e0_p"].ignore_plotting = False
    var_list["e0_V"].ignore_plotting = False
    var_list["e0_H_L_set"].expand_horizon([5000], [1.5])
    # var_list["e0_H_L_set"].expand_horizon([5000], [1.5])
    var_list["e0_p_set"].expand_horizon([10000], [0.82])
    var_list["e0_T_set"].expand_horizon([15000], [85+273.15])
    var_list["e0_F"].expand_horizon([20000], [85])
    var_list["e0_x_F_c1"].expand_horizon([25000], [0.16])

    time_grid = np.linspace(0, 30000, 1000)

    # Create simulation Object
    sim = mopeds.Simulator(m, time_grid, var_list)
    res = sim.generate_exp_data(algebraic=True)
    print(res.dataframe)
    res.plot(algebraic=True)
