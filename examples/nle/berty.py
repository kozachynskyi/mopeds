import copy

import casadi as ca

import mopeds


def initialize_problem():  # noqa: C901

    variable_list = mopeds.VariableList()
    # fmt:off
    def fun_113237__enthalpyFNC(std_T,std_A_par1,std_A_par2,std_A_par3,std_A_par4,std_A_par5,std_T_f,std_h_f):  # noqa: E501,E231,E306
        std_h = (((((std_h_f+(std_A_par1*((((std_T))**(1.0*5.0)-((std_T_f))**(1.0*5.0))/5.0)))+(std_A_par2*((((std_T))**(1.0*4.0)-((std_T_f))**(1.0*4.0))/4.0)))+(std_A_par3*((((std_T))**(1.0*3.0)-((std_T_f))**(1.0*3.0))/3.0)))+(std_A_par4*((((std_T))**(1.0*2.0)-((std_T_f))**(1.0*2.0))/2.0)))+(std_A_par5*((((std_T))**(1.0*1.0)-((std_T_f))**(1.0*1.0))/1.0)))  # noqa: E501,E226
        return std_h
    def fun_113291__ArrheniusFNC(std_T,std_A,std_E,std_R):  # noqa: E501,E231,E306
        std_k = (std_A*ca.exp((std_E/(std_R*std_T))))  # noqa: E501,E226
        return std_k

    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i2", -1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i3", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i4", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i5", -1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i6", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r1_i7", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i2", -3.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i3", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i4", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i5", -1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i6", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r2_i7", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i3", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i4", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i5", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i6", -2.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_greek_nu_r3_i7", 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i1", 5.45E-12))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i2", -7.52E-12))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i3", 1.07368E-10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i4", -1.54836E-11))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i5", 1.93148E-11))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i6", 2.65055E-10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par1_i7", 3.87561E-10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i1", -2.44E-8))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i2", 2.7E-8))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i3", -2.50993E-7))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i4", 2.10847E-8))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i5", -1.32293E-8))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i6", -6.5285E-7))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par2_i7", -8.56412E-7))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i1", 3.7E-5))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i2", -3.17E-5))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i3", 2.2E-4))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i4", 6.07E-7))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i5", -4.18E-5))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i6", 5.57E-4))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par3_i7", 6.29E-4))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i1", -0.0169))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i2", 0.0162))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i3", -0.074))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i4", -0.00322))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i5", 0.0718))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i6", -0.121))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par4_i7", -0.0526))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i1", 31.5))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i2", 26.2))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i3", 42.04061276))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i4", 29.59614774))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i5", 19.6353642))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i6", 45.78156724))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_par5_i7", 45.14491836))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_eq_r1_par1", -2073.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_eq_r1_par2", 2.029))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_eq_r2_par1", 3066.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_eq_r2_par2", -10.592))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_R", 8.314))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_T_f", 298.15))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i3", -241830.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i4", -110530.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i5", -393520.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i6", -200940.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_h_f_i7", -184100.0))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_K_i3", 6.486448694167284E9, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_K_i6", 7.155393165647614E13, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i1", -145.74383548430316, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i2", -144.29466839718256, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i3", -241998.63319940824, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i4", -110675.6190793542, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i5", -393705.0747228285, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i6", -201159.24396779988, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j1_i7", -184427.01699715087, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i1", 5235.161082479042, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i2", 5208.955153408527, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i3", -235692.84107834837, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i4", -105276.82434759167, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i5", -386226.31987338874, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i6", -191876.49529865966, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_o_j2_i7", -170120.58877928287, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_par1", 25611.10018365477, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_par2", 3453.38, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_par3", 38.085637891091515, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_par4", 2564.03213761555, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_par5", 0.2231027466231292, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_k_r3", 0.07991603396072855, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_K_r2", 6.8179E-5, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_X_i5", 0.050537, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_X_i2", 0.042196, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_Y_i7", 1.0609E-7, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_K_r3", 30.1167, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i1_j2", 0.20306, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i2_j2", 0.5819, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i3_j2", 0.010345, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i4_j2", 0.0026993, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i5_j1", 0.20161, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i5_j2", 0.19435, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i6_j2", 0.0076454, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_i7_j2", 2.1718E-8, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_j1", -0.079491, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_h_j2", -0.075159, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_K_r1", 0.0048356, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i1_j2", 2.1321, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i2_j2", 6.11, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i3_j2", 0.10862, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i4_j2", 0.028343, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i5_j2", 2.0407, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i6_j2", 0.080277, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_i7_j2", 2.2804E-7, 0.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU", 1.3234E-5, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_U", -0.0010472, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_u_j2", -0.079126, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i1", 2.6873E-6, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i2", 7.7009E-6, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i3", 1.369E-7, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i4", 3.5723E-8, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i5", 2.5721E-6, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i6", 1.0118E-7, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_HU_i7", 2.8741E-13, 0.0, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_v_j2", 0.0037781, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_r3", 3.9948E-12, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_r2", 1.4064E-6, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_r_r1", 4.9655E-7, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_F_j2", -1.8395E-4, -1.0E9, 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_Q", 0.0010206, -1.0E9, 1.0E9))  # noqa: E501

    variable_list.add_variable(mopeds.VariableControl("e0_T_j2", 477.15, 300.0, 800.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_A_r3_par4", 3.8E-7, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_T_j1", 293.15, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_A_r3_par5", -65610.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_A_r3_par6", -26.64, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_i1_j1", 0.2, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_i2_j1", 0.598387352, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_i3_j1", 1.0E-12, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_i4_j1", 1.0E-12, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_i6_j1", 1.0E-12, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_i7_j1", 1.0E-12, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_p_j2", 10.5, 1.0, 100.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_V", 5.0E-5, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_rho_r3", 1.0E-6, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_A_r3_par1", 4019.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_rho_r1", 10.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_F_j1", 1.86765E-4, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_A_r3_par2", 3.707, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_A_r3_par3", -0.002783, -1.0E9, 1.0E9))  # noqa: E501

    variable_list.add_variable(mopeds.VariableParameter("e0_A_par1_r1", 1.07, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_par2_r1", 3453.38, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_par3_r1", 0.499, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_par4_r1", 6.62E-11, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_par5_r1", 1.22E10, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_r3_i3", 0.5498, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_r3_i6", 223.2, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_r3", 85190.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_par1_r1", 40000.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_par2_r1", 0.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_par3_r1", 17197.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_par4_r1", 124119.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_par5_r1", -98084.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_r3_i3", 92000.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_r3_i6", 105100.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_E_r3", -55060.0, -1.0E9, 1.0E9))  # noqa: E501

    m = mopeds.Model(variable_list)

    e0_T_j2 = m.varlist_all["e0_T_j2"].casadi_var  # noqa: E501
    e0_A_r3_par4 = m.varlist_all["e0_A_r3_par4"].casadi_var  # noqa: E501
    e0_T_j1 = m.varlist_all["e0_T_j1"].casadi_var  # noqa: E501
    e0_A_r3_par5 = m.varlist_all["e0_A_r3_par5"].casadi_var  # noqa: E501
    e0_A_r3_par6 = m.varlist_all["e0_A_r3_par6"].casadi_var  # noqa: E501
    e0_x_i1_j1 = m.varlist_all["e0_x_i1_j1"].casadi_var  # noqa: E501
    e0_x_i2_j1 = m.varlist_all["e0_x_i2_j1"].casadi_var  # noqa: E501
    e0_x_i3_j1 = m.varlist_all["e0_x_i3_j1"].casadi_var  # noqa: E501
    e0_x_i4_j1 = m.varlist_all["e0_x_i4_j1"].casadi_var  # noqa: E501
    e0_x_i6_j1 = m.varlist_all["e0_x_i6_j1"].casadi_var  # noqa: E501
    e0_x_i7_j1 = m.varlist_all["e0_x_i7_j1"].casadi_var  # noqa: E501
    e0_p_j2 = m.varlist_all["e0_p_j2"].casadi_var  # noqa: E501
    e0_V = m.varlist_all["e0_V"].casadi_var  # noqa: E501
    e0_greek_rho_r3 = m.varlist_all["e0_greek_rho_r3"].casadi_var  # noqa: E501
    e0_A_r3_par1 = m.varlist_all["e0_A_r3_par1"].casadi_var  # noqa: E501
    e0_greek_rho_r1 = m.varlist_all["e0_greek_rho_r1"].casadi_var  # noqa: E501
    e0_F_j1 = m.varlist_all["e0_F_j1"].casadi_var  # noqa: E501
    e0_A_r3_par2 = m.varlist_all["e0_A_r3_par2"].casadi_var  # noqa: E501
    e0_A_r3_par3 = m.varlist_all["e0_A_r3_par3"].casadi_var  # noqa: E501
    e0_K_i3 = m.varlist_all["e0_K_i3"].casadi_var  # noqa: E501
    e0_K_i6 = m.varlist_all["e0_K_i6"].casadi_var  # noqa: E501
    e0_h_o_j1_i1 = m.varlist_all["e0_h_o_j1_i1"].casadi_var  # noqa: E501
    e0_h_o_j1_i2 = m.varlist_all["e0_h_o_j1_i2"].casadi_var  # noqa: E501
    e0_h_o_j1_i3 = m.varlist_all["e0_h_o_j1_i3"].casadi_var  # noqa: E501
    e0_h_o_j1_i4 = m.varlist_all["e0_h_o_j1_i4"].casadi_var  # noqa: E501
    e0_h_o_j1_i5 = m.varlist_all["e0_h_o_j1_i5"].casadi_var  # noqa: E501
    e0_h_o_j1_i6 = m.varlist_all["e0_h_o_j1_i6"].casadi_var  # noqa: E501
    e0_h_o_j1_i7 = m.varlist_all["e0_h_o_j1_i7"].casadi_var  # noqa: E501
    e0_h_o_j2_i1 = m.varlist_all["e0_h_o_j2_i1"].casadi_var  # noqa: E501
    e0_h_o_j2_i2 = m.varlist_all["e0_h_o_j2_i2"].casadi_var  # noqa: E501
    e0_h_o_j2_i3 = m.varlist_all["e0_h_o_j2_i3"].casadi_var  # noqa: E501
    e0_h_o_j2_i4 = m.varlist_all["e0_h_o_j2_i4"].casadi_var  # noqa: E501
    e0_h_o_j2_i5 = m.varlist_all["e0_h_o_j2_i5"].casadi_var  # noqa: E501
    e0_h_o_j2_i6 = m.varlist_all["e0_h_o_j2_i6"].casadi_var  # noqa: E501
    e0_h_o_j2_i7 = m.varlist_all["e0_h_o_j2_i7"].casadi_var  # noqa: E501
    e0_k_par1 = m.varlist_all["e0_k_par1"].casadi_var  # noqa: E501
    e0_k_par2 = m.varlist_all["e0_k_par2"].casadi_var  # noqa: E501
    e0_k_par3 = m.varlist_all["e0_k_par3"].casadi_var  # noqa: E501
    e0_k_par4 = m.varlist_all["e0_k_par4"].casadi_var  # noqa: E501
    e0_k_par5 = m.varlist_all["e0_k_par5"].casadi_var  # noqa: E501
    e0_k_r3 = m.varlist_all["e0_k_r3"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i1 = m.varlist_all["e0_greek_nu_r1_i1"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i2 = m.varlist_all["e0_greek_nu_r1_i2"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i3 = m.varlist_all["e0_greek_nu_r1_i3"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i4 = m.varlist_all["e0_greek_nu_r1_i4"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i5 = m.varlist_all["e0_greek_nu_r1_i5"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i6 = m.varlist_all["e0_greek_nu_r1_i6"].casadi_var  # noqa: E501
    e0_greek_nu_r1_i7 = m.varlist_all["e0_greek_nu_r1_i7"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i1 = m.varlist_all["e0_greek_nu_r2_i1"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i2 = m.varlist_all["e0_greek_nu_r2_i2"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i3 = m.varlist_all["e0_greek_nu_r2_i3"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i4 = m.varlist_all["e0_greek_nu_r2_i4"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i5 = m.varlist_all["e0_greek_nu_r2_i5"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i6 = m.varlist_all["e0_greek_nu_r2_i6"].casadi_var  # noqa: E501
    e0_greek_nu_r2_i7 = m.varlist_all["e0_greek_nu_r2_i7"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i1 = m.varlist_all["e0_greek_nu_r3_i1"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i2 = m.varlist_all["e0_greek_nu_r3_i2"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i3 = m.varlist_all["e0_greek_nu_r3_i3"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i4 = m.varlist_all["e0_greek_nu_r3_i4"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i5 = m.varlist_all["e0_greek_nu_r3_i5"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i6 = m.varlist_all["e0_greek_nu_r3_i6"].casadi_var  # noqa: E501
    e0_greek_nu_r3_i7 = m.varlist_all["e0_greek_nu_r3_i7"].casadi_var  # noqa: E501
    e0_A_par1_i1 = m.varlist_all["e0_A_par1_i1"].casadi_var  # noqa: E501
    e0_A_par1_i2 = m.varlist_all["e0_A_par1_i2"].casadi_var  # noqa: E501
    e0_A_par1_i3 = m.varlist_all["e0_A_par1_i3"].casadi_var  # noqa: E501
    e0_A_par1_i4 = m.varlist_all["e0_A_par1_i4"].casadi_var  # noqa: E501
    e0_A_par1_i5 = m.varlist_all["e0_A_par1_i5"].casadi_var  # noqa: E501
    e0_A_par1_i6 = m.varlist_all["e0_A_par1_i6"].casadi_var  # noqa: E501
    e0_A_par1_i7 = m.varlist_all["e0_A_par1_i7"].casadi_var  # noqa: E501
    e0_A_par1_r1 = m.varlist_all["e0_A_par1_r1"].casadi_var  # noqa: E501
    e0_A_par2_i1 = m.varlist_all["e0_A_par2_i1"].casadi_var  # noqa: E501
    e0_A_par2_i2 = m.varlist_all["e0_A_par2_i2"].casadi_var  # noqa: E501
    e0_A_par2_i3 = m.varlist_all["e0_A_par2_i3"].casadi_var  # noqa: E501
    e0_A_par2_i4 = m.varlist_all["e0_A_par2_i4"].casadi_var  # noqa: E501
    e0_A_par2_i5 = m.varlist_all["e0_A_par2_i5"].casadi_var  # noqa: E501
    e0_A_par2_i6 = m.varlist_all["e0_A_par2_i6"].casadi_var  # noqa: E501
    e0_A_par2_i7 = m.varlist_all["e0_A_par2_i7"].casadi_var  # noqa: E501
    e0_A_par2_r1 = m.varlist_all["e0_A_par2_r1"].casadi_var  # noqa: E501
    e0_A_par3_i1 = m.varlist_all["e0_A_par3_i1"].casadi_var  # noqa: E501
    e0_A_par3_i2 = m.varlist_all["e0_A_par3_i2"].casadi_var  # noqa: E501
    e0_A_par3_i3 = m.varlist_all["e0_A_par3_i3"].casadi_var  # noqa: E501
    e0_A_par3_i4 = m.varlist_all["e0_A_par3_i4"].casadi_var  # noqa: E501
    e0_A_par3_i5 = m.varlist_all["e0_A_par3_i5"].casadi_var  # noqa: E501
    e0_A_par3_i6 = m.varlist_all["e0_A_par3_i6"].casadi_var  # noqa: E501
    e0_A_par3_i7 = m.varlist_all["e0_A_par3_i7"].casadi_var  # noqa: E501
    e0_A_par3_r1 = m.varlist_all["e0_A_par3_r1"].casadi_var  # noqa: E501
    e0_A_par4_i1 = m.varlist_all["e0_A_par4_i1"].casadi_var  # noqa: E501
    e0_A_par4_i2 = m.varlist_all["e0_A_par4_i2"].casadi_var  # noqa: E501
    e0_A_par4_i3 = m.varlist_all["e0_A_par4_i3"].casadi_var  # noqa: E501
    e0_A_par4_i4 = m.varlist_all["e0_A_par4_i4"].casadi_var  # noqa: E501
    e0_A_par4_i5 = m.varlist_all["e0_A_par4_i5"].casadi_var  # noqa: E501
    e0_A_par4_i6 = m.varlist_all["e0_A_par4_i6"].casadi_var  # noqa: E501
    e0_A_par4_i7 = m.varlist_all["e0_A_par4_i7"].casadi_var  # noqa: E501
    e0_A_par4_r1 = m.varlist_all["e0_A_par4_r1"].casadi_var  # noqa: E501
    e0_A_par5_i1 = m.varlist_all["e0_A_par5_i1"].casadi_var  # noqa: E501
    e0_A_par5_i2 = m.varlist_all["e0_A_par5_i2"].casadi_var  # noqa: E501
    e0_A_par5_i3 = m.varlist_all["e0_A_par5_i3"].casadi_var  # noqa: E501
    e0_A_par5_i4 = m.varlist_all["e0_A_par5_i4"].casadi_var  # noqa: E501
    e0_A_par5_i5 = m.varlist_all["e0_A_par5_i5"].casadi_var  # noqa: E501
    e0_A_par5_i6 = m.varlist_all["e0_A_par5_i6"].casadi_var  # noqa: E501
    e0_A_par5_i7 = m.varlist_all["e0_A_par5_i7"].casadi_var  # noqa: E501
    e0_A_par5_r1 = m.varlist_all["e0_A_par5_r1"].casadi_var  # noqa: E501
    e0_A_eq_r1_par1 = m.varlist_all["e0_A_eq_r1_par1"].casadi_var  # noqa: E501
    e0_A_eq_r1_par2 = m.varlist_all["e0_A_eq_r1_par2"].casadi_var  # noqa: E501
    e0_A_eq_r2_par1 = m.varlist_all["e0_A_eq_r2_par1"].casadi_var  # noqa: E501
    e0_A_eq_r2_par2 = m.varlist_all["e0_A_eq_r2_par2"].casadi_var  # noqa: E501
    e0_A_r3_i3 = m.varlist_all["e0_A_r3_i3"].casadi_var  # noqa: E501
    e0_A_r3_i6 = m.varlist_all["e0_A_r3_i6"].casadi_var  # noqa: E501
    e0_A_r3 = m.varlist_all["e0_A_r3"].casadi_var  # noqa: E501
    e0_E_par1_r1 = m.varlist_all["e0_E_par1_r1"].casadi_var  # noqa: E501
    e0_E_par2_r1 = m.varlist_all["e0_E_par2_r1"].casadi_var  # noqa: E501
    e0_E_par3_r1 = m.varlist_all["e0_E_par3_r1"].casadi_var  # noqa: E501
    e0_E_par4_r1 = m.varlist_all["e0_E_par4_r1"].casadi_var  # noqa: E501
    e0_E_par5_r1 = m.varlist_all["e0_E_par5_r1"].casadi_var  # noqa: E501
    e0_E_r3_i3 = m.varlist_all["e0_E_r3_i3"].casadi_var  # noqa: E501
    e0_E_r3_i6 = m.varlist_all["e0_E_r3_i6"].casadi_var  # noqa: E501
    e0_E_r3 = m.varlist_all["e0_E_r3"].casadi_var  # noqa: E501
    e0_R = m.varlist_all["e0_R"].casadi_var  # noqa: E501
    e0_T_f = m.varlist_all["e0_T_f"].casadi_var  # noqa: E501
    e0_h_f_i1 = m.varlist_all["e0_h_f_i1"].casadi_var  # noqa: E501
    e0_h_f_i2 = m.varlist_all["e0_h_f_i2"].casadi_var  # noqa: E501
    e0_h_f_i3 = m.varlist_all["e0_h_f_i3"].casadi_var  # noqa: E501
    e0_h_f_i4 = m.varlist_all["e0_h_f_i4"].casadi_var  # noqa: E501
    e0_h_f_i5 = m.varlist_all["e0_h_f_i5"].casadi_var  # noqa: E501
    e0_h_f_i6 = m.varlist_all["e0_h_f_i6"].casadi_var  # noqa: E501
    e0_h_f_i7 = m.varlist_all["e0_h_f_i7"].casadi_var  # noqa: E501
    e0_K_r2 = m.varlist_all["e0_K_r2"].casadi_var  # noqa: E501
    e0_X_i5 = m.varlist_all["e0_X_i5"].casadi_var  # noqa: E501
    e0_X_i2 = m.varlist_all["e0_X_i2"].casadi_var  # noqa: E501
    e0_Y_i7 = m.varlist_all["e0_Y_i7"].casadi_var  # noqa: E501
    e0_K_r3 = m.varlist_all["e0_K_r3"].casadi_var  # noqa: E501
    e0_x_i1_j2 = m.varlist_all["e0_x_i1_j2"].casadi_var  # noqa: E501
    e0_x_i2_j2 = m.varlist_all["e0_x_i2_j2"].casadi_var  # noqa: E501
    e0_x_i3_j2 = m.varlist_all["e0_x_i3_j2"].casadi_var  # noqa: E501
    e0_x_i4_j2 = m.varlist_all["e0_x_i4_j2"].casadi_var  # noqa: E501
    e0_x_i5_j1 = m.varlist_all["e0_x_i5_j1"].casadi_var  # noqa: E501
    e0_x_i5_j2 = m.varlist_all["e0_x_i5_j2"].casadi_var  # noqa: E501
    e0_x_i6_j2 = m.varlist_all["e0_x_i6_j2"].casadi_var  # noqa: E501
    e0_x_i7_j2 = m.varlist_all["e0_x_i7_j2"].casadi_var  # noqa: E501
    e0_h_j1 = m.varlist_all["e0_h_j1"].casadi_var  # noqa: E501
    e0_h_j2 = m.varlist_all["e0_h_j2"].casadi_var  # noqa: E501
    e0_K_r1 = m.varlist_all["e0_K_r1"].casadi_var  # noqa: E501
    e0_p_i1_j2 = m.varlist_all["e0_p_i1_j2"].casadi_var  # noqa: E501
    e0_p_i2_j2 = m.varlist_all["e0_p_i2_j2"].casadi_var  # noqa: E501
    e0_p_i3_j2 = m.varlist_all["e0_p_i3_j2"].casadi_var  # noqa: E501
    e0_p_i4_j2 = m.varlist_all["e0_p_i4_j2"].casadi_var  # noqa: E501
    e0_p_i5_j2 = m.varlist_all["e0_p_i5_j2"].casadi_var  # noqa: E501
    e0_p_i6_j2 = m.varlist_all["e0_p_i6_j2"].casadi_var  # noqa: E501
    e0_p_i7_j2 = m.varlist_all["e0_p_i7_j2"].casadi_var  # noqa: E501
    e0_HU = m.varlist_all["e0_HU"].casadi_var  # noqa: E501
    e0_U = m.varlist_all["e0_U"].casadi_var  # noqa: E501
    e0_u_j2 = m.varlist_all["e0_u_j2"].casadi_var  # noqa: E501
    e0_HU_i1 = m.varlist_all["e0_HU_i1"].casadi_var  # noqa: E501
    e0_HU_i2 = m.varlist_all["e0_HU_i2"].casadi_var  # noqa: E501
    e0_HU_i3 = m.varlist_all["e0_HU_i3"].casadi_var  # noqa: E501
    e0_HU_i4 = m.varlist_all["e0_HU_i4"].casadi_var  # noqa: E501
    e0_HU_i5 = m.varlist_all["e0_HU_i5"].casadi_var  # noqa: E501
    e0_HU_i6 = m.varlist_all["e0_HU_i6"].casadi_var  # noqa: E501
    e0_HU_i7 = m.varlist_all["e0_HU_i7"].casadi_var  # noqa: E501
    e0_v_j2 = m.varlist_all["e0_v_j2"].casadi_var  # noqa: E501
    e0_r_r3 = m.varlist_all["e0_r_r3"].casadi_var  # noqa: E501
    e0_r_r2 = m.varlist_all["e0_r_r2"].casadi_var  # noqa: E501
    e0_r_r1 = m.varlist_all["e0_r_r1"].casadi_var  # noqa: E501
    e0_F_j2 = m.varlist_all["e0_F_j2"].casadi_var  # noqa: E501
    e0_Q = m.varlist_all["e0_Q"].casadi_var  # noqa: E501

    EQ_alg1 = (e0_K_r2-(((10.0))**(1.0*((e0_A_eq_r2_par1/e0_T_j2)+e0_A_eq_r2_par2))))  # noqa: E501,E226
    EQ_alg2 = (e0_K_r1-(((10.0))**(1.0*((e0_A_eq_r1_par1/e0_T_j2)+e0_A_eq_r1_par2))))  # noqa: E501,E226
    EQ_alg3 = (e0_K_r3-(ca.exp(((((((e0_A_r3_par1/e0_T_j2)+(e0_A_r3_par2*ca.log(e0_T_j2)))+(e0_A_r3_par3*e0_T_j2))+(e0_A_r3_par4*((e0_T_j2))**(1.0*2.0)))+(e0_A_r3_par5/((e0_T_j2))**(1.0*3.0)))+e0_A_r3_par6))))  # noqa: E501,E226
    EQ_alg4 = (1.0-((((((((e0_x_i1_j1+e0_x_i2_j1)+e0_x_i3_j1)+e0_x_i4_j1)+e0_x_i5_j1)+e0_x_i6_j1)+e0_x_i7_j1))))  # noqa: E501,E226
    EQ_alg5 = (1.0-((((((((e0_x_i1_j2+e0_x_i2_j2)+e0_x_i3_j2)+e0_x_i4_j2)+e0_x_i5_j2)+e0_x_i6_j2)+e0_x_i7_j2))))  # noqa: E501,E226
    EQ_alg6 = (e0_h_j1-((((((((((e0_x_i1_j1*e0_h_o_j1_i1)+(e0_x_i2_j1*e0_h_o_j1_i2))+(e0_x_i3_j1*e0_h_o_j1_i3))+(e0_x_i4_j1*e0_h_o_j1_i4))+(e0_x_i5_j1*e0_h_o_j1_i5))+(e0_x_i6_j1*e0_h_o_j1_i6))+(e0_x_i7_j1*e0_h_o_j1_i7)))/1000000.0)))  # noqa: E501,E226
    EQ_alg7 = (e0_h_j2-((((((((((e0_x_i1_j2*e0_h_o_j2_i1)+(e0_x_i2_j2*e0_h_o_j2_i2))+(e0_x_i3_j2*e0_h_o_j2_i3))+(e0_x_i4_j2*e0_h_o_j2_i4))+(e0_x_i5_j2*e0_h_o_j2_i5))+(e0_x_i6_j2*e0_h_o_j2_i6))+(e0_x_i7_j2*e0_h_o_j2_i7)))/1000000.0)))  # noqa: E501,E226
    EQ_alg8 = (e0_p_i1_j2-((e0_p_j2*e0_x_i1_j2)))  # noqa: E501,E226
    EQ_alg9 = (e0_p_i2_j2-((e0_p_j2*e0_x_i2_j2)))  # noqa: E501,E226
    EQ_alg10 = (e0_p_i3_j2-((e0_p_j2*e0_x_i3_j2)))  # noqa: E501,E226
    EQ_alg11 = (e0_p_i4_j2-((e0_p_j2*e0_x_i4_j2)))  # noqa: E501,E226
    EQ_alg12 = (e0_p_i5_j2-((e0_p_j2*e0_x_i5_j2)))  # noqa: E501,E226
    EQ_alg13 = (e0_p_i6_j2-((e0_p_j2*e0_x_i6_j2)))  # noqa: E501,E226
    EQ_alg14 = (e0_p_i7_j2-((e0_p_j2*e0_x_i7_j2)))  # noqa: E501,E226
    EQ_alg15 = (e0_U-((e0_u_j2*(e0_HU*1000.0))))  # noqa: E501,E226
    EQ_alg16 = ((e0_p_j2*(100.0*e0_V))-((e0_HU*(e0_R*e0_T_j2))))  # noqa: E501,E226
    EQ_alg17 = (e0_HU_i1-((e0_HU*e0_x_i1_j2)))  # noqa: E501,E226
    EQ_alg18 = (e0_HU_i2-((e0_HU*e0_x_i2_j2)))  # noqa: E501,E226
    EQ_alg19 = (e0_HU_i3-((e0_HU*e0_x_i3_j2)))  # noqa: E501,E226
    EQ_alg20 = (e0_HU_i4-((e0_HU*e0_x_i4_j2)))  # noqa: E501,E226
    EQ_alg21 = (e0_HU_i5-((e0_HU*e0_x_i5_j2)))  # noqa: E501,E226
    EQ_alg22 = (e0_HU_i6-((e0_HU*e0_x_i6_j2)))  # noqa: E501,E226
    EQ_alg23 = (e0_HU_i7-((e0_HU*e0_x_i7_j2)))  # noqa: E501,E226
    EQ_alg24 = (e0_v_j2-((e0_V/(e0_HU*1000.0))))  # noqa: E501,E226
    EQ_alg25 = (e0_r_r3-((e0_greek_rho_r3*(e0_k_r3*(((e0_K_i6))**(1.0*2.0)*(((((e0_p_i6_j2))**(1.0*2.0)-((e0_p_i7_j2*e0_p_i3_j2)/e0_K_r3))/((((1.0+(e0_K_i6*e0_p_i6_j2))+(e0_K_i3*e0_p_i3_j2))))**(1.0*2.0))*e0_V))))))  # noqa: E501,E226
    EQ_alg26 = (e0_r_r2-((e0_greek_rho_r1*(e0_k_par1*(((e0_p_i5_j2*(e0_p_i2_j2*((1.0-((e0_p_i3_j2*e0_p_i6_j2)/(((e0_p_i2_j2))**(1.0*3.0)*(e0_p_i5_j2*e0_K_r2)))))))/(((((1.0+(e0_k_par2*(e0_p_i3_j2/e0_p_i2_j2)))+(e0_k_par3*((e0_p_i2_j2))**(1.0*0.5)))+(e0_k_par4*(e0_p_i3_j2)))))**(1.0*3.0))*e0_V)))))  # noqa: E501,E226
    EQ_alg27 = (e0_r_r1-((e0_greek_rho_r1*(e0_k_par5*(((e0_p_i5_j2*((1.0-((e0_p_i4_j2*e0_p_i3_j2)/(e0_p_i2_j2*(e0_p_i5_j2*e0_K_r1))))))/(((1.0+(e0_k_par2*(e0_p_i3_j2/e0_p_i2_j2)))+(e0_k_par3*((e0_p_i2_j2))**(1.0*0.5)))+(e0_k_par4*(e0_p_i3_j2))))*e0_V)))))  # noqa: E501,E226
    EQ_alg28 = ((e0_u_j2*10.0)-(((e0_h_j2*10.0)-(e0_p_j2*e0_v_j2))))  # noqa: E501,E226
    EQ_alg29 = (0.0-(((((e0_F_j1*(1000.0*e0_h_j1))+(e0_F_j2*(1000.0*e0_h_j2))))+e0_Q)))  # noqa: E501,E226
    EQ_alg30 = (0.0-(((((e0_F_j1*e0_x_i1_j1)+(e0_F_j2*e0_x_i1_j2)))+((((e0_greek_nu_r1_i1*e0_r_r1)+(e0_greek_nu_r2_i1*e0_r_r2))+(e0_greek_nu_r3_i1*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg31 = (0.0-(((((e0_F_j1*e0_x_i2_j1)+(e0_F_j2*e0_x_i2_j2)))+((((e0_greek_nu_r1_i2*e0_r_r1)+(e0_greek_nu_r2_i2*e0_r_r2))+(e0_greek_nu_r3_i2*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg32 = (0.0-(((((e0_F_j1*e0_x_i3_j1)+(e0_F_j2*e0_x_i3_j2)))+((((e0_greek_nu_r1_i3*e0_r_r1)+(e0_greek_nu_r2_i3*e0_r_r2))+(e0_greek_nu_r3_i3*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg33 = (0.0-(((((e0_F_j1*e0_x_i4_j1)+(e0_F_j2*e0_x_i4_j2)))+((((e0_greek_nu_r1_i4*e0_r_r1)+(e0_greek_nu_r2_i4*e0_r_r2))+(e0_greek_nu_r3_i4*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg34 = (0.0-(((((e0_F_j1*e0_x_i5_j1)+(e0_F_j2*e0_x_i5_j2)))+((((e0_greek_nu_r1_i5*e0_r_r1)+(e0_greek_nu_r2_i5*e0_r_r2))+(e0_greek_nu_r3_i5*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg35 = (0.0-(((((e0_F_j1*e0_x_i6_j1)+(e0_F_j2*e0_x_i6_j2)))+((((e0_greek_nu_r1_i6*e0_r_r1)+(e0_greek_nu_r2_i6*e0_r_r2))+(e0_greek_nu_r3_i6*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg36 = (0.0-(((((e0_F_j1*e0_x_i7_j1)+(e0_F_j2*e0_x_i7_j2)))+((((e0_greek_nu_r1_i7*e0_r_r1)+(e0_greek_nu_r2_i7*e0_r_r2))+(e0_greek_nu_r3_i7*e0_r_r3))))))  # noqa: E501,E226
    EQ_alg37 = (e0_X_i5-((((e0_F_j1*e0_x_i5_j1)+(e0_F_j2*e0_x_i5_j2))/(e0_F_j1*e0_x_i5_j1))))  # noqa: E501,E226
    EQ_alg38 = (e0_X_i2-((((e0_F_j1*e0_x_i2_j1)+(e0_F_j2*e0_x_i2_j2))/(e0_F_j1*e0_x_i2_j1))))  # noqa: E501,E226
    EQ_alg39 = (e0_Y_i7-((-(1.0*((((e0_F_j1*e0_x_i7_j1)+(e0_F_j2*e0_x_i7_j2))/(e0_F_j1*e0_x_i5_j1))*(1.0/2.0))))))  # noqa: E501,E226

    list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6, EQ_alg7, EQ_alg8, EQ_alg9, EQ_alg10, EQ_alg11, EQ_alg12, EQ_alg13, EQ_alg14, EQ_alg15, EQ_alg16, EQ_alg17, EQ_alg18, EQ_alg19, EQ_alg20, EQ_alg21, EQ_alg22, EQ_alg23, EQ_alg24, EQ_alg25, EQ_alg26, EQ_alg27, EQ_alg28, EQ_alg29, EQ_alg30, EQ_alg31, EQ_alg32, EQ_alg33, EQ_alg34, EQ_alg35, EQ_alg36, EQ_alg37, EQ_alg38, EQ_alg39, ]  # noqa: E501
    try:
        Eq_fun_e0_h_o_j2_i4 = m.varlist_all["e0_h_o_j2_i4"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i4,e0_A_par2_i4,e0_A_par3_i4,e0_A_par4_i4,e0_A_par5_i4,e0_T_f,e0_h_f_i4)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i4)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_K_i3 = m.varlist_all["e0_K_i3"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_r3_i3,e0_E_r3_i3,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_K_i3)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_K_i6 = m.varlist_all["e0_K_i6"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_r3_i6,e0_E_r3_i6,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_K_i6)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i1 = m.varlist_all["e0_h_o_j1_i1"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i1,e0_A_par2_i1,e0_A_par3_i1,e0_A_par4_i1,e0_A_par5_i1,e0_T_f,e0_h_f_i1)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i1)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j2_i1 = m.varlist_all["e0_h_o_j2_i1"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i1,e0_A_par2_i1,e0_A_par3_i1,e0_A_par4_i1,e0_A_par5_i1,e0_T_f,e0_h_f_i1)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i1)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i4 = m.varlist_all["e0_h_o_j1_i4"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i4,e0_A_par2_i4,e0_A_par3_i4,e0_A_par4_i4,e0_A_par5_i4,e0_T_f,e0_h_f_i4)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i4)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j2_i5 = m.varlist_all["e0_h_o_j2_i5"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i5,e0_A_par2_i5,e0_A_par3_i5,e0_A_par4_i5,e0_A_par5_i5,e0_T_f,e0_h_f_i5)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i5)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_par4 = m.varlist_all["e0_k_par4"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_par4_r1,e0_E_par4_r1,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_par4)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_par5 = m.varlist_all["e0_k_par5"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_par5_r1,e0_E_par5_r1,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_par5)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i5 = m.varlist_all["e0_h_o_j1_i5"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i5,e0_A_par2_i5,e0_A_par3_i5,e0_A_par4_i5,e0_A_par5_i5,e0_T_f,e0_h_f_i5)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i5)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j2_i6 = m.varlist_all["e0_h_o_j2_i6"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i6,e0_A_par2_i6,e0_A_par3_i6,e0_A_par4_i6,e0_A_par5_i6,e0_T_f,e0_h_f_i6)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i6)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j2_i7 = m.varlist_all["e0_h_o_j2_i7"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i7,e0_A_par2_i7,e0_A_par3_i7,e0_A_par4_i7,e0_A_par5_i7,e0_T_f,e0_h_f_i7)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i7)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i2 = m.varlist_all["e0_h_o_j1_i2"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i2,e0_A_par2_i2,e0_A_par3_i2,e0_A_par4_i2,e0_A_par5_i2,e0_T_f,e0_h_f_i2)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i2)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j2_i3 = m.varlist_all["e0_h_o_j2_i3"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i3,e0_A_par2_i3,e0_A_par3_i3,e0_A_par4_i3,e0_A_par5_i3,e0_T_f,e0_h_f_i3)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i3)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_par3 = m.varlist_all["e0_k_par3"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_par3_r1,e0_E_par3_r1,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_par3)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i6 = m.varlist_all["e0_h_o_j1_i6"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i6,e0_A_par2_i6,e0_A_par3_i6,e0_A_par4_i6,e0_A_par5_i6,e0_T_f,e0_h_f_i6)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i6)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_par1 = m.varlist_all["e0_k_par1"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_par1_r1,e0_E_par1_r1,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_par1)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j2_i2 = m.varlist_all["e0_h_o_j2_i2"].casadi_var - fun_113237__enthalpyFNC(e0_T_j2,e0_A_par1_i2,e0_A_par2_i2,e0_A_par3_i2,e0_A_par4_i2,e0_A_par5_i2,e0_T_f,e0_h_f_i2)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j2_i2)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i3 = m.varlist_all["e0_h_o_j1_i3"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i3,e0_A_par2_i3,e0_A_par3_i3,e0_A_par4_i3,e0_A_par5_i3,e0_T_f,e0_h_f_i3)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i3)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_r3 = m.varlist_all["e0_k_r3"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_r3,e0_E_r3,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_r3)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_h_o_j1_i7 = m.varlist_all["e0_h_o_j1_i7"].casadi_var - fun_113237__enthalpyFNC(e0_T_j1,e0_A_par1_i7,e0_A_par2_i7,e0_A_par3_i7,e0_A_par4_i7,e0_A_par5_i7,e0_T_f,e0_h_f_i7)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_h_o_j1_i7)  # noqa: E501
    except KeyError:
        pass
    try:
        Eq_fun_e0_k_par2 = m.varlist_all["e0_k_par2"].casadi_var - fun_113291__ArrheniusFNC(e0_T_j2,e0_A_par2_r1,e0_E_par2_r1,e0_R)  # noqa: E501,E231
        list_algebraic_equations.append(Eq_fun_e0_k_par2)  # noqa: E501
    except KeyError:
        pass

    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


if __name__ == "__main__":

    # 1. Intialize model with corresponding variable_list
    variable_list, m = initialize_problem()

    # 1.1 Define which parameters are to be unfixed and optimized here. First fix all parametere
    for var in variable_list.values():
        var.fixed = True

    # 1.2 Unfix them here
    variable_list["e0_A_r3_i3"].fixed = False
    variable_list["e0_E_par2_r1"].fixed = False
    variable_list["e0_A_par5_r1"].fixed = False
    variable_list["e0_A_par4_r1"].fixed = False
    variable_list["e0_E_par4_r1"].fixed = False

    # 1.3 Set some consistent bounds to parameters (either in a model when you generate variables) or just set bounds to be +- 50% from the guess. Also sett guess to be not the real value, but somewhere on lower bound (else your optimizer guess is a local optima)
    variable_list.set_bounds(0.5, emerg_val=1)

    # 2. Create a new variable list that would be used to "generate" experimental data and set all Parameter Variables to fixed.
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # 3. In this list change set values of control variables and parameter variables as you want
    var_list_fixed["e0_T_j2"].value = 500

    # 4. Create simulation Object that would generate experimental data
    sim_fixed = mopeds.SimulatorNLE(m, var_list_fixed, solver_name="ipopt")

    """
    5. Generate data - what it does: 
    - Run Simulator to solve NLE at given values of Controls and Parameters
    - Create a new VariableList with all Algebraic Variables
    - Set results of the simulation to var.value of the Variable Objects
    (So you can see the simulation results by printng print(res.values())
    """
    res = sim_fixed.generate_exp_data()

    # 6. Now create a variable list that would tell PE that there is "experimental data to be optimized". It's done by setting algebraic variables in variable_optimizer_dictionary to variable in res dictionary
    variable_list_optimizer = copy.deepcopy(variable_list)
    for key, var in res.items():
        # You also have to supply guess for NLE for algebraic Variables. You eiter set it to var.value.value (to real answer), which makes everything easy for solver
        var.guess = var.value[0]

        # Or Set it to initial guess that was imported from MosaicModeling
        # var.guess = variable_list[var.name].guess

        # Here you set the variable
        variable_list_optimizer[key] = var

    """
    7. So you changed COntroller variables in Step3 only in fixed variable list, now you need last things:
    - Set values of controller variables to the values that you defined in step 3
    """
    for var in variable_list_optimizer.values():
        if isinstance(var, mopeds.VariableControl):
            var.value = var_list_fixed[var.name].value[0]

    # 7. This variable_list_optimizer conssists now of one experiment. If you want to create another experiment, repeat steps 3 till 7: example below. Of course you can create a function to repeat steps 4-7 for you
    var_list_fixed["e0_T_j2"].value = 500
    sim_fixed = mopeds.SimulatorNLE(m, var_list_fixed, solver_name="ipopt")
    res = sim_fixed.generate_exp_data()
    variable_list_optimizer_2 = copy.deepcopy(variable_list)
    for key, var in res.items():
        var.guess = var.value[0]
        variable_list_optimizer_2[key] = var
    for var in variable_list_optimizer_2.values():
        if isinstance(var, mopeds.VariableControl):
            var.value = var_list_fixed[var.name].value[0]

    # 8. Now just supply all your experiments as a list to ptimizer
    pe = mopeds.ParameterEstimationNLE(
        m, [variable_list_optimizer, variable_list_optimizer_2], simulator_name="ipopt"
    )
    res = pe.optimize(True)
