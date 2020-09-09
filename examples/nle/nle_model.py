import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est
from par_est.simulation import SimulatorNLE
from par_est.optimization import ParameterEstimationNLE

def initialize_problem():

    variable_list = par_est.VariableList()

    # fmt: off 
    variable_list.add_variable(par_est.VariableState("e0_x_i2_j2", 0.52397))
    variable_list.add_variable(par_est.VariableState("e0_x_i1_j2", 0.21455))
    variable_list.add_variable(par_est.VariableState("e0_HU_i3", 0.0074497))
    variable_list.add_variable(par_est.VariableState("e0_HU_i2", 0.063022))
    variable_list.add_variable(par_est.VariableState("e0_K_r3", 23.0699))
    variable_list.add_variable(par_est.VariableState("e0_HU_i1", 0.025805))
    variable_list.add_variable(par_est.VariableState("e0_u_j2", -0.083762))
    variable_list.add_variable(par_est.VariableState("e0_U", -10.0748))
    variable_list.add_variable(par_est.VariableState("e0_HU", 0.12028))
    variable_list.add_variable(par_est.VariableState("e0_p_i7_j2", 0.7498))
    variable_list.add_variable(par_est.VariableState("e0_p_i6_j2", 0.31874))
    variable_list.add_variable(par_est.VariableState("e0_p_i5_j2", 8.3803))
    variable_list.add_variable(par_est.VariableState("e0_p_i4_j2", 0.52869))
    variable_list.add_variable(par_est.VariableState("e0_p_i3_j2", 3.0968))
    variable_list.add_variable(par_est.VariableState("e0_p_i2_j2", 26.1983))
    variable_list.add_variable(par_est.VariableState("e0_p_i1_j2", 10.7273))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i7", -168023.1327725076))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i7", -184427.01699715087))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i6", -190537.8047525016))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i6", -201159.24396779988))
    variable_list.add_variable(par_est.VariableState("e0_Q", 4.6492))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i5", -385214.2159451514))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i5", -393705.0747228285))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i4", -104597.15130755084))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i4", -110675.6190793542))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i3", -234887.44250385024))
    variable_list.add_variable(par_est.VariableState("e0_F_j2", -0.9322))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i3", -241998.63319940824))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i2", 5877.734015746798))
    variable_list.add_variable(par_est.VariableState("e0_r_r1", 0.0098568))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i2", -144.29466839718256))
    variable_list.add_variable(par_est.VariableState("e0_k_par5", 0.6905973165554351))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j2_i1", 5910.160589771112))
    variable_list.add_variable(par_est.VariableState("e0_r_r2", 0.033901))
    variable_list.add_variable(par_est.VariableState("e0_h_o_j1_i1", -145.74383548430316))
    variable_list.add_variable(par_est.VariableState("e0_k_par4", 613.6897535075899))
    variable_list.add_variable(par_est.VariableState("e0_k_par3", 31.24090746419173))
    variable_list.add_variable(par_est.VariableState("e0_k_par2", 3453.38))
    variable_list.add_variable(par_est.VariableState("e0_k_par1", 16154.96060463936))
    variable_list.add_variable(par_est.VariableState("e0_h_j2", -0.079605))
    variable_list.add_variable(par_est.VariableState("e0_K_r1", 0.0076384))
    variable_list.add_variable(par_est.VariableState("e0_h_j1", -0.078857))
    variable_list.add_variable(par_est.VariableState("e0_x_i7_j2", 0.014996))
    variable_list.add_variable(par_est.VariableState("e0_x_i7_j1", 4.2463E-17))
    variable_list.add_variable(par_est.VariableState("e0_x_i6_j2", 0.0063748))
    variable_list.add_variable(par_est.VariableState("e0_r_r3", 0.013979))
    variable_list.add_variable(par_est.VariableState("e0_K_r2", 3.4674E-5))
    variable_list.add_variable(par_est.VariableState("e0_k_r3", 0.15069620388555666))
    variable_list.add_variable(par_est.VariableState("e0_x_i5_j2", 0.16761))
    variable_list.add_variable(par_est.VariableState("e0_K_i6", 2.132122117309319E13))
    variable_list.add_variable(par_est.VariableState("e0_K_i3", 2.2476341621984563E9))
    variable_list.add_variable(par_est.VariableState("e0_x_i4_j2", 0.010574))
    variable_list.add_variable(par_est.VariableState("e0_v_j2", 8.314E-4))
    variable_list.add_variable(par_est.VariableState("e0_HU_i7", 0.0018037))
    variable_list.add_variable(par_est.VariableState("e0_HU_i6", 7.6676E-4))
    variable_list.add_variable(par_est.VariableState("e0_HU_i5", 0.02016))
    variable_list.add_variable(par_est.VariableState("e0_HU_i4", 0.0012718))
    variable_list.add_variable(par_est.VariableState("e0_X_i2", 0.18593))
    variable_list.add_variable(par_est.VariableState("e0_Y_i7", 0.069896))
    variable_list.add_variable(par_est.VariableState("e0_x_i3_j2", 0.061937))
    variable_list.add_variable(par_est.VariableState("e0_X_i5", 0.21879))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_par4", 3.8E-7))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_par5", -65610.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_par6", -26.64))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_x_i1_j1", 0.2))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i2", -1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i2", -3.0))
    variable_list.add_variable(par_est.VariableParameter("e0_x_i2_j1", 0.6))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i2", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_x_i3_j1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i1", 3.7E-5))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i1", -0.0169))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i7", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i7", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i1", 5.45E-12))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i1", -2.44E-8))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i4", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i4", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i4", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i5", -1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i5", -1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i5", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i6", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i6", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i3", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r2_i3", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i3", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i3", -2.50993E-7))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i3", 2.2E-4))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_par2", 3.707))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_par3", -0.002783))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_par1", 4019.0))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i3", -241830.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i4", -1.54836E-11))
    variable_list.add_variable(par_est.VariableParameter("e0_R", 8.314))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i3", -0.074))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i3", 42.04061276))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i4", -0.00322))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i4", 29.59614774))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i4", 2.10847E-8))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i4", 6.07E-7))
    variable_list.add_variable(par_est.VariableParameter("e0_x_i4_j1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_x_i5_j1", 0.2))
    variable_list.add_variable(par_est.VariableParameter("e0_x_i6_j1", 1.0E-16))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i1", 31.5))
    variable_list.add_variable(par_est.VariableParameter("e0_T_j2", 500.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i3", 1.07368E-10))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i2", -7.52E-12))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i2", 2.7E-8))
    variable_list.add_variable(par_est.VariableParameter("e0_T_f", 298.15))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i2", 26.2))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i2", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i2", -3.17E-5))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i2", 0.0162))
    variable_list.add_variable(par_est.VariableParameter("e0_A_eq_r1_par1", -2073.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_eq_r1_par2", 2.029))
    variable_list.add_variable(par_est.VariableParameter("e0_V", 0.1))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i6", -200940.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_eq_r2_par1", 3066.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i7", 3.87561E-10))
    variable_list.add_variable(par_est.VariableParameter("e0_A_eq_r2_par2", -10.592))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i6", -0.121))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i6", 45.78156724))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i7", -0.0526))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i7", 45.14491836))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i7", -8.56412E-7))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i7", 6.29E-4))
    variable_list.add_variable(par_est.VariableParameter("e0_p_j2", 50.0))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i7", -184100.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_r1", 1.07))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i5", 1.93148E-11))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i5", -1.32293E-8))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i4", -110530.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_i5", 19.6353642))
    variable_list.add_variable(par_est.VariableParameter("e0_h_f_i5", -393520.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i5", -4.18E-5))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_i5", 0.0718))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_i6", 5.57E-4))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par1_i6", 2.65055E-10))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_i6", -6.5285E-7))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r3_i3", 92000.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_i6", 223.2))
    variable_list.add_variable(par_est.VariableParameter("e0_F_j1", 1.0))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r3", -55060.0))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r3_i6", 105100.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3", 85190.0))
    variable_list.add_variable(par_est.VariableParameter("e0_T_j1", 293.15))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r3_i6", -2.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_nu_r1_i7", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_rho_r1", 1775.0))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_rho_r3", 100.0))
    variable_list.add_variable(par_est.VariableParameter("e0_E_par1_r1", 40000.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par3_r1", 0.499))
    variable_list.add_variable(par_est.VariableParameter("e0_E_par3_r1", 17197.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par2_r1", 3453.38))
    variable_list.add_variable(par_est.VariableParameter("e0_E_par2_r1", 0.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par5_r1", 1.22E10))
    variable_list.add_variable(par_est.VariableParameter("e0_E_par5_r1", -98084.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_par4_r1", 6.62E-11))
    variable_list.add_variable(par_est.VariableParameter("e0_E_par4_r1", 124119.0))
    variable_list.add_variable(par_est.VariableParameter("e0_A_r3_i3", 0.5498))

    m = par_est.Model(variable_list)

    dydx1 = 	 m.varlist_all['e0_K_r2'].casadi_var - ((10.0) ** ((m.varlist_all['e0_A_eq_r2_par1'].casadi_var)/(m.varlist_all['e0_T_j2'].casadi_var) + m.varlist_all['e0_A_eq_r2_par2'].casadi_var)) 

    dydx2 = 	 m.varlist_all['e0_K_r1'].casadi_var - ((10.0) ** ((m.varlist_all['e0_A_eq_r1_par1'].casadi_var)/(m.varlist_all['e0_T_j2'].casadi_var) + m.varlist_all['e0_A_eq_r1_par2'].casadi_var)) 

    dydx3 = 	 m.varlist_all['e0_K_r3'].casadi_var - (ca.exp((m.varlist_all['e0_A_r3_par1'].casadi_var)/(m.varlist_all['e0_T_j2'].casadi_var) + m.varlist_all['e0_A_r3_par2'].casadi_var * ca.log(m.varlist_all['e0_T_j2'].casadi_var) + m.varlist_all['e0_A_r3_par3'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var + m.varlist_all['e0_A_r3_par4'].casadi_var * (m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) + (m.varlist_all['e0_A_r3_par5'].casadi_var)/((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0)) + m.varlist_all['e0_A_r3_par6'].casadi_var)) 

    dydx4 = 	 1.0 - ((m.varlist_all['e0_x_i1_j1'].casadi_var + m.varlist_all['e0_x_i2_j1'].casadi_var + m.varlist_all['e0_x_i3_j1'].casadi_var + m.varlist_all['e0_x_i4_j1'].casadi_var + m.varlist_all['e0_x_i5_j1'].casadi_var + m.varlist_all['e0_x_i6_j1'].casadi_var + m.varlist_all['e0_x_i7_j1'].casadi_var)) 

    dydx5 = 	 1.0 - ((m.varlist_all['e0_x_i1_j2'].casadi_var + m.varlist_all['e0_x_i2_j2'].casadi_var + m.varlist_all['e0_x_i3_j2'].casadi_var + m.varlist_all['e0_x_i4_j2'].casadi_var + m.varlist_all['e0_x_i5_j2'].casadi_var + m.varlist_all['e0_x_i6_j2'].casadi_var + m.varlist_all['e0_x_i7_j2'].casadi_var)) 

    dydx6 = 	 m.varlist_all['e0_h_j1'].casadi_var - (((m.varlist_all['e0_x_i1_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i1'].casadi_var + m.varlist_all['e0_x_i2_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i2'].casadi_var + m.varlist_all['e0_x_i3_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i3'].casadi_var + m.varlist_all['e0_x_i4_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i4'].casadi_var + m.varlist_all['e0_x_i5_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i5'].casadi_var + m.varlist_all['e0_x_i6_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i6'].casadi_var + m.varlist_all['e0_x_i7_j1'].casadi_var * m.varlist_all['e0_h_o_j1_i7'].casadi_var))/(1000000.0)) 

    dydx7 = 	 m.varlist_all['e0_h_j2'].casadi_var - (((m.varlist_all['e0_x_i1_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i1'].casadi_var + m.varlist_all['e0_x_i2_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i2'].casadi_var + m.varlist_all['e0_x_i3_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i3'].casadi_var + m.varlist_all['e0_x_i4_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i4'].casadi_var + m.varlist_all['e0_x_i5_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i5'].casadi_var + m.varlist_all['e0_x_i6_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i6'].casadi_var + m.varlist_all['e0_x_i7_j2'].casadi_var * m.varlist_all['e0_h_o_j2_i7'].casadi_var))/(1000000.0)) 

    dydx8 = 	 m.varlist_all['e0_p_i1_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i1_j2'].casadi_var) 

    dydx9 = 	 m.varlist_all['e0_p_i2_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i2_j2'].casadi_var) 

    dydx10 = 	 m.varlist_all['e0_p_i3_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i3_j2'].casadi_var) 

    dydx11 = 	 m.varlist_all['e0_p_i4_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i4_j2'].casadi_var) 

    dydx12 = 	 m.varlist_all['e0_p_i5_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i5_j2'].casadi_var) 

    dydx13 = 	 m.varlist_all['e0_p_i6_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i6_j2'].casadi_var) 

    dydx14 = 	 m.varlist_all['e0_p_i7_j2'].casadi_var - (m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_x_i7_j2'].casadi_var) 

    dydx15 = 	 m.varlist_all['e0_U'].casadi_var - (m.varlist_all['e0_u_j2'].casadi_var * m.varlist_all['e0_HU'].casadi_var * 1000.0) 

    dydx16 = 	 m.varlist_all['e0_p_j2'].casadi_var * 100.0 * m.varlist_all['e0_V'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var) 

    dydx17 = 	 m.varlist_all['e0_HU_i1'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i1_j2'].casadi_var) 

    dydx18 = 	 m.varlist_all['e0_HU_i2'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i2_j2'].casadi_var) 

    dydx19 = 	 m.varlist_all['e0_HU_i3'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i3_j2'].casadi_var) 

    dydx20 = 	 m.varlist_all['e0_HU_i4'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i4_j2'].casadi_var) 

    dydx21 = 	 m.varlist_all['e0_HU_i5'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i5_j2'].casadi_var) 

    dydx22 = 	 m.varlist_all['e0_HU_i6'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i6_j2'].casadi_var) 

    dydx23 = 	 m.varlist_all['e0_HU_i7'].casadi_var - (m.varlist_all['e0_HU'].casadi_var * m.varlist_all['e0_x_i7_j2'].casadi_var) 

    dydx24 = 	 m.varlist_all['e0_v_j2'].casadi_var - ((m.varlist_all['e0_V'].casadi_var)/(m.varlist_all['e0_HU'].casadi_var * 1000.0)) 

    dydx25 = 	 m.varlist_all['e0_r_r3'].casadi_var - (m.varlist_all['e0_greek_rho_r3'].casadi_var * m.varlist_all['e0_k_r3'].casadi_var * (m.varlist_all['e0_K_i6'].casadi_var) ** (2.0) * ((m.varlist_all['e0_p_i6_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_p_i7_j2'].casadi_var * m.varlist_all['e0_p_i3_j2'].casadi_var)/(m.varlist_all['e0_K_r3'].casadi_var))/((1.0 + m.varlist_all['e0_K_i6'].casadi_var * m.varlist_all['e0_p_i6_j2'].casadi_var + m.varlist_all['e0_K_i3'].casadi_var * m.varlist_all['e0_p_i3_j2'].casadi_var) ** (2.0)) * m.varlist_all['e0_V'].casadi_var) 

    dydx26 = 	 m.varlist_all['e0_r_r2'].casadi_var - (m.varlist_all['e0_greek_rho_r1'].casadi_var * m.varlist_all['e0_k_par1'].casadi_var * (m.varlist_all['e0_p_i5_j2'].casadi_var * m.varlist_all['e0_p_i2_j2'].casadi_var * (1.0 - (m.varlist_all['e0_p_i3_j2'].casadi_var * m.varlist_all['e0_p_i6_j2'].casadi_var)/((m.varlist_all['e0_p_i2_j2'].casadi_var) ** (3.0) * m.varlist_all['e0_p_i5_j2'].casadi_var * m.varlist_all['e0_K_r2'].casadi_var)))/((1.0 + m.varlist_all['e0_k_par2'].casadi_var * (m.varlist_all['e0_p_i3_j2'].casadi_var)/(m.varlist_all['e0_p_i2_j2'].casadi_var) + m.varlist_all['e0_k_par3'].casadi_var * (m.varlist_all['e0_p_i2_j2'].casadi_var) ** (0.5) + m.varlist_all['e0_k_par4'].casadi_var * (m.varlist_all['e0_p_i3_j2'].casadi_var)) ** (3.0)) * m.varlist_all['e0_V'].casadi_var) 

    dydx27 = 	 m.varlist_all['e0_r_r1'].casadi_var - (m.varlist_all['e0_greek_rho_r1'].casadi_var * m.varlist_all['e0_k_par5'].casadi_var * (m.varlist_all['e0_p_i5_j2'].casadi_var * (1.0 - (m.varlist_all['e0_p_i4_j2'].casadi_var * m.varlist_all['e0_p_i3_j2'].casadi_var)/(m.varlist_all['e0_p_i2_j2'].casadi_var * m.varlist_all['e0_p_i5_j2'].casadi_var * m.varlist_all['e0_K_r1'].casadi_var)))/(1.0 + m.varlist_all['e0_k_par2'].casadi_var * (m.varlist_all['e0_p_i3_j2'].casadi_var)/(m.varlist_all['e0_p_i2_j2'].casadi_var) + m.varlist_all['e0_k_par3'].casadi_var * (m.varlist_all['e0_p_i2_j2'].casadi_var) ** (0.5) + m.varlist_all['e0_k_par4'].casadi_var * (m.varlist_all['e0_p_i3_j2'].casadi_var)) * m.varlist_all['e0_V'].casadi_var) 

    dydx28 = 	 m.varlist_all['e0_u_j2'].casadi_var * 10.0 - (m.varlist_all['e0_h_j2'].casadi_var * 10.0 - m.varlist_all['e0_p_j2'].casadi_var * m.varlist_all['e0_v_j2'].casadi_var) 

    dydx29 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * 1000.0 * m.varlist_all['e0_h_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * 1000.0 * m.varlist_all['e0_h_j2'].casadi_var) + m.varlist_all['e0_Q'].casadi_var) 

    dydx30 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i1_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i1_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i1'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i1'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i1'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx31 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i2_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i2_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i2'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i2'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i2'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx32 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i3_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i3_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i3'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i3'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i3'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx33 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i4_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i4_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i4'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i4'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i4'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx34 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i5_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i5_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i5'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i5'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i5'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx35 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i6_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i6_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i6'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i6'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i6'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx36 = 	 0.0 - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i7_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i7_j2'].casadi_var) + (m.varlist_all['e0_greek_nu_r1_i7'].casadi_var * m.varlist_all['e0_r_r1'].casadi_var + m.varlist_all['e0_greek_nu_r2_i7'].casadi_var * m.varlist_all['e0_r_r2'].casadi_var + m.varlist_all['e0_greek_nu_r3_i7'].casadi_var * m.varlist_all['e0_r_r3'].casadi_var)) 

    dydx37 = 	 m.varlist_all['e0_X_i5'].casadi_var - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i5_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i5_j2'].casadi_var)/(m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i5_j1'].casadi_var)) 

    dydx38 = 	 m.varlist_all['e0_X_i2'].casadi_var - ((m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i2_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i2_j2'].casadi_var)/(m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i2_j1'].casadi_var)) 

    dydx39 = 	 m.varlist_all['e0_Y_i7'].casadi_var - ( - 1.0 * (m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i7_j1'].casadi_var + m.varlist_all['e0_F_j2'].casadi_var * m.varlist_all['e0_x_i7_j2'].casadi_var)/(m.varlist_all['e0_F_j1'].casadi_var * m.varlist_all['e0_x_i5_j1'].casadi_var)) 

    dydx40 = 	 m.varlist_all['e0_h_f_i1'].casadi_var + m.varlist_all['e0_A_par1_i1'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i1'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i1'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i1'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i1'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i1'].casadi_var) 

    dydx41 = 	 m.varlist_all['e0_h_f_i2'].casadi_var + m.varlist_all['e0_A_par1_i2'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i2'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i2'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i2'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i2'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i2'].casadi_var) 

    dydx42 = 	 m.varlist_all['e0_h_f_i3'].casadi_var + m.varlist_all['e0_A_par1_i3'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i3'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i3'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i3'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i3'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i3'].casadi_var) 

    dydx43 = 	 m.varlist_all['e0_h_f_i4'].casadi_var + m.varlist_all['e0_A_par1_i4'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i4'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i4'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i4'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i4'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i4'].casadi_var) 

    dydx44 = 	 m.varlist_all['e0_h_f_i5'].casadi_var + m.varlist_all['e0_A_par1_i5'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i5'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i5'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i5'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i5'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i5'].casadi_var) 

    dydx45 = 	 m.varlist_all['e0_h_f_i6'].casadi_var + m.varlist_all['e0_A_par1_i6'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i6'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i6'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i6'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i6'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i6'].casadi_var) 

    dydx46 = 	 m.varlist_all['e0_h_f_i7'].casadi_var + m.varlist_all['e0_A_par1_i7'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i7'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i7'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i7'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i7'].casadi_var * ((m.varlist_all['e0_T_j1'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j1_i7'].casadi_var) 

    dydx47 = 	 m.varlist_all['e0_h_f_i1'].casadi_var + m.varlist_all['e0_A_par1_i1'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i1'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i1'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i1'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i1'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i1'].casadi_var) 

    dydx48 = 	 m.varlist_all['e0_h_f_i2'].casadi_var + m.varlist_all['e0_A_par1_i2'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i2'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i2'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i2'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i2'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i2'].casadi_var) 

    dydx49 = 	 m.varlist_all['e0_h_f_i3'].casadi_var + m.varlist_all['e0_A_par1_i3'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i3'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i3'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i3'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i3'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i3'].casadi_var) 

    dydx50 = 	 m.varlist_all['e0_h_f_i4'].casadi_var + m.varlist_all['e0_A_par1_i4'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i4'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i4'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i4'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i4'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i4'].casadi_var) 

    dydx51 = 	 m.varlist_all['e0_h_f_i5'].casadi_var + m.varlist_all['e0_A_par1_i5'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i5'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i5'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i5'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i5'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i5'].casadi_var) 

    dydx52 = 	 m.varlist_all['e0_h_f_i6'].casadi_var + m.varlist_all['e0_A_par1_i6'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i6'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i6'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i6'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i6'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i6'].casadi_var) 

    dydx53 = 	 m.varlist_all['e0_h_f_i7'].casadi_var + m.varlist_all['e0_A_par1_i7'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (5.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (5.0))/(5.0) + m.varlist_all['e0_A_par2_i7'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (4.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (4.0))/(4.0) + m.varlist_all['e0_A_par3_i7'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (3.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (3.0))/(3.0) + m.varlist_all['e0_A_par4_i7'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (2.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (2.0))/(2.0) + m.varlist_all['e0_A_par5_i7'].casadi_var * ((m.varlist_all['e0_T_j2'].casadi_var) ** (1.0) - (m.varlist_all['e0_T_f'].casadi_var) ** (1.0))/(1.0) - (m.varlist_all['e0_h_o_j2_i7'].casadi_var) 

    dydx54 = 	 m.varlist_all['e0_A_par1_r1'].casadi_var * ca.exp((m.varlist_all['e0_E_par1_r1'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_k_par1'].casadi_var) 

    dydx55 = 	 m.varlist_all['e0_A_par2_r1'].casadi_var * ca.exp((m.varlist_all['e0_E_par2_r1'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_k_par2'].casadi_var) 

    dydx56 = 	 m.varlist_all['e0_A_par3_r1'].casadi_var * ca.exp((m.varlist_all['e0_E_par3_r1'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_k_par3'].casadi_var) 

    dydx57 = 	 m.varlist_all['e0_A_par4_r1'].casadi_var * ca.exp((m.varlist_all['e0_E_par4_r1'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_k_par4'].casadi_var) 

    dydx58 = 	 m.varlist_all['e0_A_par5_r1'].casadi_var * ca.exp((m.varlist_all['e0_E_par5_r1'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_k_par5'].casadi_var) 

    dydx59 = 	 m.varlist_all['e0_A_r3_i3'].casadi_var * ca.exp((m.varlist_all['e0_E_r3_i3'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_K_i3'].casadi_var) 

    dydx60 = 	 m.varlist_all['e0_A_r3_i6'].casadi_var * ca.exp((m.varlist_all['e0_E_r3_i6'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_K_i6'].casadi_var) 

    dydx61 = 	 m.varlist_all['e0_A_r3'].casadi_var * ca.exp((m.varlist_all['e0_E_r3'].casadi_var)/(m.varlist_all['e0_R'].casadi_var * m.varlist_all['e0_T_j2'].casadi_var)) - (m.varlist_all['e0_k_r3'].casadi_var) 


    # fmt: on"

    m.add_equations_algebraic([dydx1 ,dydx2 ,dydx3 ,dydx4 ,dydx5 ,dydx6 ,dydx7 ,dydx8 ,dydx9 ,dydx10 ,dydx11 ,dydx12 ,dydx13 ,dydx14 ,dydx15 ,dydx16 ,dydx17 ,dydx18 ,dydx19 ,dydx20 ,dydx21 ,dydx22 ,dydx23 ,dydx24 ,dydx25 ,dydx26 ,dydx27 ,dydx28 ,dydx29 ,dydx30 ,dydx31 ,dydx32 ,dydx33 ,dydx34 ,dydx35 ,dydx36 ,dydx37 ,dydx38 ,dydx39 ,dydx40 ,dydx41 ,dydx42 ,dydx43 ,dydx44 ,dydx45 ,dydx46 ,dydx47 ,dydx48 ,dydx49 ,dydx50 ,dydx51 ,dydx52 ,dydx53 ,dydx54 ,dydx55 ,dydx56 ,dydx57 ,dydx58 ,dydx59 ,dydx60 ,dydx61 ,])

    return variable_list, m

if __name__ == "__main__":

    variable_list, m = initialize_problem()

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    var_list_fixed["e0_x_i2_j2"].value.value = 0
    # Create simulation Object
    sim_fixed = SimulatorNLE(m, var_list_fixed)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate_sym()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data()
    print(res["e0_x_i2_j2"].value.value)

    """Erstelle Startwerte"""
    variable_list_optimizer = copy.deepcopy(variable_list)
    for key, var in res.items():
        var.guess = var.value.value
        variable_list_optimizer[key] = var

    """Obere und untere Grenze"""
    for var in variable_list_optimizer.values():
        var.fixed = True
        if isinstance(var, par_est.VariableParameter):
            var.lower_bound = var.value - var.value*0.05
            var.upper_bound = var.value + var.value*0.05
            var.guess = var.lower_bound

    variable_list_optimizer["e0_A_r3_i3"].fixed = False
    variable_list_optimizer["e0_E_par2_r1"].fixed = False
    variable_list_optimizer["e0_A_par5_r1"].fixed = False
    variable_list_optimizer["e0_A_par4_r1"].fixed = False
    variable_list_optimizer["e0_E_par4_r1"].fixed = False

    pe = ParameterEstimationNLE(m, [variable_list_optimizer, variable_list_optimizer])
    res = pe.optimize(False)
    # , 0.5498))

    #res.plot_states()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="	")

