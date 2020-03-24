import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est


def initialize_problem():

    variable_list = par_est.VariableList()

    # fmt: off 
    variable_list.add_variable(par_est.State_variable("e0_c_i1", 2.2564697))
    variable_list.add_variable(par_est.State_variable("e0_c_i2", 0.08573363))
    variable_list.add_variable(par_est.State_variable("e0_c_i3", 1.0E-7))
    variable_list.add_variable(par_est.State_variable("e0_c_i4", 0.00824684))
    variable_list.add_variable(par_est.State_variable("e0_c_i5", 1.0E-7))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i1"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i2"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i3"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i4"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i5"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_L"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_alpha"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_gamma"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_X"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_x_i6"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_x_i7"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_c_i6"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_c_i7"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_psi_cat"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_DeltaG_r3"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_K_eq_r3"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_K_eq_r1"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r1"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r2"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r3"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r4"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r5"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r6"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i1"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i2"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i3"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i4"))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i5"))
    variable_list.add_variable(par_est.Parameter_variable("e0_V_Reactor", 0.872033))
    variable_list.add_variable(par_est.Parameter_variable("e0_n_Surfactant", 0.117))
    variable_list.add_variable(par_est.Parameter_variable("e0_n_Water", 19.1561))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_i1", 168.32))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_i2", 168.32))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_i3", 198.34))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_i4", 170.34))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_i5", 198.34))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_Water", 18.0153))
    variable_list.add_variable(par_est.Parameter_variable("e0_M_Surfactant", 513.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_T", 368.15))
    variable_list.add_variable(par_est.Parameter_variable("e0_p_Reactor", 15.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol1", -6.4909E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol2", 1.1885E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol3", 0.0010631))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol4", -0.027378))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol5", 1.7599E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol6", 0.17476))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol7", 9.2954E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol8", 2.8881E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol9", 2.9467E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol10", 3.7274E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol11", -4.1033E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol12", -9.9645E-6))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol13", -3.8368E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol14", -6.9782E-6))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol15", -8.2558E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol1", -1.7718E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol2", 1.7692E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol3", 0.0016934))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol4", -0.047302))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol5", 4.3746E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol6", 0.28638))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol7", 0.001592))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol8", -1.7107E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol9", 6.5328E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol10", 5.3043E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol11", -7.299E-6))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol12", -1.4868E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol13", -3.0261E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol14", -1.2455E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol15", -1.1598E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_cat_e1", 45087.07))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_cat_e2", 189.31375))
    variable_list.add_variable(par_est.Parameter_variable("e0_c_cat", 0.25682598))
    variable_list.add_variable(par_est.Parameter_variable("e0_R", 8.314))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_DeltaG_r1", 38165.484))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r1", 40749.277))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r1_e1", 0.72770315))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r1_e2", 4.05E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_LM", 2.7251527))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_Surfactant", 1.0315819))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_trig_r1", 14.282191))
    variable_list.add_variable(par_est.Parameter_variable("e0_T_ref", 363.15))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_LM_r1", 66.92345))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r1", 4.242135))
    variable_list.add_variable(par_est.Parameter_variable("e0_n_Cat", 8.58E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_n_Lig", 0.0043))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r2", 6285.8706))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r2", 0.005641299))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r3", 104496.37))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r3_e1", 0.47820565))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r3_e2", 13262.677))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r3_e3", 1028.9795))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r3", 17428.53))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r4", 107045.41))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_trig_Hyfo", 11.312137))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r4", 15349.087))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_LM_Hyfo", 1.0487578))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r5", 57858.113))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r5_e1", 0.023340752))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r5_e2", 895.06036))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r5_e3", 44226.242))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r5", 9.94E7))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r6", 32422.021))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r6", 0.010987442))

    m = par_est.Model(variable_list)

    dydx1 =  m._all_variables["e0_r_i1"].casadi_var  *  60.0 
    dydx2 =  m._all_variables["e0_r_i2"].casadi_var  *  60.0 
    dydx3 =  m._all_variables["e0_r_i3"].casadi_var  *  60.0 
    dydx4 =  m._all_variables["e0_r_i4"].casadi_var  *  60.0 
    dydx5 =  m._all_variables["e0_r_i5"].casadi_var  *  60.0 
    dydx6 =  m._all_variables["e0_n_i1"].casadi_var  -  (  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  ) 
    dydx7 =  m._all_variables["e0_n_i2"].casadi_var  -  (  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  ) 
    dydx8 =  m._all_variables["e0_n_i3"].casadi_var  -  (  m._all_variables["e0_c_i3"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  ) 
    dydx9 =  m._all_variables["e0_n_i4"].casadi_var  -  (  m._all_variables["e0_c_i4"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  ) 
    dydx10 =  m._all_variables["e0_n_i5"].casadi_var  -  (  m._all_variables["e0_c_i5"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  ) 
    dydx11 =  m._all_variables["e0_n_L"].casadi_var  -  (  (  m._all_variables["e0_n_i1"].casadi_var  +  m._all_variables["e0_n_i2"].casadi_var  +  m._all_variables["e0_n_i3"].casadi_var  +  m._all_variables["e0_n_i4"].casadi_var  +  m._all_variables["e0_n_i5"].casadi_var  )  +  m._all_variables["e0_n_Water"].casadi_var  +  m._all_variables["e0_n_Surfactant"].casadi_var  ) 
    dydx12 =  m._all_variables["e0_greek_alpha"].casadi_var  -  (  (  (  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i1"].casadi_var  +  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i2"].casadi_var  +  m._all_variables["e0_c_i3"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i3"].casadi_var  +  m._all_variables["e0_c_i4"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i4"].casadi_var  +  m._all_variables["e0_c_i5"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i5"].casadi_var  )  )/(  (  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i1"].casadi_var  +  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i2"].casadi_var  +  m._all_variables["e0_c_i3"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i3"].casadi_var  +  m._all_variables["e0_c_i4"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i4"].casadi_var  +  m._all_variables["e0_c_i5"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i5"].casadi_var  )  +  m._all_variables["e0_n_Water"].casadi_var  *  m._all_variables["e0_M_Water"].casadi_var  )  ) 
    dydx13 =  m._all_variables["e0_greek_gamma"].casadi_var  -  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  *  m._all_variables["e0_M_Surfactant"].casadi_var  )/(  (  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i1"].casadi_var  +  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i2"].casadi_var  +  m._all_variables["e0_c_i3"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i3"].casadi_var  +  m._all_variables["e0_c_i4"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i4"].casadi_var  +  m._all_variables["e0_c_i5"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i5"].casadi_var  )  +  m._all_variables["e0_n_Water"].casadi_var  *  m._all_variables["e0_M_Water"].casadi_var  +  m._all_variables["e0_n_Surfactant"].casadi_var  *  m._all_variables["e0_M_Surfactant"].casadi_var  )  ) 
    dydx14 =  m._all_variables["e0_X"].casadi_var  -  (  (  (  m._all_variables["e0_c_i3"].casadi_var  *  m._all_variables["e0_M_i3"].casadi_var  +  m._all_variables["e0_c_i5"].casadi_var  *  m._all_variables["e0_M_i5"].casadi_var  )  *  m._all_variables["e0_V_Reactor"].casadi_var  )/(  (  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i1"].casadi_var  +  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i2"].casadi_var  +  m._all_variables["e0_c_i3"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i3"].casadi_var  +  m._all_variables["e0_c_i4"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i4"].casadi_var  +  m._all_variables["e0_c_i5"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  *  m._all_variables["e0_M_i5"].casadi_var  )  )  ) 
    dydx15 =  m._all_variables["e0_x_i6"].casadi_var  -  (  (  m._all_variables["e0_p_Reactor"].casadi_var  *  m._all_variables["e0_P_i6_Sol1"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_P_i6_Sol2"].casadi_var  +  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_P_i6_Sol3"].casadi_var  +  m._all_variables["e0_greek_gamma"].casadi_var  *  m._all_variables["e0_P_i6_Sol4"].casadi_var  +  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i6_Sol5"].casadi_var  +  (  (  m._all_variables["e0_greek_gamma"].casadi_var  )  )**(  2.0  )  *  m._all_variables["e0_P_i6_Sol6"].casadi_var  +  (  (  m._all_variables["e0_X"].casadi_var  )  )**(  2.0  )  *  m._all_variables["e0_P_i6_Sol7"].casadi_var  +  m._all_variables["e0_p_Reactor"].casadi_var  *  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_P_i6_Sol8"].casadi_var  +  m._all_variables["e0_p_Reactor"].casadi_var  *  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_P_i6_Sol9"].casadi_var  +  m._all_variables["e0_p_Reactor"].casadi_var  *  m._all_variables["e0_greek_gamma"].casadi_var  *  m._all_variables["e0_P_i6_Sol10"].casadi_var  +  m._all_variables["e0_p_Reactor"].casadi_var  *  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i6_Sol11"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_P_i6_Sol12"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_greek_gamma"].casadi_var  *  m._all_variables["e0_P_i6_Sol13"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i6_Sol14"].casadi_var  +  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i6_Sol15"].casadi_var  )  -  m._all_variables["e0_x_i7"].casadi_var  ) 
    dydx16 =  m._all_variables["e0_x_i7"].casadi_var  -  (  (  m._all_variables["e0_p_Reactor"].casadi_var  )/(  2.0  )  *  m._all_variables["e0_P_i7_Sol1"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_P_i7_Sol2"].casadi_var  +  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_P_i7_Sol3"].casadi_var  +  m._all_variables["e0_greek_gamma"].casadi_var  *  m._all_variables["e0_P_i7_Sol4"].casadi_var  +  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i7_Sol5"].casadi_var  +  (  (  m._all_variables["e0_greek_gamma"].casadi_var  )  )**(  2.0  )  *  m._all_variables["e0_P_i7_Sol6"].casadi_var  +  (  (  m._all_variables["e0_X"].casadi_var  )  )**(  2.0  )  *  m._all_variables["e0_P_i7_Sol7"].casadi_var  +  (  m._all_variables["e0_p_Reactor"].casadi_var  )/(  2.0  )  *  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_P_i7_Sol8"].casadi_var  +  (  m._all_variables["e0_p_Reactor"].casadi_var  )/(  2.0  )  *  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_P_i7_Sol9"].casadi_var  +  (  m._all_variables["e0_p_Reactor"].casadi_var  )/(  2.0  )  *  m._all_variables["e0_greek_gamma"].casadi_var  *  m._all_variables["e0_P_i7_Sol10"].casadi_var  +  (  m._all_variables["e0_p_Reactor"].casadi_var  )/(  2.0  )  *  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i7_Sol11"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_P_i7_Sol12"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_greek_gamma"].casadi_var  *  m._all_variables["e0_P_i7_Sol13"].casadi_var  +  (  m._all_variables["e0_T"].casadi_var  -  273.15  )  *  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i7_Sol14"].casadi_var  +  m._all_variables["e0_greek_alpha"].casadi_var  *  m._all_variables["e0_X"].casadi_var  *  m._all_variables["e0_P_i7_Sol15"].casadi_var  ) 
    dydx17 =  m._all_variables["e0_c_i6"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  -  (  (  m._all_variables["e0_n_L"].casadi_var  *  m._all_variables["e0_x_i6"].casadi_var  )/(  1.0  -  m._all_variables["e0_x_i6"].casadi_var  )  ) 
    dydx18 =  m._all_variables["e0_c_i7"].casadi_var  *  m._all_variables["e0_V_Reactor"].casadi_var  -  (  (  m._all_variables["e0_n_L"].casadi_var  *  m._all_variables["e0_x_i7"].casadi_var  )/(  1.0  -  m._all_variables["e0_x_i7"].casadi_var  )  ) 
    dydx19 =  m._all_variables["e0_greek_psi_cat"].casadi_var  *  (  1.0  +  m._all_variables["e0_K_cat_e1"].casadi_var  *  m._all_variables["e0_c_i7"].casadi_var  +  m._all_variables["e0_K_cat_e2"].casadi_var  *  (  m._all_variables["e0_c_i7"].casadi_var  )/(  m._all_variables["e0_c_i6"].casadi_var  )  )  -  (  m._all_variables["e0_c_cat"].casadi_var  ) 
    dydx20 =  m._all_variables["e0_greek_DeltaG_r3"].casadi_var  -  (  (  -  126.28  +  0.13  *  m._all_variables["e0_T"].casadi_var  +  6.8  *  (  (  10.0  )  )**(  -  6.0  )  *  (  (  m._all_variables["e0_T"].casadi_var  )  )**(  2.0  )  )  *  (  (  10.0  )  )**(  3.0  )  ) 
    dydx21 =  m._all_variables["e0_K_eq_r3"].casadi_var  -  (  ca.exp(  -  (  m._all_variables["e0_greek_DeltaG_r3"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  *  m._all_variables["e0_T"].casadi_var  )  )  ) 
    dydx22 =  m._all_variables["e0_K_eq_r1"].casadi_var  -  (  ca.exp(  (  m._all_variables["e0_greek_DeltaG_r1"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  *  m._all_variables["e0_T"].casadi_var  )  )  ) 
    dydx23 =  m._all_variables["e0_r_r1"].casadi_var  *  (  1.0  +  m._all_variables["e0_K_r1_e1"].casadi_var  *  m._all_variables["e0_c_i1"].casadi_var  +  m._all_variables["e0_K_r1_e2"].casadi_var  *  m._all_variables["e0_c_i2"].casadi_var  )  -  (  (  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  )/(  m._all_variables["e0_V_Reactor"].casadi_var  )  )  )**(  m._all_variables["e0_P_Surfactant"].casadi_var  )  *  (  1.0  +  (  m._all_variables["e0_k_LM_r1"].casadi_var  )/(  1.0  +  ca.exp(  -  (  m._all_variables["e0_K_LM"].casadi_var  -  (  m._all_variables["e0_n_Lig"].casadi_var  )/(  m._all_variables["e0_n_Cat"].casadi_var  )  )  *  m._all_variables["e0_P_trig_r1"].casadi_var  )  )  )  *  m._all_variables["e0_greek_psi_cat"].casadi_var  *  m._all_variables["e0_k_ref_r1"].casadi_var  *  ca.exp(  -  (  m._all_variables["e0_E_r1"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  )  *  (  (  1.0  )/(  m._all_variables["e0_T"].casadi_var  )  -  (  1.0  )/(  m._all_variables["e0_T_ref"].casadi_var  )  )  )  *  (  m._all_variables["e0_c_i1"].casadi_var  -  (  m._all_variables["e0_c_i2"].casadi_var  )/(  m._all_variables["e0_K_eq_r1"].casadi_var  )  )  ) 
    dydx24 =  m._all_variables["e0_r_r2"].casadi_var  -  (  (  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  )/(  m._all_variables["e0_V_Reactor"].casadi_var  )  )  )**(  m._all_variables["e0_P_Surfactant"].casadi_var  )  *  m._all_variables["e0_greek_psi_cat"].casadi_var  *  m._all_variables["e0_k_ref_r2"].casadi_var  *  ca.exp(  -  (  m._all_variables["e0_E_r2"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  )  *  (  (  1.0  )/(  m._all_variables["e0_T"].casadi_var  )  -  (  1.0  )/(  m._all_variables["e0_T_ref"].casadi_var  )  )  )  *  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  ) 
    dydx25 =  m._all_variables["e0_r_r3"].casadi_var  *  (  1.0  +  m._all_variables["e0_K_r3_e1"].casadi_var  *  m._all_variables["e0_c_i1"].casadi_var  +  m._all_variables["e0_K_r3_e2"].casadi_var  *  m._all_variables["e0_c_i4"].casadi_var  +  m._all_variables["e0_K_r3_e3"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  )  -  (  (  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  )/(  m._all_variables["e0_V_Reactor"].casadi_var  )  )  )**(  m._all_variables["e0_P_Surfactant"].casadi_var  )  *  m._all_variables["e0_greek_psi_cat"].casadi_var  *  m._all_variables["e0_k_ref_r3"].casadi_var  *  ca.exp(  -  (  m._all_variables["e0_E_r3"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  )  *  (  (  1.0  )/(  m._all_variables["e0_T"].casadi_var  )  -  (  1.0  )/(  m._all_variables["e0_T_ref"].casadi_var  )  )  )  *  (  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  -  (  m._all_variables["e0_c_i4"].casadi_var  )/(  m._all_variables["e0_K_eq_r3"].casadi_var  )  )  ) 
    dydx26 =  m._all_variables["e0_r_r4"].casadi_var  -  (  (  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  )/(  m._all_variables["e0_V_Reactor"].casadi_var  )  )  )**(  m._all_variables["e0_P_Surfactant"].casadi_var  )  *  (  1.0  +  (  m._all_variables["e0_k_LM_Hyfo"].casadi_var  )/(  1.0  +  ca.exp(  -  (  m._all_variables["e0_K_LM"].casadi_var  -  (  m._all_variables["e0_n_Lig"].casadi_var  )/(  m._all_variables["e0_n_Cat"].casadi_var  )  )  *  m._all_variables["e0_P_trig_Hyfo"].casadi_var  )  )  )  *  m._all_variables["e0_greek_psi_cat"].casadi_var  *  m._all_variables["e0_k_ref_r4"].casadi_var  *  ca.exp(  -  (  m._all_variables["e0_E_r4"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  )  *  (  (  1.0  )/(  m._all_variables["e0_T"].casadi_var  )  -  (  1.0  )/(  m._all_variables["e0_T_ref"].casadi_var  )  )  )  *  m._all_variables["e0_c_i2"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  *  m._all_variables["e0_c_i7"].casadi_var  ) 
    dydx27 =  m._all_variables["e0_r_r5"].casadi_var  *  (  1.0  +  m._all_variables["e0_K_r5_e1"].casadi_var  *  m._all_variables["e0_c_i1"].casadi_var  +  m._all_variables["e0_K_r5_e2"].casadi_var  *  m._all_variables["e0_c_i5"].casadi_var  +  m._all_variables["e0_K_r5_e3"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  )  -  (  (  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  )/(  m._all_variables["e0_V_Reactor"].casadi_var  )  )  )**(  m._all_variables["e0_P_Surfactant"].casadi_var  )  *  (  1.0  +  (  m._all_variables["e0_k_LM_Hyfo"].casadi_var  )/(  1.0  +  ca.exp(  -  (  m._all_variables["e0_K_LM"].casadi_var  -  (  m._all_variables["e0_n_Lig"].casadi_var  )/(  m._all_variables["e0_n_Cat"].casadi_var  )  )  *  m._all_variables["e0_P_trig_Hyfo"].casadi_var  )  )  )  *  m._all_variables["e0_greek_psi_cat"].casadi_var  *  m._all_variables["e0_k_ref_r5"].casadi_var  *  ca.exp(  -  (  m._all_variables["e0_E_r5"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  )  *  (  (  1.0  )/(  m._all_variables["e0_T"].casadi_var  )  -  (  1.0  )/(  m._all_variables["e0_T_ref"].casadi_var  )  )  )  *  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  *  m._all_variables["e0_c_i7"].casadi_var  ) 
    dydx28 =  m._all_variables["e0_r_r6"].casadi_var  -  (  (  (  (  m._all_variables["e0_n_Surfactant"].casadi_var  )/(  m._all_variables["e0_V_Reactor"].casadi_var  )  )  )**(  m._all_variables["e0_P_Surfactant"].casadi_var  )  *  (  1.0  +  (  m._all_variables["e0_k_LM_Hyfo"].casadi_var  )/(  1.0  +  ca.exp(  -  (  m._all_variables["e0_K_LM"].casadi_var  -  (  m._all_variables["e0_n_Lig"].casadi_var  )/(  m._all_variables["e0_n_Cat"].casadi_var  )  )  *  m._all_variables["e0_P_trig_Hyfo"].casadi_var  )  )  )  *  m._all_variables["e0_greek_psi_cat"].casadi_var  *  m._all_variables["e0_k_ref_r6"].casadi_var  *  ca.exp(  -  (  m._all_variables["e0_E_r6"].casadi_var  )/(  m._all_variables["e0_R"].casadi_var  )  *  (  (  1.0  )/(  m._all_variables["e0_T"].casadi_var  )  -  (  1.0  )/(  m._all_variables["e0_T_ref"].casadi_var  )  )  )  *  m._all_variables["e0_c_i1"].casadi_var  *  m._all_variables["e0_c_i6"].casadi_var  *  m._all_variables["e0_c_i7"].casadi_var  ) 
    dydx29 =  m._all_variables["e0_r_i1"].casadi_var  -  (  -  m._all_variables["e0_r_r1"].casadi_var  -  m._all_variables["e0_r_r3"].casadi_var  -  m._all_variables["e0_r_r5"].casadi_var  -  m._all_variables["e0_r_r6"].casadi_var  ) 
    dydx30 =  m._all_variables["e0_r_i2"].casadi_var  -  (  m._all_variables["e0_r_r1"].casadi_var  -  m._all_variables["e0_r_r2"].casadi_var  -  m._all_variables["e0_r_r4"].casadi_var  ) 
    dydx31 =  m._all_variables["e0_r_i3"].casadi_var  -  (  m._all_variables["e0_r_r4"].casadi_var  +  m._all_variables["e0_r_r6"].casadi_var  ) 
    dydx32 =  m._all_variables["e0_r_i4"].casadi_var  -  (  m._all_variables["e0_r_r2"].casadi_var  +  m._all_variables["e0_r_r3"].casadi_var  ) 
    dydx33 =  m._all_variables["e0_r_i5"].casadi_var  -  (  m._all_variables["e0_r_r5"].casadi_var  ) 

    # fmt: on"

    m.add_differential_equations([dydx1 ,dydx2 ,dydx3 ,dydx4 ,dydx5 ,])
    m.add_algebraic_equations([dydx6 ,dydx7 ,dydx8 ,dydx9 ,dydx10 ,dydx11 ,dydx12 ,dydx13 ,dydx14 ,dydx15 ,dydx16 ,dydx17 ,dydx18 ,dydx19 ,dydx20 ,dydx21 ,dydx22 ,dydx23 ,dydx24 ,dydx25 ,dydx26 ,dydx27 ,dydx28 ,dydx29 ,dydx30 ,dydx31 ,dydx32 ,dydx33 ,])

    return variable_list, m

if __name__ == "__main__":

    variable_list, m = initialize_problem()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 40)
    time_grid = np.insert(time_grid, 0, 0)

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data()
    res.plot_states()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="	")