import copy

import casadi as ca
import numpy as np

import par_est


def initialize_problem():

    variable_list = par_est.VariableList()

    variable_list.add_variable(par_est.State_variable("e0_c_i1", 2.2564697))
    variable_list.add_variable(par_est.State_variable("e0_c_i2", 0.08573363))
    variable_list.add_variable(par_est.State_variable("e0_c_i3", 1.0e-7))
    variable_list.add_variable(par_est.State_variable("e0_c_i4", 0.00824684))
    variable_list.add_variable(par_est.State_variable("e0_c_i5", 1.0e-7))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i1", 1.967716))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i2", 0.07476256))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i3", 8.72e-8))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i4", 0.007191517))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_i5", 8.72e-8))
    variable_list.add_variable(par_est.Algebraic_variable("e0_n_L", 21.322771))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_alpha", 0.49993634))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_gamma", 0.08001318))
    variable_list.add_variable(par_est.Algebraic_variable("e0_X", 1.0e-7))
    variable_list.add_variable(par_est.Algebraic_variable("e0_x_i6", 9.62e-4))
    variable_list.add_variable(par_est.Algebraic_variable("e0_x_i7", 9.57e-4))
    variable_list.add_variable(par_est.Algebraic_variable("e0_c_i6", 0.023546033))
    variable_list.add_variable(par_est.Algebraic_variable("e0_c_i7", 0.023419946))
    variable_list.add_variable(par_est.Algebraic_variable("e0_greek_psi_cat", 2.06e-4))
    variable_list.add_variable(
        par_est.Algebraic_variable("e0_greek_DeltaG_r3", -77498.87)
    )
    variable_list.add_variable(par_est.Algebraic_variable("e0_K_eq_r3", 9.91e10))
    variable_list.add_variable(par_est.Algebraic_variable("e0_K_eq_r1", 260178.73))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r1", 1.13e-4))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r2", 3.04e-10))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r3", 1.08e-5))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r4", 3.05e-5))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r5", 0.003999636))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_r6", 4.11e-10))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i1", -0.004123436))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i2", 8.25e-5))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i3", 3.05e-5))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i4", 1.08e-5))
    variable_list.add_variable(par_est.Algebraic_variable("e0_r_i5", 0.003999636))
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
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol1", -6.4909e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol2", 1.1885e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol3", 0.0010631))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol4", -0.027378))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol5", 1.7599e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol6", 0.17476))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol7", 9.2954e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol8", 2.8881e-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol9", 2.9467e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol10", 3.7274e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol11", -4.1033e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol12", -9.9645e-6))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol13", -3.8368e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol14", -6.9782e-6))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i6_Sol15", -8.2558e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol1", -1.7718e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol2", 1.7692e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol3", 0.0016934))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol4", -0.047302))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol5", 4.3746e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol6", 0.28638))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol7", 0.001592))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol8", -1.7107e-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol9", 6.5328e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol10", 5.3043e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol11", -7.299e-6))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol12", -1.4868e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol13", -3.0261e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol14", -1.2455e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_i7_Sol15", -1.1598e-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_cat_e1", 45087.07))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_cat_e2", 189.31375))
    variable_list.add_variable(par_est.Parameter_variable("e0_c_cat", 0.25682598))
    variable_list.add_variable(par_est.Parameter_variable("e0_R", 8.314))
    variable_list.add_variable(
        par_est.Parameter_variable("e0_greek_DeltaG_r1", 38165.484)
    )
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r1", 40749.277))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r1_e1", 0.72770315))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_r1_e2", 4.05e-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_K_LM", 2.7251527))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_Surfactant", 1.0315819))
    variable_list.add_variable(par_est.Parameter_variable("e0_P_trig_r1", 14.282191))
    variable_list.add_variable(par_est.Parameter_variable("e0_T_ref", 363.15))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_LM_r1", 66.92345))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r1", 4.242135))
    variable_list.add_variable(par_est.Parameter_variable("e0_n_Cat", 8.58e-4))
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
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r5", 9.94e7))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r6", 32422.021))
    variable_list.add_variable(par_est.Parameter_variable("e0_k_ref_r6", 0.010987442))

    m = par_est.Model(variable_list)

    dydx1 = m._all_variables["e0_r_i1"].casadi_var * 60.0
    dydx2 = m._all_variables["e0_r_i2"].casadi_var * 60.0
    dydx3 = m._all_variables["e0_r_i3"].casadi_var * 60.0
    dydx4 = m._all_variables["e0_r_i4"].casadi_var * 60.0
    dydx5 = m._all_variables["e0_r_i5"].casadi_var * 60.0
    dydx6 = m._all_variables["e0_n_i1"].casadi_var - (
        m._all_variables["e0_c_i1"].casadi_var
        * m._all_variables["e0_V_Reactor"].casadi_var
    )
    dydx7 = m._all_variables["e0_n_i2"].casadi_var - (
        m._all_variables["e0_c_i2"].casadi_var
        * m._all_variables["e0_V_Reactor"].casadi_var
    )
    dydx8 = m._all_variables["e0_n_i3"].casadi_var - (
        m._all_variables["e0_c_i3"].casadi_var
        * m._all_variables["e0_V_Reactor"].casadi_var
    )
    dydx9 = m._all_variables["e0_n_i4"].casadi_var - (
        m._all_variables["e0_c_i4"].casadi_var
        * m._all_variables["e0_V_Reactor"].casadi_var
    )
    dydx10 = m._all_variables["e0_n_i5"].casadi_var - (
        m._all_variables["e0_c_i5"].casadi_var
        * m._all_variables["e0_V_Reactor"].casadi_var
    )
    dydx11 = m._all_variables["e0_n_L"].casadi_var - (
        (
            m._all_variables["e0_n_i1"].casadi_var
            + m._all_variables["e0_n_i2"].casadi_var
            + m._all_variables["e0_n_i3"].casadi_var
            + m._all_variables["e0_n_i4"].casadi_var
            + m._all_variables["e0_n_i5"].casadi_var
        )
        + m._all_variables["e0_n_Water"].casadi_var
        + m._all_variables["e0_n_Surfactant"].casadi_var
    )
    dydx12 = m._all_variables["e0_greek_alpha"].casadi_var - (
        (
            (
                m._all_variables["e0_c_i1"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i1"].casadi_var
                + m._all_variables["e0_c_i2"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i2"].casadi_var
                + m._all_variables["e0_c_i3"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i3"].casadi_var
                + m._all_variables["e0_c_i4"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i4"].casadi_var
                + m._all_variables["e0_c_i5"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i5"].casadi_var
            )
        )
        / (
            (
                m._all_variables["e0_c_i1"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i1"].casadi_var
                + m._all_variables["e0_c_i2"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i2"].casadi_var
                + m._all_variables["e0_c_i3"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i3"].casadi_var
                + m._all_variables["e0_c_i4"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i4"].casadi_var
                + m._all_variables["e0_c_i5"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i5"].casadi_var
            )
            + m._all_variables["e0_n_Water"].casadi_var
            * m._all_variables["e0_M_Water"].casadi_var
        )
    )
    dydx13 = m._all_variables["e0_greek_gamma"].casadi_var - (
        (
            m._all_variables["e0_n_Surfactant"].casadi_var
            * m._all_variables["e0_M_Surfactant"].casadi_var
        )
        / (
            (
                m._all_variables["e0_c_i1"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i1"].casadi_var
                + m._all_variables["e0_c_i2"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i2"].casadi_var
                + m._all_variables["e0_c_i3"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i3"].casadi_var
                + m._all_variables["e0_c_i4"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i4"].casadi_var
                + m._all_variables["e0_c_i5"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i5"].casadi_var
            )
            + m._all_variables["e0_n_Water"].casadi_var
            * m._all_variables["e0_M_Water"].casadi_var
            + m._all_variables["e0_n_Surfactant"].casadi_var
            * m._all_variables["e0_M_Surfactant"].casadi_var
        )
    )
    dydx14 = m._all_variables["e0_X"].casadi_var - (
        (
            (
                m._all_variables["e0_c_i3"].casadi_var
                * m._all_variables["e0_M_i3"].casadi_var
                + m._all_variables["e0_c_i5"].casadi_var
                * m._all_variables["e0_M_i5"].casadi_var
            )
            * m._all_variables["e0_V_Reactor"].casadi_var
        )
        / (
            (
                m._all_variables["e0_c_i1"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i1"].casadi_var
                + m._all_variables["e0_c_i2"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i2"].casadi_var
                + m._all_variables["e0_c_i3"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i3"].casadi_var
                + m._all_variables["e0_c_i4"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i4"].casadi_var
                + m._all_variables["e0_c_i5"].casadi_var
                * m._all_variables["e0_V_Reactor"].casadi_var
                * m._all_variables["e0_M_i5"].casadi_var
            )
        )
    )
    dydx15 = m._all_variables["e0_x_i6"].casadi_var - (
        (
            m._all_variables["e0_p_Reactor"].casadi_var
            * m._all_variables["e0_P_i6_Sol1"].casadi_var
            + (m._all_variables["e0_T"].casadi_var - 273.15)
            * m._all_variables["e0_P_i6_Sol2"].casadi_var
            + m._all_variables["e0_greek_alpha"].casadi_var
            * m._all_variables["e0_P_i6_Sol3"].casadi_var
            + m._all_variables["e0_greek_gamma"].casadi_var
            * m._all_variables["e0_P_i6_Sol4"].casadi_var
            + m._all_variables["e0_X"].casadi_var
            * m._all_variables["e0_P_i6_Sol5"].casadi_var
            + ((m._all_variables["e0_greek_gamma"].casadi_var)) ** (2.0)
            * m._all_variables["e0_P_i6_Sol6"].casadi_var
            + ((m._all_variables["e0_X"].casadi_var)) ** (2.0)
            * m._all_variables["e0_P_i6_Sol7"].casadi_var
            + m._all_variables["e0_p_Reactor"].casadi_var
            * (m._all_variables["e0_T"].casadi_var - 273.15)
            * m._all_variables["e0_P_i6_Sol8"].casadi_var
            + m._all_variables["e0_p_Reactor"].casadi_var
            * m._all_variables["e0_greek_alpha"].casadi_var
            * m._all_variables["e0_P_i6_Sol9"].casadi_var
            + m._all_variables["e0_p_Reactor"].casadi_var
            * m._all_variables["e0_greek_gamma"].casadi_var
            * m._all_variables["e0_P_i6_Sol10"].casadi_var
            + m._all_variables["e0_p_Reactor"].casadi_var
            * m._all_variables["e0_X"].casadi_var
            * m._all_variables["e0_P_i6_Sol11"].casadi_var
            + (m._all_variables["e0_T"].casadi_var - 273.15)
            * m._all_variables["e0_greek_alpha"].casadi_var
            * m._all_variables["e0_P_i6_Sol12"].casadi_var
            + (m._all_variables["e0_T"].casadi_var - 273.15)
            * m._all_variables["e0_greek_gamma"].casadi_var
            * m._all_variables["e0_P_i6_Sol13"].casadi_var
            + (m._all_variables["e0_T"].casadi_var - 273.15)
            * m._all_variables["e0_X"].casadi_var
            * m._all_variables["e0_P_i6_Sol14"].casadi_var
            + m._all_variables["e0_greek_alpha"].casadi_var
            * m._all_variables["e0_X"].casadi_var
            * m._all_variables["e0_P_i6_Sol15"].casadi_var
        )
        - m._all_variables["e0_x_i7"].casadi_var
    )
    dydx16 = m._all_variables["e0_x_i7"].casadi_var - (
        (m._all_variables["e0_p_Reactor"].casadi_var)
        / (2.0)
        * m._all_variables["e0_P_i7_Sol1"].casadi_var
        + (m._all_variables["e0_T"].casadi_var - 273.15)
        * m._all_variables["e0_P_i7_Sol2"].casadi_var
        + m._all_variables["e0_greek_alpha"].casadi_var
        * m._all_variables["e0_P_i7_Sol3"].casadi_var
        + m._all_variables["e0_greek_gamma"].casadi_var
        * m._all_variables["e0_P_i7_Sol4"].casadi_var
        + m._all_variables["e0_X"].casadi_var
        * m._all_variables["e0_P_i7_Sol5"].casadi_var
        + ((m._all_variables["e0_greek_gamma"].casadi_var)) ** (2.0)
        * m._all_variables["e0_P_i7_Sol6"].casadi_var
        + ((m._all_variables["e0_X"].casadi_var)) ** (2.0)
        * m._all_variables["e0_P_i7_Sol7"].casadi_var
        + (m._all_variables["e0_p_Reactor"].casadi_var)
        / (2.0)
        * (m._all_variables["e0_T"].casadi_var - 273.15)
        * m._all_variables["e0_P_i7_Sol8"].casadi_var
        + (m._all_variables["e0_p_Reactor"].casadi_var)
        / (2.0)
        * m._all_variables["e0_greek_alpha"].casadi_var
        * m._all_variables["e0_P_i7_Sol9"].casadi_var
        + (m._all_variables["e0_p_Reactor"].casadi_var)
        / (2.0)
        * m._all_variables["e0_greek_gamma"].casadi_var
        * m._all_variables["e0_P_i7_Sol10"].casadi_var
        + (m._all_variables["e0_p_Reactor"].casadi_var)
        / (2.0)
        * m._all_variables["e0_X"].casadi_var
        * m._all_variables["e0_P_i7_Sol11"].casadi_var
        + (m._all_variables["e0_T"].casadi_var - 273.15)
        * m._all_variables["e0_greek_alpha"].casadi_var
        * m._all_variables["e0_P_i7_Sol12"].casadi_var
        + (m._all_variables["e0_T"].casadi_var - 273.15)
        * m._all_variables["e0_greek_gamma"].casadi_var
        * m._all_variables["e0_P_i7_Sol13"].casadi_var
        + (m._all_variables["e0_T"].casadi_var - 273.15)
        * m._all_variables["e0_X"].casadi_var
        * m._all_variables["e0_P_i7_Sol14"].casadi_var
        + m._all_variables["e0_greek_alpha"].casadi_var
        * m._all_variables["e0_X"].casadi_var
        * m._all_variables["e0_P_i7_Sol15"].casadi_var
    )
    dydx17 = m._all_variables["e0_c_i6"].casadi_var * m._all_variables[
        "e0_V_Reactor"
    ].casadi_var - (
        (m._all_variables["e0_n_L"].casadi_var * m._all_variables["e0_x_i6"].casadi_var)
        / (1.0 - m._all_variables["e0_x_i6"].casadi_var)
    )
    dydx18 = m._all_variables["e0_c_i7"].casadi_var * m._all_variables[
        "e0_V_Reactor"
    ].casadi_var - (
        (m._all_variables["e0_n_L"].casadi_var * m._all_variables["e0_x_i7"].casadi_var)
        / (1.0 - m._all_variables["e0_x_i7"].casadi_var)
    )
    dydx19 = m._all_variables["e0_greek_psi_cat"].casadi_var * (
        1.0
        + m._all_variables["e0_K_cat_e1"].casadi_var
        * m._all_variables["e0_c_i7"].casadi_var
        + m._all_variables["e0_K_cat_e2"].casadi_var
        * (m._all_variables["e0_c_i7"].casadi_var)
        / (m._all_variables["e0_c_i6"].casadi_var)
    ) - (m._all_variables["e0_c_cat"].casadi_var)
    dydx20 = m._all_variables["e0_greek_DeltaG_r3"].casadi_var - (
        (
            -126.28
            + 0.13 * m._all_variables["e0_T"].casadi_var
            + 6.8
            * ((10.0)) ** (-6.0)
            * ((m._all_variables["e0_T"].casadi_var)) ** (2.0)
        )
        * ((10.0)) ** (3.0)
    )
    dydx21 = m._all_variables["e0_K_eq_r3"].casadi_var - (
        ca.exp(
            -(m._all_variables["e0_greek_DeltaG_r3"].casadi_var)
            / (
                m._all_variables["e0_R"].casadi_var
                * m._all_variables["e0_T"].casadi_var
            )
        )
    )
    dydx22 = m._all_variables["e0_K_eq_r1"].casadi_var - (
        ca.exp(
            (m._all_variables["e0_greek_DeltaG_r1"].casadi_var)
            / (
                m._all_variables["e0_R"].casadi_var
                * m._all_variables["e0_T"].casadi_var
            )
        )
    )
    dydx23 = m._all_variables["e0_r_r1"].casadi_var * (
        1.0
        + m._all_variables["e0_K_r1_e1"].casadi_var
        * m._all_variables["e0_c_i1"].casadi_var
        + m._all_variables["e0_K_r1_e2"].casadi_var
        * m._all_variables["e0_c_i2"].casadi_var
    ) - (
        (
            (
                (m._all_variables["e0_n_Surfactant"].casadi_var)
                / (m._all_variables["e0_V_Reactor"].casadi_var)
            )
        )
        ** (m._all_variables["e0_P_Surfactant"].casadi_var)
        * (
            1.0
            + (m._all_variables["e0_k_LM_r1"].casadi_var)
            / (
                1.0
                + ca.exp(
                    -(
                        m._all_variables["e0_K_LM"].casadi_var
                        - (m._all_variables["e0_n_Lig"].casadi_var)
                        / (m._all_variables["e0_n_Cat"].casadi_var)
                    )
                    * m._all_variables["e0_P_trig_r1"].casadi_var
                )
            )
        )
        * m._all_variables["e0_greek_psi_cat"].casadi_var
        * m._all_variables["e0_k_ref_r1"].casadi_var
        * ca.exp(
            -(m._all_variables["e0_E_r1"].casadi_var)
            / (m._all_variables["e0_R"].casadi_var)
            * (
                (1.0) / (m._all_variables["e0_T"].casadi_var)
                - (1.0) / (m._all_variables["e0_T_ref"].casadi_var)
            )
        )
        * (
            m._all_variables["e0_c_i1"].casadi_var
            - (m._all_variables["e0_c_i2"].casadi_var)
            / (m._all_variables["e0_K_eq_r1"].casadi_var)
        )
    )
    dydx24 = m._all_variables["e0_r_r2"].casadi_var - (
        (
            (
                (m._all_variables["e0_n_Surfactant"].casadi_var)
                / (m._all_variables["e0_V_Reactor"].casadi_var)
            )
        )
        ** (m._all_variables["e0_P_Surfactant"].casadi_var)
        * m._all_variables["e0_greek_psi_cat"].casadi_var
        * m._all_variables["e0_k_ref_r2"].casadi_var
        * ca.exp(
            -(m._all_variables["e0_E_r2"].casadi_var)
            / (m._all_variables["e0_R"].casadi_var)
            * (
                (1.0) / (m._all_variables["e0_T"].casadi_var)
                - (1.0) / (m._all_variables["e0_T_ref"].casadi_var)
            )
        )
        * m._all_variables["e0_c_i2"].casadi_var
        * m._all_variables["e0_c_i6"].casadi_var
    )
    dydx25 = m._all_variables["e0_r_r3"].casadi_var * (
        1.0
        + m._all_variables["e0_K_r3_e1"].casadi_var
        * m._all_variables["e0_c_i1"].casadi_var
        + m._all_variables["e0_K_r3_e2"].casadi_var
        * m._all_variables["e0_c_i4"].casadi_var
        + m._all_variables["e0_K_r3_e3"].casadi_var
        * m._all_variables["e0_c_i6"].casadi_var
    ) - (
        (
            (
                (m._all_variables["e0_n_Surfactant"].casadi_var)
                / (m._all_variables["e0_V_Reactor"].casadi_var)
            )
        )
        ** (m._all_variables["e0_P_Surfactant"].casadi_var)
        * m._all_variables["e0_greek_psi_cat"].casadi_var
        * m._all_variables["e0_k_ref_r3"].casadi_var
        * ca.exp(
            -(m._all_variables["e0_E_r3"].casadi_var)
            / (m._all_variables["e0_R"].casadi_var)
            * (
                (1.0) / (m._all_variables["e0_T"].casadi_var)
                - (1.0) / (m._all_variables["e0_T_ref"].casadi_var)
            )
        )
        * (
            m._all_variables["e0_c_i2"].casadi_var
            * m._all_variables["e0_c_i6"].casadi_var
            - (m._all_variables["e0_c_i4"].casadi_var)
            / (m._all_variables["e0_K_eq_r3"].casadi_var)
        )
    )
    dydx26 = m._all_variables["e0_r_r4"].casadi_var - (
        (
            (
                (m._all_variables["e0_n_Surfactant"].casadi_var)
                / (m._all_variables["e0_V_Reactor"].casadi_var)
            )
        )
        ** (m._all_variables["e0_P_Surfactant"].casadi_var)
        * (
            1.0
            + (m._all_variables["e0_k_LM_Hyfo"].casadi_var)
            / (
                1.0
                + ca.exp(
                    -(
                        m._all_variables["e0_K_LM"].casadi_var
                        - (m._all_variables["e0_n_Lig"].casadi_var)
                        / (m._all_variables["e0_n_Cat"].casadi_var)
                    )
                    * m._all_variables["e0_P_trig_Hyfo"].casadi_var
                )
            )
        )
        * m._all_variables["e0_greek_psi_cat"].casadi_var
        * m._all_variables["e0_k_ref_r4"].casadi_var
        * ca.exp(
            -(m._all_variables["e0_E_r4"].casadi_var)
            / (m._all_variables["e0_R"].casadi_var)
            * (
                (1.0) / (m._all_variables["e0_T"].casadi_var)
                - (1.0) / (m._all_variables["e0_T_ref"].casadi_var)
            )
        )
        * m._all_variables["e0_c_i2"].casadi_var
        * m._all_variables["e0_c_i6"].casadi_var
        * m._all_variables["e0_c_i7"].casadi_var
    )
    dydx27 = m._all_variables["e0_r_r5"].casadi_var * (
        1.0
        + m._all_variables["e0_K_r5_e1"].casadi_var
        * m._all_variables["e0_c_i1"].casadi_var
        + m._all_variables["e0_K_r5_e2"].casadi_var
        * m._all_variables["e0_c_i5"].casadi_var
        + m._all_variables["e0_K_r5_e3"].casadi_var
        * m._all_variables["e0_c_i6"].casadi_var
    ) - (
        (
            (
                (m._all_variables["e0_n_Surfactant"].casadi_var)
                / (m._all_variables["e0_V_Reactor"].casadi_var)
            )
        )
        ** (m._all_variables["e0_P_Surfactant"].casadi_var)
        * (
            1.0
            + (m._all_variables["e0_k_LM_Hyfo"].casadi_var)
            / (
                1.0
                + ca.exp(
                    -(
                        m._all_variables["e0_K_LM"].casadi_var
                        - (m._all_variables["e0_n_Lig"].casadi_var)
                        / (m._all_variables["e0_n_Cat"].casadi_var)
                    )
                    * m._all_variables["e0_P_trig_Hyfo"].casadi_var
                )
            )
        )
        * m._all_variables["e0_greek_psi_cat"].casadi_var
        * m._all_variables["e0_k_ref_r5"].casadi_var
        * ca.exp(
            -(m._all_variables["e0_E_r5"].casadi_var)
            / (m._all_variables["e0_R"].casadi_var)
            * (
                (1.0) / (m._all_variables["e0_T"].casadi_var)
                - (1.0) / (m._all_variables["e0_T_ref"].casadi_var)
            )
        )
        * m._all_variables["e0_c_i1"].casadi_var
        * m._all_variables["e0_c_i6"].casadi_var
        * m._all_variables["e0_c_i7"].casadi_var
    )
    dydx28 = m._all_variables["e0_r_r6"].casadi_var - (
        (
            (
                (m._all_variables["e0_n_Surfactant"].casadi_var)
                / (m._all_variables["e0_V_Reactor"].casadi_var)
            )
        )
        ** (m._all_variables["e0_P_Surfactant"].casadi_var)
        * (
            1.0
            + (m._all_variables["e0_k_LM_Hyfo"].casadi_var)
            / (
                1.0
                + ca.exp(
                    -(
                        m._all_variables["e0_K_LM"].casadi_var
                        - (m._all_variables["e0_n_Lig"].casadi_var)
                        / (m._all_variables["e0_n_Cat"].casadi_var)
                    )
                    * m._all_variables["e0_P_trig_Hyfo"].casadi_var
                )
            )
        )
        * m._all_variables["e0_greek_psi_cat"].casadi_var
        * m._all_variables["e0_k_ref_r6"].casadi_var
        * ca.exp(
            -(m._all_variables["e0_E_r6"].casadi_var)
            / (m._all_variables["e0_R"].casadi_var)
            * (
                (1.0) / (m._all_variables["e0_T"].casadi_var)
                - (1.0) / (m._all_variables["e0_T_ref"].casadi_var)
            )
        )
        * m._all_variables["e0_c_i1"].casadi_var
        * m._all_variables["e0_c_i6"].casadi_var
        * m._all_variables["e0_c_i7"].casadi_var
    )
    dydx29 = m._all_variables["e0_r_i1"].casadi_var - (
        -m._all_variables["e0_r_r1"].casadi_var
        - m._all_variables["e0_r_r3"].casadi_var
        - m._all_variables["e0_r_r5"].casadi_var
        - m._all_variables["e0_r_r6"].casadi_var
    )
    dydx30 = m._all_variables["e0_r_i2"].casadi_var - (
        m._all_variables["e0_r_r1"].casadi_var
        - m._all_variables["e0_r_r2"].casadi_var
        - m._all_variables["e0_r_r4"].casadi_var
    )
    dydx31 = m._all_variables["e0_r_i3"].casadi_var - (
        m._all_variables["e0_r_r4"].casadi_var + m._all_variables["e0_r_r6"].casadi_var
    )
    dydx32 = m._all_variables["e0_r_i4"].casadi_var - (
        m._all_variables["e0_r_r2"].casadi_var + m._all_variables["e0_r_r3"].casadi_var
    )
    dydx33 = m._all_variables["e0_r_i5"].casadi_var - (
        m._all_variables["e0_r_r5"].casadi_var
    )

    m.add_differential_equations(
        [dydx1, dydx2, dydx3, dydx4, dydx5,]
    )
    m.add_algebraic_equations(
        [
            dydx6,
            dydx7,
            dydx8,
            dydx9,
            dydx10,
            dydx11,
            dydx12,
            dydx13,
            dydx14,
            dydx15,
            dydx16,
            dydx17,
            dydx18,
            dydx19,
            dydx20,
            dydx21,
            dydx22,
            dydx23,
            dydx24,
            dydx25,
            dydx26,
            dydx27,
            dydx28,
            dydx29,
            dydx30,
            dydx31,
            dydx32,
            dydx33,
        ]
    )

    return variable_list, m


if __name__ == "__main__":

    variable_list, m = initialize_problem()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(0, 10000, 11)
    # time_grid = np.insert(time_grid, 0, 0)
    # time_grid = np.array([0, 1000])

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
    r = sim_fixed.analyze()

    print(r)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data()

    with open("outfile.csv", "w") as outfile:
        outfile.write("")
    for name, r in res.items():
        for num, rr in enumerate(r.value.value):
            with open("outfile.csv", "a") as outfile:
                outfile.write(f"{name}, {r.value.time[num]} , {rr} \n")

        with open("outfile.csv", "a") as outfile:
            outfile.write(f"{name}, {r.value.time[10]} , {r.value.value[10]} \n")


    # res.plot_states()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="	")
