from __future__ import annotations

import copy
import numpy as np

import casadi as ca
import mopeds

SUPPORTED_COMPOUNDS = [
    "1_propanol,water",
    "1_propanol,acetic_acid",
    "1_propanol,propyl_acetate",
    "water,acetic_acid",
    "propyl_acetate,acetic_acid",
    "water,propyl_acetate",
    "methanol,water",
    "phenol,butyl_acetate",
]


def get_model_e1_1(
    mode: str = "Px",
    compounds: str = "1_propanol,water",
    antoine: bool = False,
    dippr: bool = False,
) -> tuple[mopeds.VariableList, mopeds.Model]:  # noqa: C901
    """Model E1.1.
    - Mode "Px" - e0_P and e0_x_L_i1 are control variables. Used for parameter estimation.
    - Mode "PT" - e0_P and e0_T are control variables. Used to plot Pxy diagram.
    """
    available_mods = ["PT", "Px"]

    if mode not in available_mods:
        raise TypeError("Mode is not supported")

    if antoine and dippr:
        raise NotImplementedError

    if antoine:
        if not (
            compounds == "methanol,water" or compounds == "1_propanol,propyl_acetate"
        ):
            antoine = False
            print("antoine flag is ignored")
    elif dippr is False:
        if compounds == "phenol,butyl_acetate":
            raise NotImplementedError

    if compounds not in SUPPORTED_COMPOUNDS:
        raise NotImplementedError

    variable_list = mopeds.VariableList()
    # fmt:off

    A_NRTL_lb = -20
    A_NRTL_ub = 20
    B_NRTL_lb = -12000
    B_NRTL_ub = 12000
    C_NRTL_lb = 0.1
    C_NRTL_ub = 0.6

    variable_list.add_variable(mopeds.VariableConstant("e0_A_NRTL_i1_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_NRTL_i1_j2", -1.7411, A_NRTL_lb, A_NRTL_ub))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_w25_i1", -8.26933))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_NRTL_i2_j1", 5.4486, A_NRTL_lb, A_NRTL_ub))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_NRTL_i2_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_w25_i2", -7.90299))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_NRTL_i1_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_B_NRTL_i1_j2", 576.446, B_NRTL_lb, B_NRTL_ub))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_w25_i1", 1.28113))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_B_NRTL_i2_j1", -861.179, B_NRTL_lb, B_NRTL_ub))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_NRTL_i2_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_w25_i2", 2.01149))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_NRTL_i1_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_C_NRTL_i1_j2", 0.3, C_NRTL_lb, C_NRTL_ub))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_w25_i1", -6.84101))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_C_NRTL_i2_j1", 0.3, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_NRTL_i2_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_w25_i2", -2.46727))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_D_NRTL_i1_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_D_NRTL_i1_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_D_w25_i1", 1.36585))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_D_NRTL_i2_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_D_NRTL_i2_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_D_w25_i2", -1.85798))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_E_NRTL_i1_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_E_NRTL_i1_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_E_NRTL_i2_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_E_NRTL_i2_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_F_NRTL_i1_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_F_NRTL_i1_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_F_NRTL_i2_j1", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_F_NRTL_i2_j2", 0.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_w25_i1", 15.4515))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_P_w25_i2", 16.9098))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_T_cr_i1", 536.765))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_T_cr_i2", 647.108))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_gamma_i1", 1.1037020390316847, 1e-20, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_gamma_i2", 1.9917551738165233, 1e-20, 1.0E9))  # noqa: E501
    if mode == "Px":
        variable_list.add_variable(mopeds.VariableControl("e0_x_L_i1", 0.7254665499859735, 1e-8, 0.99999))  # noqa: E501
    elif mode == "PT":
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_L_i1", 0.7254665499859735, 1e-8, 0.99999))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_L_i2", 0.27453345001402657, 1e-8, 0.99999))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_f_V_i1", 61677.74932881687, 1e3, 1.0E6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_f_V_i2", 38322.25067118313, 1e3, 1.0E6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_V_i1", 0.6167774932881687, 1e-8, 0.99999))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_V_i2", 0.3832225067118313, 1e-8, 0.99999))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_alpha_i1_j1", 0.0, 0, 0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_alpha_i1_j2", 0.3, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_alpha_i2_j1", 0.3, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_alpha_i2_j2", 0.0, 0, 0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_tau_i1_j1", 0.0, 0, 0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_tau_i1_j2", -0.15375041993666527, -50, 50))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_tau_i2_j1", 3.07718598375327, -50, 50))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_tau_i2_j2", 0.0, -0.001, 0.001))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_G_i1_j1", 1.0, 1, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_G_i1_j2", 1.0472054353668958, 1e-20, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_G_i2_j1", 0.3972633789853145, 1e-20, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_G_i2_j2", 1.0, 1, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P_s_o_i1", 77029.89041499408, 1e4, 1e6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P_s_o_i2", 70084.14598001179, 1e4, 1e6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_f_L_i1", 61677.74932881687, 1e4, 1e6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_f_L_i2", 38322.25067118313, 1e4, 1e6))  # noqa: E501

    variable_list.add_variable(mopeds.VariableControl("e0_greek_phiv_i1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_phiv_i2", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_P", 100000.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_phiv_s_o_i1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_greek_phiv_s_o_i2", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    if mode == "Px":
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_T", 363.15, 300, 500))  # noqa: E501
    elif mode == "PT":
        variable_list.add_variable(mopeds.VariableControl("e0_T", 363.15, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_Pe_o_i1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_Pe_o_i2", 1.0, -1.0E9, 1.0E9))  # noqa: E501

    data_1_propanol = {
        "e0_A_w25": -8.26933,
        "e0_B_w25": 1.28113,
        "e0_C_w25": -6.84101,
        "e0_D_w25": 1.36585,
        "e0_P_w25": 15.4515,
        "e0_T_cr": 536.765,
    }

    # own data
    data_1_propanol_antoine = {
        "e0_A_w25": 4.65413,
        "e0_B_w25": 1292.869,
        "e0_C_w25": -91.992,
        "e0_D_w25": 0,
        "e0_P_w25": 0,
        "e0_T_cr": 0,
    }

    # 10.1016/j.fluid.2017.05.012
    # data_1_propanol_antoine = {
    #     'e0_A_w25': 4.87601,
    #     'e0_B_w25': 1441.629,
    #     'e0_C_w25': -74.299,
    #     "e0_D_w25": 0,
    #     "e0_P_w25": 0,
    #     "e0_T_cr": 0,
    # }

    data_water = {
        "e0_A_w25": -7.90299,
        "e0_B_w25": 2.01149,
        "e0_C_w25": -2.46727,
        "e0_D_w25": -1.85798,
        "e0_P_w25": 16.9098,
        "e0_T_cr": 647.108,
    }

    data_acetic_acid = {
        "e0_A_w25": -8.47304,
        "e0_B_w25": 1.47485,
        "e0_C_w25": -0.702165,
        "e0_D_w25": -5.64878,
        "e0_P_w25": 15.5709,
        "e0_T_cr": 592.998,
    }

    data_propyl_acetate = {
        "e0_A_w25": -8.63721,
        "e0_B_w25": 3.80217,
        "e0_C_w25": -5.46836,
        "e0_D_w25": -2.20478,
        "e0_P_w25": 15.0363,
        "e0_T_cr": 549.691,
    }

    # own data
    data_propyl_acetate_antoine = {
        "e0_A_w25": 3.84871,
        "e0_B_w25": 1088.392,
        "e0_C_w25": -90.571,
        "e0_D_w25": 0,
        "e0_P_w25": 0,
        "e0_T_cr": 0,
    }

    # 10.1016/j.fluid.2017.05.012
    # data_propyl_acetate_antoine = {
    #     'e0_A_w25': 4.14386,
    #     'e0_B_w25': 1283.861,
    #     'e0_C_w25': -64.378,
    #     "e0_D_w25": 0,
    #     "e0_P_w25": 0,
    #     "e0_T_cr": 0,
    # }

    data_methanol = {
        "e0_A_w25": -8.68171,
        "e0_B_w25": 1.42196,
        "e0_C_w25": -2.93922,
        "e0_D_w25": -0.348628,
        "e0_P_w25": 15.9016,
        "e0_T_cr": 512.68,
    }

    data_methanol_antoine = {
        "e0_A_w25": 82.718,
        "e0_B_w25": -6904.5,
        "e0_C_w25": -8.8622,
        "e0_D_w25": 0,
        "e0_P_w25": 7.4664e-6,
        "e0_T_cr": 0,
    }

    data_water_antoine = {
        "e0_A_w25": 73.649,
        "e0_B_w25": -7258.2,
        "e0_C_w25": -7.3037,
        "e0_D_w25": 0,
        "e0_P_w25": 4.1653e-6,
        "e0_T_cr": 0,
    }

    data_phenol_dippr = {
        "e0_A_w25": 95.444,
        "e0_B_w25": -10113,
        "e0_C_w25": -10.09,
        "e0_D_w25": 6.76e-18,
        "e0_P_w25": 6,
        "e0_T_cr": 0,
    }

    data_butyl_acetate_dippr = {
        "e0_A_w25": 122.82,
        "e0_B_w25": -9253.2,
        "e0_C_w25": -14.99,
        "e0_D_w25": 1.05e-5,
        "e0_P_w25": 2,
        "e0_T_cr": 0,
    }

    # https://pubs.acs.org/doi/pdf/10.1021/acs.jced.5b01015
    if compounds == "1_propanol,water":
        # data = {
        #     "e0_A_NRTL_i1_j2": 1648.8,
        #     "e0_A_NRTL_i2_j1": 7896.7,
        #     # "e0_B_NRTL_i1_j2": 799.35,
        #     # "e0_B_NRTL_i2_j1": -238.29,
        #     "e0_C_NRTL_i1_j2": 0.477,
        # }

        # 10.1021/je950237o
        #         g12 - g22 ) 163.680 K
        # g21 - g11 ) 912.580 K
        # alpha 0.4548

        # 10.1021/je990069q
        # Ri,j ∆g i,j/(J‚mol-1 ) ∆g j,i/(J‚mol-1 ) ∆y a ∆T b 1-propanol water 0.477 1 648.8 7 896.7

        # 10.1021/je9601467 (parameters for different pressures)

        # 10.1016/j.jct.2016.10.001
        data = {
            "e0_A_NRTL_i1_j2": -1.7387,
            "e0_A_NRTL_i2_j1": 3.2932 ,
            "e0_B_NRTL_i1_j2": 799.35,
            "e0_B_NRTL_i2_j1": -238.29,
            "e0_C_NRTL_i1_j2": 0.47,
        }
        mapping = {"_i1": data_1_propanol, "_i2": data_water}

        guess = {'e0_greek_gamma_i1': 1.9116419142498342, 'e0_greek_gamma_i2': 1.3520187168265763, 'e0_x_L_i2': 0.7, 'e0_f_V_i1': 39862.98213733318, 'e0_f_V_i2': 60137.01786266681, 'e0_x_V_i1': 0.39862982137333186, 'e0_x_V_i2': 0.6013701786266681, 'e0_greek_alpha_i1_j1': 5.912534467629902e-35, 'e0_greek_alpha_i1_j2': 0.44999981013450896, 'e0_greek_alpha_i2_j1': 0.4499999969994455, 'e0_greek_alpha_i2_j2': -3.76158192263132e-37, 'e0_greek_tau_i1_j1': 6.509643212028857e-30, 'e0_greek_tau_i1_j2': 0.42196832497944764, 'e0_greek_tau_i2_j1': 2.614874578832609, 'e0_greek_tau_i2_j2': 8.192833761050387e-30, 'e0_G_i1_j1': 1.0, 'e0_G_i1_j2': 0.827053688792138, 'e0_G_i2_j1': 0.3082964211083324, 'e0_G_i2_j2': 1.0, 'e0_P_s_o_i1': 69509.15133945853, 'e0_P_s_o_i2': 63542.03861891888, 'e0_f_L_i1': 39862.98213733318, 'e0_f_L_i2': 60137.01786266681, 'e0_T': 360.59089103398713}  # noqa:E501

    elif compounds == "1_propanol,acetic_acid":
        data = {
            "e0_A_NRTL_i1_j2": 1.72661,
            "e0_A_NRTL_i2_j1": -1.3422,
            "e0_B_NRTL_i1_j2": -14.741,
            "e0_B_NRTL_i2_j1": 13.0943,
            "e0_C_NRTL_i1_j2": 0.5,
        }
        mapping = {"_i1": data_1_propanol, "_i2": data_acetic_acid}

        guess = {}

    elif compounds == "1_propanol,propyl_acetate":
        data = {
            "e0_A_NRTL_i1_j2": 0.860697,
            "e0_A_NRTL_i2_j1": -0.226867,
            "e0_B_NRTL_i1_j2": -1.15726,
            "e0_B_NRTL_i2_j1": 7.1614,
            "e0_C_NRTL_i1_j2": 0.1,
        }

        if antoine:
            mapping = {"_i1": data_1_propanol_antoine, "_i2": data_propyl_acetate_antoine}
        else:
            mapping = {"_i1": data_1_propanol, "_i2": data_propyl_acetate}

        guess = {'e0_greek_gamma_i1': 1.353241274690732, 'e0_greek_gamma_i2': 1.051092623897063, 'e0_x_L_i2': 0.7, 'e0_f_V_i1': 38888.25980164735, 'e0_f_V_i2': 62431.74019835266, 'e0_x_V_i1': 0.3838162238615016, 'e0_x_V_i2': 0.6161837761384984, 'e0_greek_alpha_i1_j1': -1.1402042850996157e-35, 'e0_greek_alpha_i1_j2': 0.1, 'e0_greek_alpha_i2_j1': 0.1, 'e0_greek_alpha_i2_j2': 1.3598824722637043e-36, 'e0_greek_tau_i1_j1': -2.8166725436663324e-33, 'e0_greek_tau_i1_j2': 0.8575585588962009, 'e0_greek_tau_i2_j1': -0.20744558033566593, 'e0_greek_tau_i2_j2': -1.2574590005112332e-33, 'e0_G_i1_j1': 1.0, 'e0_G_i1_j2': 0.9178182837968675, 'e0_G_i2_j1': 1.0209612219833315, 'e0_G_i2_j2': 1.0, 'e0_P_s_o_i1': 95790.40714803999, 'e0_P_s_o_i2': 84852.84574891608, 'e0_f_L_i1': 38888.25980164735, 'e0_f_L_i2': 62431.74019835266, 'e0_T': 368.73720478587626}  # noqa:E501

    elif compounds == "water,acetic_acid":
        data = {
            "e0_A_NRTL_i1_j2": 3.3293,
            "e0_A_NRTL_i2_j1": -1.9763,
            "e0_B_NRTL_i1_j2": -723.888,
            "e0_B_NRTL_i2_j1": 609.889,
            "e0_C_NRTL_i1_j2": 0.3,
        }

        mapping = {"_i1": data_water, "_i2": data_acetic_acid}

        guess = {'e0_greek_gamma_i1': 1.4550329121681969, 'e0_greek_gamma_i2': 1.037299873134237, 'e0_x_L_i2': 0.7, 'e0_f_V_i1': 52525.66080176381, 'e0_f_V_i2': 48804.339198236186, 'e0_x_V_i1': 0.5183623882538617, 'e0_x_V_i2': 0.48163761174613823, 'e0_greek_alpha_i1_j1': 0.0, 'e0_greek_alpha_i1_j2': 0.3, 'e0_greek_alpha_i2_j1': 0.3, 'e0_greek_alpha_i2_j2': 4.814824864399629e-35, 'e0_greek_tau_i1_j1': 1.398779314498214e-32, 'e0_greek_tau_i1_j2': 1.4144988583293105, 'e0_greek_tau_i2_j1': -0.3630446987760605, 'e0_greek_tau_i2_j2': 1.396299209680746e-32, 'e0_G_i1_j1': 1.0, 'e0_G_i1_j2': 0.6541950976915494, 'e0_G_i2_j1': 1.115065792204627, 'e0_G_i2_j2': 1.0, 'e0_P_s_o_i1': 120330.97982985013, 'e0_P_s_o_i2': 67213.43207942945, 'e0_f_L_i1': 52525.66080176381, 'e0_f_L_i2': 48804.339198236186, 'e0_T': 378.04865698398226}  # noqa:E501

    elif compounds == "propyl_acetate,acetic_acid":
        data = {
            "e0_A_NRTL_i1_j2": -0.245896,
            "e0_A_NRTL_i2_j1": 0.415517,
            "e0_B_NRTL_i1_j2": 2.17433,
            "e0_B_NRTL_i2_j1": 1.12531,
            "e0_C_NRTL_i1_j2": 0.1,
        }

        mapping = {"_i1": data_propyl_acetate, "_i2": data_acetic_acid}

        guess = {}

    elif compounds == "water,propyl_acetate":
        data = {
            "e0_A_NRTL_i1_j2": 14.1308,
            "e0_A_NRTL_i2_j1": -6.57523,
            "e0_B_NRTL_i1_j2": -2298.07,
            "e0_B_NRTL_i2_j1": 1730.49,
            "e0_C_NRTL_i1_j2": 0.1,
        }

        mapping = {"_i1": data_water, "_i2": data_propyl_acetate}

        guess = {'e0_greek_gamma_i1': 2.476438948012118, 'e0_greek_gamma_i2': 1.9806906932218222, 'e0_x_L_i2': 0.4, 'e0_f_V_i1': 64476.92919200069, 'e0_f_V_i2': 36848.07080799932, 'e0_x_V_i1': 0.6363378158598636, 'e0_x_V_i2': 0.3636621841401364, 'e0_greek_alpha_i1_j1': 0.0, 'e0_greek_alpha_i1_j2': 0.1, 'e0_greek_alpha_i2_j1': 0.1, 'e0_greek_alpha_i2_j2': 1.2753613028178902e-34, 'e0_greek_tau_i1_j1': -2.2710565904214285e-30, 'e0_greek_tau_i1_j2': 7.583790630451274, 'e0_greek_tau_i2_j1': -1.6452086151377565, 'e0_greek_tau_i2_j2': 8.628166150854817e-32, 'e0_G_i1_j1': 1.0, 'e0_G_i1_j2': 0.4684250995180417, 'e0_G_i2_j1': 1.1788281614350231, 'e0_G_i2_j2': 1.0, 'e0_P_s_o_i1': 43393.57880782637, 'e0_P_s_o_i2': 46509.117922977755, 'e0_f_L_i1': 64476.92919200069, 'e0_f_L_i2': 36848.07080799932, 'e0_T': 351.01064780641997}  # noqa:E501

    elif compounds == "methanol,water":
        data = {
            "e0_A_NRTL_i1_j2": -2.96824,
            "e0_A_NRTL_i2_j1": 6.11368,
            "e0_B_NRTL_i1_j2": 509.265,
            "e0_B_NRTL_i2_j1": -1232.35,
            "e0_C_NRTL_i1_j2": 0.1,
        }

        if antoine:
            mapping = {"_i1": data_methanol_antoine, "_i2": data_water_antoine}
        else:
            mapping = {"_i1": data_methanol, "_i2": data_water}

        guess = {'e0_greek_gamma_i1': 1.3353422992879826, 'e0_greek_gamma_i2': 1.0922475962652898, 'e0_x_L_i2': 0.7, 'e0_f_V_i1': 67771.3763043795, 'e0_f_V_i2': 33528.6236956205, 'e0_x_V_i1': 0.6690165479208243, 'e0_x_V_i2': 0.3309834520791758, 'e0_greek_alpha_i1_j1': 0.0, 'e0_greek_alpha_i1_j2': 0.1, 'e0_greek_alpha_i2_j1': 0.1, 'e0_greek_alpha_i2_j2': 0.0, 'e0_greek_tau_i1_j1': 2.164602317378193e-33, 'e0_greek_tau_i1_j2': -1.518446413085901, 'e0_greek_tau_i2_j1': 2.6053824990258705, 'e0_greek_tau_i2_j2': -1.4060793226795874e-33, 'e0_G_i1_j1': 1.0, 'e0_G_i1_j2': 1.1639793880697349, 'e0_G_i2_j1': 0.7706366790343989, 'e0_G_i2_j2': 1.0, 'e0_P_s_o_i1': 169173.54284494658, 'e0_P_s_o_i2': 43852.72534786404, 'e0_f_L_i1': 67771.3763043795, 'e0_f_L_i2': 33528.6236956205, 'e0_T': 351.2672456249278}  # noqa: E501

    elif compounds == "phenol,butyl_acetate":
        data = {
            "e0_A_NRTL_i1_j2": 1.9532,
            "e0_A_NRTL_i2_j1": 0.0854,
            "e0_B_NRTL_i1_j2": -873.997,
            "e0_B_NRTL_i2_j1": -302.802,
            "e0_C_NRTL_i1_j2": 0.3,
        }

        if dippr:
            mapping = {"_i1": data_phenol_dippr, "_i2": data_butyl_acetate_dippr}
        else:
            raise NotImplementedError

        guess = {'e0_C_NRTL_i2_j1': 0.3, 'e0_greek_gamma_i1': 0.9405572473707309, 'e0_greek_gamma_i2': 0.6974943478887967, 'e0_x_L_i2': 0.2745334500140265, 'e0_f_V_i1': 44531.668748590906, 'e0_f_V_i2': 55468.331251409094, 'e0_x_V_i1': 0.4453166874859091, 'e0_x_V_i2': 0.554683312514091, 'e0_greek_alpha_i1_j1': 2.177322440987124e-43, 'e0_greek_alpha_i1_j2': 0.3, 'e0_greek_alpha_i2_j1': 0.3, 'e0_greek_alpha_i2_j2': 2.1773461751377496e-43, 'e0_greek_tau_i1_j1': 2.617442297629857e-32, 'e0_greek_tau_i1_j2': -0.03424494666986845, 'e0_greek_tau_i2_j1': -0.6031633528965541, 'e0_greek_tau_i2_j2': -1.6971008373735278e-31, 'e0_G_i1_j1': 1.0, 'e0_G_i1_j2': 1.0103264374210694, 'e0_G_i2_j1': 1.1983540687106677, 'e0_G_i2_j2': 1.0, 'e0_P_s_o_i1': 65262.89624152961, 'e0_P_s_o_i2': 289673.74892955023, 'e0_f_L_i1': 44531.668748590906, 'e0_f_L_i2': 55468.331251409094, 'e0_T': 439.7590994731476}
    for suffix, pure_data in mapping.items():
        for key, value in pure_data.items():
            data[key+suffix] = value

    for var_name, value in data.items():
        variable = variable_list[var_name]
        variable.dataframe.iloc[0] = value
        # breakpoint()
        if isinstance(variable, mopeds.VariableParameter):
            variable_list[var_name].guess = value

    for var_name, value in guess.items():
        variable_list[var_name].guess = value

    m = mopeds.Model(variable_list)

    e0_greek_phiv_i1 = m.varlist_all["e0_greek_phiv_i1"].casadi_var  # noqa: E501
    e0_greek_phiv_i2 = m.varlist_all["e0_greek_phiv_i2"].casadi_var  # noqa: E501
    e0_P = m.varlist_all["e0_P"].casadi_var  # noqa: E501
    e0_greek_phiv_s_o_i1 = m.varlist_all["e0_greek_phiv_s_o_i1"].casadi_var  # noqa: E501
    e0_greek_phiv_s_o_i2 = m.varlist_all["e0_greek_phiv_s_o_i2"].casadi_var  # noqa: E501
    e0_T = m.varlist_all["e0_T"].casadi_var  # noqa: E501
    e0_Pe_o_i1 = m.varlist_all["e0_Pe_o_i1"].casadi_var  # noqa: E501
    e0_Pe_o_i2 = m.varlist_all["e0_Pe_o_i2"].casadi_var  # noqa: E501
    e0_A_NRTL_i1_j1 = m.varlist_all["e0_A_NRTL_i1_j1"].casadi_var  # noqa: E501
    e0_A_NRTL_i1_j2 = m.varlist_all["e0_A_NRTL_i1_j2"].casadi_var  # noqa: E501
    e0_A_w25_i1 = m.varlist_all["e0_A_w25_i1"].casadi_var  # noqa: E501
    e0_A_NRTL_i2_j1 = m.varlist_all["e0_A_NRTL_i2_j1"].casadi_var  # noqa: E501
    e0_A_NRTL_i2_j2 = m.varlist_all["e0_A_NRTL_i2_j2"].casadi_var  # noqa: E501
    e0_A_w25_i2 = m.varlist_all["e0_A_w25_i2"].casadi_var  # noqa: E501
    e0_B_NRTL_i1_j1 = m.varlist_all["e0_B_NRTL_i1_j1"].casadi_var  # noqa: E501
    e0_B_NRTL_i1_j2 = m.varlist_all["e0_B_NRTL_i1_j2"].casadi_var  # noqa: E501
    e0_B_w25_i1 = m.varlist_all["e0_B_w25_i1"].casadi_var  # noqa: E501
    e0_B_NRTL_i2_j1 = m.varlist_all["e0_B_NRTL_i2_j1"].casadi_var  # noqa: E501
    e0_B_NRTL_i2_j2 = m.varlist_all["e0_B_NRTL_i2_j2"].casadi_var  # noqa: E501
    e0_B_w25_i2 = m.varlist_all["e0_B_w25_i2"].casadi_var  # noqa: E501
    e0_C_NRTL_i1_j1 = m.varlist_all["e0_C_NRTL_i1_j1"].casadi_var  # noqa: E501
    e0_C_NRTL_i1_j2 = m.varlist_all["e0_C_NRTL_i1_j2"].casadi_var  # noqa: E501
    e0_C_w25_i1 = m.varlist_all["e0_C_w25_i1"].casadi_var  # noqa: E501
    e0_C_NRTL_i2_j1 = m.varlist_all["e0_C_NRTL_i2_j1"].casadi_var  # noqa: E501
    e0_C_NRTL_i2_j2 = m.varlist_all["e0_C_NRTL_i2_j2"].casadi_var  # noqa: E501
    e0_C_w25_i2 = m.varlist_all["e0_C_w25_i2"].casadi_var  # noqa: E501
    e0_D_NRTL_i1_j1 = m.varlist_all["e0_D_NRTL_i1_j1"].casadi_var  # noqa: E501
    e0_D_NRTL_i1_j2 = m.varlist_all["e0_D_NRTL_i1_j2"].casadi_var  # noqa: E501
    e0_D_w25_i1 = m.varlist_all["e0_D_w25_i1"].casadi_var  # noqa: E501
    e0_D_NRTL_i2_j1 = m.varlist_all["e0_D_NRTL_i2_j1"].casadi_var  # noqa: E501
    e0_D_NRTL_i2_j2 = m.varlist_all["e0_D_NRTL_i2_j2"].casadi_var  # noqa: E501
    e0_D_w25_i2 = m.varlist_all["e0_D_w25_i2"].casadi_var  # noqa: E501
    e0_E_NRTL_i1_j1 = m.varlist_all["e0_E_NRTL_i1_j1"].casadi_var  # noqa: E501
    e0_E_NRTL_i1_j2 = m.varlist_all["e0_E_NRTL_i1_j2"].casadi_var  # noqa: E501
    e0_E_NRTL_i2_j1 = m.varlist_all["e0_E_NRTL_i2_j1"].casadi_var  # noqa: E501
    e0_E_NRTL_i2_j2 = m.varlist_all["e0_E_NRTL_i2_j2"].casadi_var  # noqa: E501
    e0_F_NRTL_i1_j1 = m.varlist_all["e0_F_NRTL_i1_j1"].casadi_var  # noqa: E501
    e0_F_NRTL_i1_j2 = m.varlist_all["e0_F_NRTL_i1_j2"].casadi_var  # noqa: E501
    e0_F_NRTL_i2_j1 = m.varlist_all["e0_F_NRTL_i2_j1"].casadi_var  # noqa: E501
    e0_F_NRTL_i2_j2 = m.varlist_all["e0_F_NRTL_i2_j2"].casadi_var  # noqa: E501
    e0_P_w25_i1 = m.varlist_all["e0_P_w25_i1"].casadi_var  # noqa: E501
    e0_P_w25_i2 = m.varlist_all["e0_P_w25_i2"].casadi_var  # noqa: E501
    e0_T_cr_i1 = m.varlist_all["e0_T_cr_i1"].casadi_var  # noqa: E501
    e0_T_cr_i2 = m.varlist_all["e0_T_cr_i2"].casadi_var  # noqa: E501
    e0_greek_gamma_i1 = m.varlist_all["e0_greek_gamma_i1"].casadi_var  # noqa: E501
    e0_greek_gamma_i2 = m.varlist_all["e0_greek_gamma_i2"].casadi_var  # noqa: E501
    e0_x_L_i1 = m.varlist_all["e0_x_L_i1"].casadi_var  # noqa: E501
    e0_x_L_i2 = m.varlist_all["e0_x_L_i2"].casadi_var  # noqa: E501
    e0_f_V_i1 = m.varlist_all["e0_f_V_i1"].casadi_var  # noqa: E501
    e0_f_V_i2 = m.varlist_all["e0_f_V_i2"].casadi_var  # noqa: E501
    e0_x_V_i1 = m.varlist_all["e0_x_V_i1"].casadi_var  # noqa: E501
    e0_x_V_i2 = m.varlist_all["e0_x_V_i2"].casadi_var  # noqa: E501
    e0_greek_alpha_i1_j1 = m.varlist_all["e0_greek_alpha_i1_j1"].casadi_var  # noqa: E501
    e0_greek_alpha_i1_j2 = m.varlist_all["e0_greek_alpha_i1_j2"].casadi_var  # noqa: E501
    e0_greek_alpha_i2_j1 = m.varlist_all["e0_greek_alpha_i2_j1"].casadi_var  # noqa: E501
    e0_greek_alpha_i2_j2 = m.varlist_all["e0_greek_alpha_i2_j2"].casadi_var  # noqa: E501
    e0_greek_tau_i1_j1 = m.varlist_all["e0_greek_tau_i1_j1"].casadi_var  # noqa: E501
    e0_greek_tau_i1_j2 = m.varlist_all["e0_greek_tau_i1_j2"].casadi_var  # noqa: E501
    e0_greek_tau_i2_j1 = m.varlist_all["e0_greek_tau_i2_j1"].casadi_var  # noqa: E501
    e0_greek_tau_i2_j2 = m.varlist_all["e0_greek_tau_i2_j2"].casadi_var  # noqa: E501
    e0_G_i1_j1 = m.varlist_all["e0_G_i1_j1"].casadi_var  # noqa: E501
    e0_G_i1_j2 = m.varlist_all["e0_G_i1_j2"].casadi_var  # noqa: E501
    e0_G_i2_j1 = m.varlist_all["e0_G_i2_j1"].casadi_var  # noqa: E501
    e0_G_i2_j2 = m.varlist_all["e0_G_i2_j2"].casadi_var  # noqa: E501
    e0_P_s_o_i1 = m.varlist_all["e0_P_s_o_i1"].casadi_var  # noqa: E501
    e0_P_s_o_i2 = m.varlist_all["e0_P_s_o_i2"].casadi_var  # noqa: E501
    e0_f_L_i1 = m.varlist_all["e0_f_L_i1"].casadi_var  # noqa: E501
    e0_f_L_i2 = m.varlist_all["e0_f_L_i2"].casadi_var  # noqa: E501

    EQ_alg1 = (e0_f_L_i1-((e0_x_L_i1*(e0_greek_gamma_i1*(e0_P_s_o_i1*(e0_greek_phiv_s_o_i1*e0_Pe_o_i1))))))  # noqa: E501,E226
    EQ_alg2 = (e0_f_L_i2-((e0_x_L_i2*(e0_greek_gamma_i2*(e0_P_s_o_i2*(e0_greek_phiv_s_o_i2*e0_Pe_o_i2))))))  # noqa: E501,E226
    EQ_alg3 = (e0_f_V_i1-((e0_x_V_i1*(e0_greek_phiv_i1*e0_P))))  # noqa: E501,E226
    EQ_alg4 = (e0_f_V_i2-((e0_x_V_i2*(e0_greek_phiv_i2*e0_P))))  # noqa: E501,E226
    EQ_alg5 = (e0_f_V_i1-(e0_f_L_i1))  # noqa: E501,E226
    EQ_alg6 = (e0_f_V_i2-(e0_f_L_i2))  # noqa: E501,E226
    EQ_alg7 = (e0_G_i1_j1-(ca.exp((-(e0_greek_alpha_i1_j1*e0_greek_tau_i1_j1)))))  # noqa: E501,E226
    EQ_alg8 = (e0_G_i1_j2-(ca.exp((-(e0_greek_alpha_i1_j2*e0_greek_tau_i1_j2)))))  # noqa: E501,E226
    EQ_alg9 = (e0_G_i2_j1-(ca.exp((-(e0_greek_alpha_i2_j1*e0_greek_tau_i2_j1)))))  # noqa: E501,E226
    EQ_alg10 = (e0_G_i2_j2-(ca.exp((-(e0_greek_alpha_i2_j2*e0_greek_tau_i2_j2)))))  # noqa: E501,E226
    EQ_alg11 = (e0_greek_tau_i1_j1-((((e0_A_NRTL_i1_j1+(e0_B_NRTL_i1_j1/e0_T))+(e0_E_NRTL_i1_j1*ca.log(e0_T)))+(e0_F_NRTL_i1_j1*e0_T))))  # noqa: E501,E226
    EQ_alg12 = (e0_greek_tau_i1_j2-((((e0_A_NRTL_i1_j2+(e0_B_NRTL_i1_j2/e0_T))+(e0_E_NRTL_i1_j2*ca.log(e0_T)))+(e0_F_NRTL_i1_j2*e0_T))))  # noqa: E501,E226
    EQ_alg13 = (e0_greek_tau_i2_j1-((((e0_A_NRTL_i2_j1+(e0_B_NRTL_i2_j1/e0_T))+(e0_E_NRTL_i2_j1*ca.log(e0_T)))+(e0_F_NRTL_i2_j1*e0_T))))  # noqa: E501,E226
    EQ_alg14 = (e0_greek_tau_i2_j2-((((e0_A_NRTL_i2_j2+(e0_B_NRTL_i2_j2/e0_T))+(e0_E_NRTL_i2_j2*ca.log(e0_T)))+(e0_F_NRTL_i2_j2*e0_T))))  # noqa: E501,E226
    EQ_alg15 = (e0_greek_gamma_i1-(ca.exp((((((e0_x_L_i1*(e0_greek_tau_i1_j1*e0_G_i1_j1))+(e0_x_L_i2*(e0_greek_tau_i2_j1*e0_G_i2_j1))))/(((e0_x_L_i1*e0_G_i1_j1)+(e0_x_L_i2*e0_G_i2_j1))))+(((((e0_x_L_i1*e0_G_i1_j1)/(((e0_x_L_i1*e0_G_i1_j1)+(e0_x_L_i2*e0_G_i2_j1))))*((e0_greek_tau_i1_j1-((((e0_x_L_i1*(e0_G_i1_j1*e0_greek_tau_i1_j1))+(e0_x_L_i2*(e0_G_i2_j1*e0_greek_tau_i2_j1))))/(((e0_x_L_i1*e0_G_i1_j1)+(e0_x_L_i2*e0_G_i2_j1)))))))+(((e0_x_L_i2*e0_G_i1_j2)/(((e0_x_L_i1*e0_G_i1_j2)+(e0_x_L_i2*e0_G_i2_j2))))*((e0_greek_tau_i1_j2-((((e0_x_L_i1*(e0_G_i1_j2*e0_greek_tau_i1_j2))+(e0_x_L_i2*(e0_G_i2_j2*e0_greek_tau_i2_j2))))/(((e0_x_L_i1*e0_G_i1_j2)+(e0_x_L_i2*e0_G_i2_j2)))))))))))))  # noqa: E501,E226
    EQ_alg16 = (e0_greek_gamma_i2-(ca.exp((((((e0_x_L_i1*(e0_greek_tau_i1_j2*e0_G_i1_j2))+(e0_x_L_i2*(e0_greek_tau_i2_j2*e0_G_i2_j2))))/(((e0_x_L_i1*e0_G_i1_j2)+(e0_x_L_i2*e0_G_i2_j2))))+(((((e0_x_L_i1*e0_G_i2_j1)/(((e0_x_L_i1*e0_G_i1_j1)+(e0_x_L_i2*e0_G_i2_j1))))*((e0_greek_tau_i2_j1-((((e0_x_L_i1*(e0_G_i1_j1*e0_greek_tau_i1_j1))+(e0_x_L_i2*(e0_G_i2_j1*e0_greek_tau_i2_j1))))/(((e0_x_L_i1*e0_G_i1_j1)+(e0_x_L_i2*e0_G_i2_j1)))))))+(((e0_x_L_i2*e0_G_i2_j2)/(((e0_x_L_i1*e0_G_i1_j2)+(e0_x_L_i2*e0_G_i2_j2))))*((e0_greek_tau_i2_j2-((((e0_x_L_i1*(e0_G_i1_j2*e0_greek_tau_i1_j2))+(e0_x_L_i2*(e0_G_i2_j2*e0_greek_tau_i2_j2))))/(((e0_x_L_i1*e0_G_i1_j2)+(e0_x_L_i2*e0_G_i2_j2)))))))))))))  # noqa: E501,E226
    EQ_alg17 = (((e0_x_L_i1+e0_x_L_i2))-(1.0))  # noqa: E501,E226
    EQ_alg18 = (((e0_x_V_i1+e0_x_V_i2))-(1.0))  # noqa: E501,E226
    if antoine:
        EQ_alg19 = e0_P_s_o_i1 - 1e5 * 10**(e0_A_w25_i1 - e0_B_w25_i1 / (e0_T + e0_C_w25_i1))  # noqa: E501,E226
        EQ_alg20 = e0_P_s_o_i2 - 1e5 * 10**(e0_A_w25_i2 - e0_B_w25_i2 / (e0_T + e0_C_w25_i2))  # noqa: E501,E226
    elif dippr:
        EQ_alg19 = e0_P_s_o_i1 - ca.exp(e0_A_w25_i1 + e0_B_w25_i1 / e0_T + e0_C_w25_i1 * ca.log(e0_T) + e0_D_w25_i1 * (e0_T ** e0_P_w25_i1))  # noqa: E501,E226
        EQ_alg20 = e0_P_s_o_i2 - ca.exp(e0_A_w25_i2 + e0_B_w25_i2 / e0_T + e0_C_w25_i2 * ca.log(e0_T) + e0_D_w25_i2 * (e0_T ** e0_P_w25_i2))  # noqa: E501,E226
    else:
        EQ_alg19 = (e0_P_s_o_i1-(ca.exp((e0_P_w25_i1+(((((e0_A_w25_i1*((1.0-(e0_T/e0_T_cr_i1))))+(e0_B_w25_i1*(((1.0-(e0_T/e0_T_cr_i1))))**(1.0*1.5)))+(e0_C_w25_i1*(((1.0-(e0_T/e0_T_cr_i1))))**(1.0*2.5)))+(e0_D_w25_i1*(((1.0-(e0_T/e0_T_cr_i1))))**(1.0*5.0)))/(e0_T/e0_T_cr_i1))))))  # noqa: E501,E226
        EQ_alg20 = (e0_P_s_o_i2-(ca.exp((e0_P_w25_i2+(((((e0_A_w25_i2*((1.0-(e0_T/e0_T_cr_i2))))+(e0_B_w25_i2*(((1.0-(e0_T/e0_T_cr_i2))))**(1.0*1.5)))+(e0_C_w25_i2*(((1.0-(e0_T/e0_T_cr_i2))))**(1.0*2.5)))+(e0_D_w25_i2*(((1.0-(e0_T/e0_T_cr_i2))))**(1.0*5.0)))/(e0_T/e0_T_cr_i2))))))  # noqa: E501,E226
    EQ_alg21 = (e0_greek_alpha_i1_j1-((e0_C_NRTL_i1_j1+(e0_D_NRTL_i1_j1*((e0_T-273.15))))))  # noqa: E501,E226
    EQ_alg22 = (e0_greek_alpha_i1_j2-((e0_C_NRTL_i1_j2+(e0_D_NRTL_i1_j2*((e0_T-273.15))))))  # noqa: E501,E226
    EQ_alg23 = (e0_greek_alpha_i2_j1-((e0_C_NRTL_i2_j1+(e0_D_NRTL_i2_j1*((e0_T-273.15))))))  # noqa: E501,E226
    EQ_alg24 = (e0_greek_alpha_i2_j2-((e0_C_NRTL_i2_j2+(e0_D_NRTL_i2_j2*((e0_T-273.15))))))  # noqa: E501,E226

    EQ_alg101 = e0_C_NRTL_i1_j2 - e0_C_NRTL_i2_j1

    list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6, EQ_alg7, EQ_alg8, EQ_alg9, EQ_alg10, EQ_alg11, EQ_alg12, EQ_alg13, EQ_alg14, EQ_alg15, EQ_alg16, EQ_alg17, EQ_alg18, EQ_alg19, EQ_alg20, EQ_alg21, EQ_alg22, EQ_alg23, EQ_alg24, EQ_alg101]  # noqa: E501

    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    m.compounds = compounds
    m.antoine = antoine
    m.dippr = dippr

    return variable_list, m


varlist, m = get_model_e1_1(compounds="1_propanol,propyl_acetate")
sim = mopeds.SimulatorNLE(m, varlist)
l_vl = []

for P in np.linspace(1e5, 2e5, 2):
    sim.change_independent_variables({"e0_P": P})
    res = sim.simulate()[2]
    vl = copy.deepcopy(varlist)
    vl["e0_x_V_i1"].value = res["e0_x_V_i1"].value[0]
    vl["e0_T"].value = res["e0_T"].value[0]
    vl["e0_P"].value = P
    vl["e0_A_NRTL_i2_j1"].fixed = False
    # vl["e0_A_NRTL_i2_j1"].guess = 3.2932
    vl["e0_A_NRTL_i2_j1"].guess = 3.1932
    vl["e0_A_NRTL_i1_j2"].fixed = False
    l_vl.append(vl)

pe = mopeds.ParameterEstimationNLE(m, l_vl)
pe.prepare_nle()
v = pe.simulate_all_mx
v1 = dict(
    zip(
        pe.varlist_decision.keys(),
        [
            -0.44155,
            7.98915,
            0.47,
            1.1797,
            1.65805,
            0.274533,
            67411.6,
            32588.4,
            0.674116,
            0.325884,
            0,
            0.47,
            0.47,
            0,
            0,
            1.75621,
            7.33399,
            0,
            1,
            0.438028,
            1.08688e-05,
            1,
            78766.8,
            71592.6,
            67411.6,
            32588.4,
            363.712,
            0.47,
            1.16673,
            1.64365,
            0.274533,
            135893,
            64106.5,
            0.679467,
            0.320533,
            0,
            0.47,
            0.47,
            0,
            0,
            1.64601,
            7.36684,
            0,
            1,
            0.461312,
            1.07328e-05,
            1,
            160550,
            142068,
            135893,
            64106.5,
            382.911,
        ],
    )
)
# breakpoint()
# pe.solver_settings["verbose_init"] = True
# pe.solver_settings["verbose"] = True
# pe.solver_settings["expand"] = True
pe.solver_settings["show_eval_warnings"] = False
# pe.solver_settings["ipopt"]["linear_solver"] = "ma57"

res = pe.optimize(
    False,
)
print(res)
breakpoint()
