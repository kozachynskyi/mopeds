import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est


def initialize_problem():

    variable_list = par_est.VariableList()

    variable_list.add_variable(par_est.VariableState("e0_HU_i1", 3.006))
    variable_list.add_variable(par_est.VariableState("e0_HU_i2", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_HU_i3", 0.038722374948782))
    variable_list.add_variable(par_est.VariableState("e0_HU_i4", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_HU_i5", 0.0))
    variable_list.add_variable(par_est.VariableState("e0_HU_i9", 11.82565))
    variable_list.add_variable(par_est.VariableState("e0_HU_i13", 1.300000000000000))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i1", 0.017095627))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i2", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i3", 6.39e-4))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i4", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i5", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i9", 0.65642333))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i12", 0.003910874))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i13", 0.017637407))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i14", 0.05283861))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i7", 0.001479374))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i6", 8.13e-4))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_HU_i6", 0.0016381))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_HU_i7", 0.041437257))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i10", 5.03e-5))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_n_i11", 2.02e-4))
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i1", 772.3244306654833)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i2", 776.3104270331723)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i3", 1090.546648999151)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i4", 910.3398759323168)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i5", 914.3602790411466)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i9", 961.0597510718158)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i10", 16547.68197801902)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i11", 41861.86377694738)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i12", 1121.3142855340068)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i13", 615.9736188554265)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i14", 690.0862705922896)
    )
    variable_list.add_variable(par_est.VariableAlgebraic("e0_V_Reactor", 0.03341221))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i1_u9", 0.5116581))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i13_u9", 0.5278731))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i14_u9", 1.5814161))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i2_u9", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i3_u9", 0.019119525))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i4_u9", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i5_u9", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i6_u9", 0.024320394))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i7_u9", 0.04427644))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i9_u9", 19.64621))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_c_i10_u9", 0.001505191))
    variable_list.add_variable(par_est.VariableParameter("e0_greek_DeltaG_r1", 6.1222))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_K_eq_r1", 1.1757618))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_r1", 1351500.2))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_r2", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_r3", 1.5e-18))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_r4", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i1", -1351500.2))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i2", 1.5e-18))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i3", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i4", 1351500.2))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i5", 0.0))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i6", -1.5e-18))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i9", 1351500.2))
    variable_list.add_variable(par_est.VariableAlgebraic("e0_r_i13", -1351500.2))
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i6", 1.9335168620499728)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i8", 27.067704397594195)
    )
    variable_list.add_variable(
        par_est.VariableAlgebraic("e0_greek_rho_u9_i7", 27.51435517638899)
    )
    variable_list.add_variable(par_est.VariableConstant("e0_M_i1", 170.29))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i2", 172.3077))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i3", 340.5784))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i4", 226.43))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i5", 227.43))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i9", 18.01528))
    variable_list.add_variable(par_est.VariableParameter("e0_HU_i12", 2.24))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i12", 572.762))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i13", 73.14))
    variable_list.add_variable(par_est.VariableParameter("e0_HU_i14", 9.0))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i14", 170.33))
    variable_list.add_variable(par_est.VariableControl("e0_T", 373.0))
    variable_list.add_variable(par_est.VariableControl("e0_p", 30.0))
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol1_SolP1", -8209.592)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol1_SolP2", 117.0208)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol1_SolP3", -0.20219)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol2_SolP1", -8209.592)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol2_SolP2", 117.0208)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol2_SolP3", -0.20219)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol3_SolP1", -8209.592)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol3_SolP2", 117.0208)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol3_SolP3", -0.20219)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol4_SolP1", -8209.592)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol4_SolP2", 117.0208)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol4_SolP3", -0.20219)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol5_SolP1", -8209.592)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol5_SolP2", 117.0208)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol5_SolP3", -0.20219)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol9_SolP1", 406189.431)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol9_SolP2", -1243.066)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol9_SolP3", 0.95534)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol12_SolP1", 406189.431)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol12_SolP2", -1243.066)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol12_SolP3", 0.95534)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol13_SolP1", 406189.431)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol13_SolP2", -1243.066)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol13_SolP3", 0.95534)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol14_SolP1", -790.7257)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol14_SolP2", 8.3594)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i7_Sol14_SolP3", -0.012265)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol1_SolP1", 196006.2531)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol1_SolP2", -4898.8245)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol1_SolP3", 30.4876)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol2_SolP1", 196006.2531)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol2_SolP2", -4898.8245)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol2_SolP3", 30.4876)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol3_SolP1", 196006.2531)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol3_SolP2", -4898.8245)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol3_SolP3", 30.4876)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol4_SolP1", 89706.8526)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol4_SolP2", -2230.584)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol4_SolP3", 13.9085)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol5_SolP1", 89706.8526)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol5_SolP2", -2230.584)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol5_SolP3", 13.9085)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol9_SolP1", -547307.7898)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol9_SolP2", 3779.8068)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol9_SolP3", -5.7113)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol12_SolP1", -547307.7898)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol12_SolP2", 3779.8068)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol12_SolP3", -5.7113)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol13_SolP1", -1086.9672)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol13_SolP2", 34.3488)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol13_SolP3", -0.21305)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol14_SolP1", 18109.9827)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol14_SolP2", -83.742)
    )
    variable_list.add_variable(
        par_est.VariableConstant("e0_P_GLE_i6_Sol14_SolP3", 0.10281)
    )
    variable_list.add_variable(par_est.VariableConstant("e0_M_i6", 2.01588))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i7", 28.01))
    variable_list.add_variable(par_est.VariableParameter("e0_HU_i10", 0.0156))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i10", 310.19))
    variable_list.add_variable(par_est.VariableParameter("e0_HU_i11", 0.15875))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i11", 784.71))
    variable_list.add_variable(par_est.VariableConstant("e0_R", 8.314))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r1", 12.0197049848239))
    variable_list.add_variable(
        par_est.VariableParameter("e0_k_ref_r1", 7783626.90046052)
    )
    variable_list.add_variable(par_est.VariableParameter("e0_E_r2", 26519.645017333))
    variable_list.add_variable(
        par_est.VariableParameter("e0_k_ref_r2", 1856644.93679767)
    )
    variable_list.add_variable(par_est.VariableParameter("e0_E_r3", 466761.0))
    variable_list.add_variable(par_est.VariableParameter("e0_k_ref_r3", 0.001))
    variable_list.add_variable(par_est.VariableParameter("e0_E_r4", 6127.65655628595))
    variable_list.add_variable(
        par_est.VariableParameter("e0_k_ref_r4", 0.693575064095538)
    )
    variable_list.add_variable(par_est.VariableConstant("e0_a_i14", 1.9385))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i14", 0.58748))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i14", 506.0108))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i14", 0.71269))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i1", 1.505))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i1", 0.50133))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i1", 598.4108))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i1", 0.52735))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i9", 0.02848))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i9", 0.02203))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i9", 447.246))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i9", 0.014091))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i2", 0.56598))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i2", 0.31081))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i2", 669.5672))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i2", 0.31266))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i10", 0.02848))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i10", 0.02203))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i10", 447.246))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i10", 0.014091))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i5", 0.25878))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i5", 0.21792))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i5", 692.8852))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i5", 0.28804))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i13", 1.9404e-10))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i13", 3.4576e-6))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i13", -757347.676))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i13", -108.8259))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i4", 0.25878))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i4", 0.21792))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i4", 692.8852))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i4", 0.28804))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i12", 0.097181))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i12", 0.20184))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i12", 1069.2551))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i12", 0.3072))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i3", 0.13876))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i3", 0.18201))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i3", 816.1091))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i3", 0.28091))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i11", 0.02848))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i11", 0.02203))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i11", 447.246))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i11", 0.014091))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i6", 0.50448))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i6", -0.0027972))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i6", 0.066184))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i6", 3.8561e-6))
    variable_list.add_variable(par_est.VariableConstant("e0_e_i6", -9.151e-5))
    variable_list.add_variable(par_est.VariableConstant("e0_M_i8", 28.0134))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i8", 0.55105))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i8", -0.0030955))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i8", 0.068661))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i8", 4.307e-6))
    variable_list.add_variable(par_est.VariableConstant("e0_e_i8", -9.7341e-5))
    variable_list.add_variable(par_est.VariableConstant("e0_a_i7", -0.5581))
    variable_list.add_variable(par_est.VariableConstant("e0_b_i7", 0.0031211))
    variable_list.add_variable(par_est.VariableConstant("e0_c_i7", 0.063578))
    variable_list.add_variable(par_est.VariableConstant("e0_d_i7", -4.3159e-6))
    variable_list.add_variable(par_est.VariableConstant("e0_e_i7", -8.3167e-5))

    m = par_est.Model(variable_list)

    dydx1 = (
        (
            -(
                m.varlist_all["e0_k_ref_r1"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r1"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    * (
                        (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    - (
                        (
                            (
                                (
                                    m.varlist_all["e0_HU_i4"].casadi_var
                                    / m.varlist_all["e0_M_i4"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                            * (
                                (
                                    m.varlist_all["e0_HU_i9"].casadi_var
                                    / m.varlist_all["e0_M_i9"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                        )
                        / (
                            (
                                ca.exp(
                                    (m.varlist_all["e0_greek_DeltaG_r1"].casadi_var)
                                    / (
                                        m.varlist_all["e0_R"].casadi_var
                                        * m.varlist_all["e0_T"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                )
            )
            - (
                (
                    (
                        m.varlist_all["e0_HU_i10"].casadi_var
                        / m.varlist_all["e0_M_i10"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * m.varlist_all["e0_k_ref_r3"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r3"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i1"].casadi_var
                        / m.varlist_all["e0_M_i1"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        (
                            (
                                m.varlist_all["e0_HU_i1"].casadi_var
                                / m.varlist_all["e0_M_i1"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol1_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i2"].casadi_var
                                / m.varlist_all["e0_M_i2"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol2_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i2"].casadi_var
                            / m.varlist_all["e0_M_i2"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i3"].casadi_var
                                / m.varlist_all["e0_M_i3"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol3_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i3"].casadi_var
                            / m.varlist_all["e0_M_i3"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i4"].casadi_var
                                / m.varlist_all["e0_M_i4"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol4_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i4"].casadi_var
                            / m.varlist_all["e0_M_i4"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i5"].casadi_var
                                / m.varlist_all["e0_M_i5"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol5_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i5"].casadi_var
                            / m.varlist_all["e0_M_i5"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i9"].casadi_var
                                / m.varlist_all["e0_M_i9"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol9_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i9"].casadi_var
                            / m.varlist_all["e0_M_i9"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i12"].casadi_var
                                / m.varlist_all["e0_M_i12"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol12_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i12"].casadi_var
                            / m.varlist_all["e0_M_i12"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i13"].casadi_var
                                / m.varlist_all["e0_M_i13"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol13_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i14"].casadi_var
                                / m.varlist_all["e0_M_i14"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol14_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i14"].casadi_var
                            / m.varlist_all["e0_M_i14"].casadi_var
                        )
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
            - (
                m.varlist_all["e0_k_ref_r4"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r4"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i1"].casadi_var
                        / m.varlist_all["e0_M_i1"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i4"].casadi_var
                        / m.varlist_all["e0_M_i4"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i1"].casadi_var
        * 60
    )

    dydx2 = (
        (
            (
                (
                    (
                        m.varlist_all["e0_HU_i10"].casadi_var
                        / m.varlist_all["e0_M_i10"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * m.varlist_all["e0_k_ref_r3"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r3"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i1"].casadi_var
                        / m.varlist_all["e0_M_i1"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        (
                            (
                                m.varlist_all["e0_HU_i1"].casadi_var
                                / m.varlist_all["e0_M_i1"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol1_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i2"].casadi_var
                                / m.varlist_all["e0_M_i2"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol2_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i2"].casadi_var
                            / m.varlist_all["e0_M_i2"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i3"].casadi_var
                                / m.varlist_all["e0_M_i3"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol3_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i3"].casadi_var
                            / m.varlist_all["e0_M_i3"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i4"].casadi_var
                                / m.varlist_all["e0_M_i4"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol4_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i4"].casadi_var
                            / m.varlist_all["e0_M_i4"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i5"].casadi_var
                                / m.varlist_all["e0_M_i5"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol5_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i5"].casadi_var
                            / m.varlist_all["e0_M_i5"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i9"].casadi_var
                                / m.varlist_all["e0_M_i9"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol9_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i9"].casadi_var
                            / m.varlist_all["e0_M_i9"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i12"].casadi_var
                                / m.varlist_all["e0_M_i12"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol12_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i12"].casadi_var
                            / m.varlist_all["e0_M_i12"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i13"].casadi_var
                                / m.varlist_all["e0_M_i13"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol13_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i14"].casadi_var
                                / m.varlist_all["e0_M_i14"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol14_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i14"].casadi_var
                            / m.varlist_all["e0_M_i14"].casadi_var
                        )
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i2"].casadi_var
        * 60
    )

    dydx3 = (
        (
            (
                m.varlist_all["e0_k_ref_r4"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r4"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i1"].casadi_var
                        / m.varlist_all["e0_M_i1"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i4"].casadi_var
                        / m.varlist_all["e0_M_i4"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i3"].casadi_var
        * 60
    )

    dydx4 = (
        (
            (
                m.varlist_all["e0_k_ref_r1"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r1"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    * (
                        (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    - (
                        (
                            (
                                (
                                    m.varlist_all["e0_HU_i4"].casadi_var
                                    / m.varlist_all["e0_M_i4"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                            * (
                                (
                                    m.varlist_all["e0_HU_i9"].casadi_var
                                    / m.varlist_all["e0_M_i9"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                        )
                        / (
                            (
                                ca.exp(
                                    (m.varlist_all["e0_greek_DeltaG_r1"].casadi_var)
                                    / (
                                        m.varlist_all["e0_R"].casadi_var
                                        * m.varlist_all["e0_T"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                )
            )
            - (
                (
                    (
                        m.varlist_all["e0_HU_i10"].casadi_var
                        / m.varlist_all["e0_M_i10"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * m.varlist_all["e0_k_ref_r2"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r2"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i4"].casadi_var
                        / m.varlist_all["e0_M_i4"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        (
                            (
                                m.varlist_all["e0_HU_i1"].casadi_var
                                / m.varlist_all["e0_M_i1"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol1_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i2"].casadi_var
                                / m.varlist_all["e0_M_i2"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol2_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i2"].casadi_var
                            / m.varlist_all["e0_M_i2"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i3"].casadi_var
                                / m.varlist_all["e0_M_i3"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol3_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i3"].casadi_var
                            / m.varlist_all["e0_M_i3"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i4"].casadi_var
                                / m.varlist_all["e0_M_i4"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol4_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i4"].casadi_var
                            / m.varlist_all["e0_M_i4"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i5"].casadi_var
                                / m.varlist_all["e0_M_i5"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol5_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i5"].casadi_var
                            / m.varlist_all["e0_M_i5"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i9"].casadi_var
                                / m.varlist_all["e0_M_i9"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol9_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i9"].casadi_var
                            / m.varlist_all["e0_M_i9"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i12"].casadi_var
                                / m.varlist_all["e0_M_i12"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol12_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i12"].casadi_var
                            / m.varlist_all["e0_M_i12"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i13"].casadi_var
                                / m.varlist_all["e0_M_i13"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol13_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i14"].casadi_var
                                / m.varlist_all["e0_M_i14"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol14_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i14"].casadi_var
                            / m.varlist_all["e0_M_i14"].casadi_var
                        )
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
            - (
                m.varlist_all["e0_k_ref_r4"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r4"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i1"].casadi_var
                        / m.varlist_all["e0_M_i1"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i4"].casadi_var
                        / m.varlist_all["e0_M_i4"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i4"].casadi_var
        * 60
    )

    dydx5 = (
        (
            (
                (
                    (
                        m.varlist_all["e0_HU_i10"].casadi_var
                        / m.varlist_all["e0_M_i10"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * m.varlist_all["e0_k_ref_r2"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r2"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i4"].casadi_var
                        / m.varlist_all["e0_M_i4"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        (
                            (
                                m.varlist_all["e0_HU_i1"].casadi_var
                                / m.varlist_all["e0_M_i1"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol1_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol1_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i2"].casadi_var
                                / m.varlist_all["e0_M_i2"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol2_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol2_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i2"].casadi_var
                            / m.varlist_all["e0_M_i2"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i3"].casadi_var
                                / m.varlist_all["e0_M_i3"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol3_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol3_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i3"].casadi_var
                            / m.varlist_all["e0_M_i3"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i4"].casadi_var
                                / m.varlist_all["e0_M_i4"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol4_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol4_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i4"].casadi_var
                            / m.varlist_all["e0_M_i4"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i5"].casadi_var
                                / m.varlist_all["e0_M_i5"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol5_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol5_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i5"].casadi_var
                            / m.varlist_all["e0_M_i5"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i9"].casadi_var
                                / m.varlist_all["e0_M_i9"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol9_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP2"].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all["e0_P_GLE_i6_Sol9_SolP3"].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i9"].casadi_var
                            / m.varlist_all["e0_M_i9"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i12"].casadi_var
                                / m.varlist_all["e0_M_i12"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol12_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol12_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i12"].casadi_var
                            / m.varlist_all["e0_M_i12"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i13"].casadi_var
                                / m.varlist_all["e0_M_i13"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol13_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol13_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        + (
                            (
                                m.varlist_all["e0_HU_i14"].casadi_var
                                / m.varlist_all["e0_M_i14"].casadi_var
                            )
                        )
                        / (
                            1.0
                            - (m.varlist_all["e0_p"].casadi_var)
                            / (
                                2.0
                                * (
                                    m.varlist_all["e0_P_GLE_i6_Sol14_SolP1"].casadi_var
                                    + m.varlist_all["e0_T"].casadi_var
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP2"
                                    ].casadi_var
                                    + (m.varlist_all["e0_T"].casadi_var) ** (2.0)
                                    * m.varlist_all[
                                        "e0_P_GLE_i6_Sol14_SolP3"
                                    ].casadi_var
                                )
                            )
                        )
                        - (
                            m.varlist_all["e0_HU_i14"].casadi_var
                            / m.varlist_all["e0_M_i14"].casadi_var
                        )
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i5"].casadi_var
        * 60
    )

    dydx6 = (
        (
            (
                m.varlist_all["e0_k_ref_r1"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r1"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    * (
                        (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    - (
                        (
                            (
                                (
                                    m.varlist_all["e0_HU_i4"].casadi_var
                                    / m.varlist_all["e0_M_i4"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                            * (
                                (
                                    m.varlist_all["e0_HU_i9"].casadi_var
                                    / m.varlist_all["e0_M_i9"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                        )
                        / (
                            (
                                ca.exp(
                                    (m.varlist_all["e0_greek_DeltaG_r1"].casadi_var)
                                    / (
                                        m.varlist_all["e0_R"].casadi_var
                                        * m.varlist_all["e0_T"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i9"].casadi_var
        * 60
    )

    dydx7 = (
        (
            -(
                m.varlist_all["e0_k_ref_r1"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r1"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        (
                            m.varlist_all["e0_HU_i1"].casadi_var
                            / m.varlist_all["e0_M_i1"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    * (
                        (
                            m.varlist_all["e0_HU_i13"].casadi_var
                            / m.varlist_all["e0_M_i13"].casadi_var
                        )
                        / (
                            (
                                (m.varlist_all["e0_HU_i1"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i1"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i1"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i1"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i1"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i2"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i2"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i2"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i2"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i2"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i3"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i3"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i3"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i3"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i3"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i4"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i4"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i4"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i4"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i4"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i5"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i5"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i5"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i5"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i5"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i9"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i9"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i9"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i9"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                                            )
                                        )
                                        * m.varlist_all["e0_M_i9"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i10"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i10"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i10"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i10"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i10"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i10"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i11"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i11"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i11"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i11"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i11"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i11"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i12"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i12"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i12"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i12"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i12"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i12"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i13"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i13"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i13"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i13"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i13"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i13"].casadi_var
                                    )
                                )
                            )
                            + (
                                (m.varlist_all["e0_HU_i14"].casadi_var)
                                / (
                                    (
                                        (m.varlist_all["e0_a_i14"].casadi_var)
                                        / (
                                            (m.varlist_all["e0_b_i14"].casadi_var)
                                            ** (
                                                1.0
                                                + (
                                                    1.0
                                                    - (m.varlist_all["e0_T"].casadi_var)
                                                    / (
                                                        m.varlist_all[
                                                            "e0_c_i14"
                                                        ].casadi_var
                                                    )
                                                )
                                                ** (
                                                    m.varlist_all["e0_d_i14"].casadi_var
                                                )
                                            )
                                        )
                                        * m.varlist_all["e0_M_i14"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                    - (
                        (
                            (
                                (
                                    m.varlist_all["e0_HU_i4"].casadi_var
                                    / m.varlist_all["e0_M_i4"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                            * (
                                (
                                    m.varlist_all["e0_HU_i9"].casadi_var
                                    / m.varlist_all["e0_M_i9"].casadi_var
                                )
                                / (
                                    (
                                        (m.varlist_all["e0_HU_i1"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i1"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i1"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i1"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i1"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i1"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i2"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i2"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i2"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i2"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i2"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i2"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i3"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i3"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i3"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i3"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i3"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i3"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i4"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i4"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i4"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i4"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i4"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i4"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i5"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i5"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i5"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i5"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i5"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i5"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i9"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i9"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i9"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i9"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i9"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i9"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i10"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i10"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i10"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i10"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i10"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i10"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i11"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i11"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i11"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i11"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i11"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i11"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i12"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i12"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i12"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i12"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i12"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i12"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i13"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i13"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i13"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i13"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i13"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i13"].casadi_var
                                            )
                                        )
                                    )
                                    + (
                                        (m.varlist_all["e0_HU_i14"].casadi_var)
                                        / (
                                            (
                                                (m.varlist_all["e0_a_i14"].casadi_var)
                                                / (
                                                    (
                                                        m.varlist_all[
                                                            "e0_b_i14"
                                                        ].casadi_var
                                                    )
                                                    ** (
                                                        1.0
                                                        + (
                                                            1.0
                                                            - (
                                                                m.varlist_all[
                                                                    "e0_T"
                                                                ].casadi_var
                                                            )
                                                            / (
                                                                m.varlist_all[
                                                                    "e0_c_i14"
                                                                ].casadi_var
                                                            )
                                                        )
                                                        ** (
                                                            m.varlist_all[
                                                                "e0_d_i14"
                                                            ].casadi_var
                                                        )
                                                    )
                                                )
                                                * m.varlist_all["e0_M_i14"].casadi_var
                                            )
                                        )
                                    )
                                )
                            )
                        )
                        / (
                            (
                                ca.exp(
                                    (m.varlist_all["e0_greek_DeltaG_r1"].casadi_var)
                                    / (
                                        m.varlist_all["e0_R"].casadi_var
                                        * m.varlist_all["e0_T"].casadi_var
                                    )
                                )
                            )
                        )
                    )
                )
            )
            + (
                m.varlist_all["e0_k_ref_r4"].casadi_var
                * ca.exp(
                    -(m.varlist_all["e0_E_r4"].casadi_var)
                    / (
                        m.varlist_all["e0_R"].casadi_var
                        * m.varlist_all["e0_T"].casadi_var
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i1"].casadi_var
                        / m.varlist_all["e0_M_i1"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
                * (
                    (
                        m.varlist_all["e0_HU_i4"].casadi_var
                        / m.varlist_all["e0_M_i4"].casadi_var
                    )
                    / (
                        (
                            (m.varlist_all["e0_HU_i1"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i1"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i1"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i1"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i1"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i1"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i2"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i2"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i2"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i2"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i2"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i2"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i3"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i3"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i3"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i3"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i3"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i3"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i4"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i4"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i4"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i4"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i4"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i4"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i5"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i5"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i5"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i5"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i5"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i5"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i9"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i9"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i9"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i9"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i9"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i9"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i10"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i10"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i10"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i10"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i10"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i10"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i11"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i11"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i11"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i11"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i11"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i11"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i12"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i12"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i12"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i12"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i12"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i12"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i13"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i13"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i13"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i13"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i13"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i13"].casadi_var
                                )
                            )
                        )
                        + (
                            (m.varlist_all["e0_HU_i14"].casadi_var)
                            / (
                                (
                                    (m.varlist_all["e0_a_i14"].casadi_var)
                                    / (
                                        (m.varlist_all["e0_b_i14"].casadi_var)
                                        ** (
                                            1.0
                                            + (
                                                1.0
                                                - (m.varlist_all["e0_T"].casadi_var)
                                                / (m.varlist_all["e0_c_i14"].casadi_var)
                                            )
                                            ** (m.varlist_all["e0_d_i14"].casadi_var)
                                        )
                                    )
                                    * m.varlist_all["e0_M_i14"].casadi_var
                                )
                            )
                        )
                    )
                )
            )
        )
        * (
            (
                (m.varlist_all["e0_HU_i1"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i1"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i1"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i1"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i1"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i1"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i2"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i2"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i2"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i2"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i2"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i2"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i3"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i3"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i3"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i3"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i3"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i3"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i4"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i4"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i4"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i4"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i4"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i4"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i5"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i5"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i5"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i5"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i5"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i5"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i9"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i9"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i9"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i9"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i9"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i9"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i10"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i10"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i10"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i10"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i10"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i10"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i11"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i11"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i11"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i11"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i11"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i11"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i12"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i12"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i12"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i12"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i12"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i12"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i13"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i13"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i13"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i13"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i13"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i13"].casadi_var
                    )
                )
            )
            + (
                (m.varlist_all["e0_HU_i14"].casadi_var)
                / (
                    (
                        (m.varlist_all["e0_a_i14"].casadi_var)
                        / (
                            (m.varlist_all["e0_b_i14"].casadi_var)
                            ** (
                                1.0
                                + (
                                    1.0
                                    - (m.varlist_all["e0_T"].casadi_var)
                                    / (m.varlist_all["e0_c_i14"].casadi_var)
                                )
                                ** (m.varlist_all["e0_d_i14"].casadi_var)
                            )
                        )
                        * m.varlist_all["e0_M_i14"].casadi_var
                    )
                )
            )
        )
        * m.varlist_all["e0_M_i13"].casadi_var
        * 60
    )
    m.add_equations_differential(
        [dydx1, dydx2, dydx3, dydx4, dydx5, dydx6, dydx7,]
    )

    return variable_list, m


if __name__ == "__main__":

    variable_list, m = initialize_problem()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(0, 4, 10)
    # time_grid = np.insert(time_grid, 0, 0)

    for var in variable_list.values():
        var.fixed = True
        if isinstance(var, par_est.VariableControl) or isinstance(
            var, par_est.ParameterEstimation
        ):
            var.lower_bound = var.value - var.value * 0.05
            var.upper_bound = var.value + var.value * 0.05
            var.guess = var.value
    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    sim_fixed = par_est.Simulator(
        m, time_grid, var_list_fixed, integrator_name="cvodes"
    )

    aa = sim_fixed.generate_exp_data()
    aa.plot_states()
    bb = sim_fixed.simulate(True)
    par_est.tools.plot_array(bb["jac_xf_p"])
