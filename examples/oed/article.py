import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from tools import unfix_parameters
import casadi as ca

import par_est
from quaglio import CriteriaD_asprey2002

def yeast_oed(mode="initial", estimate=False, plot=False):
    if mode in ["initial", "x_time_zero", "adaptive", "optimal"]:
        varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True)
    else:
        varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=False)
    time_grid = [0,2,4,6,8,10,12,14,16,18,20]

    par_initial = {
        "theta1": 0.5,
        "theta2": 0.5,
        "theta3": 0.5,
        "theta4": 0.5,
    }

    varlist["x1"].variance = 1
    varlist["x2"].variance = 1

    varlist["x1"].lower_bound = 1.0
    varlist["x1"].upper_bound = 10.0

    varlist["x2"].value = 0.1
    varlist["x2"].fixed = True

    varlist["u1"].lower_bound = 0.05
    varlist["u1"].upper_bound = 0.2
    varlist["u1"].ignore_plotting = False

    varlist["u2"].lower_bound = 5
    varlist["u2"].upper_bound = 35
    varlist["u2"].ignore_plotting = False

    for par_name, par_value in par_initial.items():
        varlist[par_name].value = par_value

    time_grid = np.linspace(0, 20 ,11)

    if mode == "initial":

        varlist["x1"].value = 5.5
        varlist["x2"].value = 0.01

        varlist["u1"].guess = 0.2
        varlist["u2"].guess = 5

        varlist["u1"].fixed = False
        varlist["u2"].fixed = False

        varlist_oed = unfix_parameters(varlist)

        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        oed.solver_settings["ipopt"]["linear_solver"] = "ma57"

        if estimate is False:
            res = {'u1_t0': 0.06546297244296971, 'u2_t0': 35.0}
        else:
            res = oed.optimize(1e-3)["x_dict"]

        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        if plot:
            sim_res.plot()
        obj = oed.calculate_objective_and_jacobian(res, "A")
        print(obj)
        # expected = 2.4e8

    elif mode == "x_time_zero":

        varlist["x1"].value = 5.5
        varlist["x1"].fixed = False

        varlist["u1"].guess = 0.2
        varlist["u2"].guess = 5

        varlist["u1"].fixed = False
        varlist["u2"].fixed = False

        varlist_oed = unfix_parameters(varlist)

        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        oed.solver_settings["ipopt"]["linear_solver"] = "ma57"

        if estimate is False:
            res = {'x1': 10.0, 'u1_t0': 0.06602044304102403, 'u2_t0': 35.0}
        else:
            res = oed.optimize(1e-3)["x_dict"]
        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        if plot:
            sim_res.plot()
        obj = oed.calculate_objective_and_jacobian(res, "A")
        print(obj["f"])
        # expected = 2.4e8

    elif mode == "adaptive":
        varlist["x1"].value = 5.5
        varlist["x1"].fixed = False

        varlist["u1"].guess = 0.125
        varlist["u2"].guess = 5
        if isinstance(varlist["u2"], par_est.VariableControlPiecewiseConstant):
            varlist["u1"].variable_list["u1_t0"].guess = 0.125
            varlist["u2"].variable_list["u2_t0"].guess = 5

        varlist["u1"].fixed = False
        varlist["u2"].fixed = False
        # a =varlist["u2"]
        # a.expand_horizon([ca.MX.sym("a")], [1])
        # breakpoint()

        varlist_oed = unfix_parameters(varlist)

        oed_setttings = par_est.AdaptiveSampling(num_control_switches=0, num_sampling_times=len(time_grid)-1, max_time_experiment=time_grid[-1], min_sampling_delay=0.2)

        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid, oed_setttings)
        oed.solver_settings["ipopt"]["linear_solver"] = "ma57"

        if estimate is False:
            res = {'x1': 10.0, 'u1_t0': 0.07015821120073494, 'u2_t0': 35.0, 'time_sp0': 1.1097082753731937, 'time_sp1': 1.30970836060938, 'time_sp2': 1.509708455666444, 'time_sp3': 7.113051205708761, 'time_sp4': 7.313052947783498, 'time_sp5': 7.51305406153896, 'time_sp6': 7.713055077904099, 'time_sp7': 7.913056254971746, 'time_sp8': 8.113058197838106, 'time_sp9': 19.999999546450784}
        else:
            res = oed.optimize(1e-3)["x_dict"]
        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        obj = oed.calculate_objective_and_jacobian(res, "A")
        print(obj["f"])
        if plot:
            sim_res.plot(marker="o")

    elif mode == "optimal":
        varlist["x1"].value = 5.5
        varlist["x1"].fixed = False

        varlist["u1"].guess = 0.2
        varlist["u2"].guess = 5

        varlist["u1"].fixed = False
        varlist["u2"].fixed = False

        varlist_oed = unfix_parameters(varlist)

        # adaptive
        extended_time_grid = np.linspace(0, time_grid[-1], 11)
        extended_time_grid = list(extended_time_grid) + [0.0, 1.109708275, 1.309708361, 1.509708456, 7.113051206, 7.313052948, 7.513054062, 7.713055078, 7.913056255, 8.113058198, 19.999999546]
        extended_time_grid = np.unique(extended_time_grid)

        # ideal
        extended_time_grid = np.linspace(0, time_grid[-1], 101)

        oed_setttings = par_est.OptimalSampling(end_time_fixed=True, num_sampling_times=len(time_grid),num_control_switches=0)

        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], extended_time_grid, oed_setttings)
        oed.solver_settings["ipopt"]["linear_solver"] = "ma57"

        if estimate is False:
            res = {'x1': 10.0, 'u1_t0': 0.07480154927217653, 'u2_t0': 35.0, 'weight_0': 4.767843917768539e-08, 'weight_1': 4.767844231443849e-08, 'weight_2': 4.767844547828658e-08, 'weight_3': 4.767844802146793e-08, 'weight_4': 4.7678449654032215e-08, 'weight_5': 0.9999999948593128, 'weight_6': 0.9999999947510796, 'weight_7': 0.9999999937192698, 'weight_8': 4.7678448591972244e-08, 'weight_9': 4.76784473568733e-08, 'weight_10': 4.7678446071095666e-08, 'weight_11': 4.767844484123872e-08, 'weight_12': 4.7678443737263066e-08, 'weight_13': 4.767844279827124e-08, 'weight_14': 4.76784420393945e-08, 'weight_15': 4.7678441458542514e-08, 'weight_16': 4.7678441042357474e-08, 'weight_17': 4.7678440771056567e-08, 'weight_18': 4.767844062208448e-08, 'weight_19': 4.767844057266792e-08, 'weight_20': 4.7678440601443255e-08, 'weight_21': 4.767844068936163e-08, 'weight_22': 4.767844082006883e-08, 'weight_23': 4.7678440979943454e-08, 'weight_24': 4.7678441157910205e-08, 'weight_25': 4.767844134515008e-08, 'weight_26': 4.7678441534747414e-08, 'weight_27': 4.7678441721442555e-08, 'weight_28': 4.767844190124956e-08, 'weight_29': 4.7678442071226875e-08, 'weight_30': 0.9999999422464696, 'weight_31': 0.9999999461887541, 'weight_32': 0.999999949366473, 'weight_33': 0.9999999519273544, 'weight_34': 0.9999999539824775, 'weight_35': 0.9999999556160946, 'weight_36': 4.767844287485024e-08, 'weight_37': 4.767844293034348e-08, 'weight_38': 4.767844297169472e-08, 'weight_39': 4.767844299946093e-08, 'weight_40': 4.767844301426934e-08, 'weight_41': 4.767844301680029e-08, 'weight_42': 4.767844300777239e-08, 'weight_43': 4.7678442987929646e-08, 'weight_44': 4.767844295803393e-08, 'weight_45': 4.7678442918854784e-08, 'weight_46': 4.767844287116134e-08, 'weight_47': 4.7678442815716485e-08, 'weight_48': 4.767844275327133e-08, 'weight_49': 4.767844268456002e-08, 'weight_50': 4.7678442610296874e-08, 'weight_51': 4.767844253117353e-08, 'weight_52': 4.7678442447856006e-08, 'weight_53': 4.7678442360982235e-08, 'weight_54': 4.767844227116085e-08, 'weight_55': 4.767844217897012e-08, 'weight_56': 4.7678442084957187e-08, 'weight_57': 4.7678441989637684e-08, 'weight_58': 4.76784418934957e-08, 'weight_59': 4.7678441796983746e-08, 'weight_60': 4.76784417005232e-08, 'weight_61': 4.7678441604505395e-08, 'weight_62': 4.767844150928901e-08, 'weight_63': 4.7678441415209685e-08, 'weight_64': 4.767844132257065e-08, 'weight_65': 4.767844123164927e-08, 'weight_66': 4.7678441142697315e-08, 'weight_67': 4.767844105594177e-08, 'weight_68': 4.7678440971585985e-08, 'weight_69': 4.767844088981083e-08, 'weight_70': 4.767844081077581e-08, 'weight_71': 4.767844073462025e-08, 'weight_72': 4.7678440661464486e-08, 'weight_73': 4.7678440591410785e-08, 'weight_74': 4.767844052454469e-08, 'weight_75': 4.767844046093577e-08, 'weight_76': 4.767844040063897e-08, 'weight_77': 4.767844034369543e-08, 'weight_78': 4.767844029013351e-08, 'weight_79': 4.7678440239969595e-08, 'weight_80': 4.7678440193209115e-08, 'weight_81': 4.767844014984721e-08, 'weight_82': 4.7678440109869676e-08, 'weight_83': 4.7678440073253626e-08, 'weight_84': 4.767844003996822e-08, 'weight_85': 4.7678445404330534e-08, 'weight_86': 0.9999983388903676, 'weight_87': 0.9999977873925899, 'weight_88': 4.7678439939273276e-08, 'weight_89': 4.767843992194607e-08, 'weight_90': 4.7678439907633446e-08, 'weight_91': 4.767843989626591e-08, 'weight_92': 4.767843988777055e-08, 'weight_93': 4.7678439882071275e-08, 'weight_94': 4.76784398790894e-08, 'weight_95': 4.767843987874383e-08, 'weight_96': 4.7678439880951504e-08, 'weight_97': 4.76784398856277e-08, 'weight_98': 4.7678439892686286e-08, 'weight_99': 4.767843990204004e-08}
            # adaptive
            # res = {'x1': 10.0, 'u1_t0': 0.07327634813201177, 'u2_t0': 35.0, 'weight_0': 0.9999999805088146, 'weight_1': 0.9999999808576362, 'weight_2': 0.9999999789348611, 'weight_3': 0.9999999639940393, 'weight_4': 5.006532840177589e-08, 'weight_5': 5.006535575797373e-08, 'weight_6': 0.9999999638180009, 'weight_7': 0.9999999638305416, 'weight_8': 0.999999963648716, 'weight_9': 0.9999999632859315, 'weight_10': 0.9999999627515819, 'weight_11': 0.9999999624672523, 'weight_12': 5.006535914318206e-08, 'weight_13': 5.006534904238574e-08, 'weight_14': 5.006533647093379e-08, 'weight_15': 5.006532788835129e-08, 'weight_16': 5.006532550644309e-08, 'weight_17': 5.006532927314581e-08, 'weight_18': 5.006554233716967e-08, 'weight_19': 0.999999915249009}
        else:
            res = oed.optimize(1e-3)["x_dict"]
        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        obj = oed.calculate_objective_and_jacobian(res, "A")
        print(obj["f"])
        if plot:
            sim_res.plot(marker="o")
    breakpoint()

if __name__ == "__main__":
    # asprey2002()

    mode = "initial"
    # mode = "x_time_zero"
    mode = "adaptive"
    # mode = "optimal"

    yeast_oed(mode, estimate=False, plot=True)
