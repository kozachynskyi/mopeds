import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from tools import unfix_parameters
import casadi as ca

import par_est
from quaglio import CriteriaD_asprey2002

def asprey2002(plot=False):
    varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True)

    par_initial = {
        "theta1": 0.5,
        "theta2": 0.5,
        "theta3": 0.5,
        "theta4": 0.5,
    }

    # varlist["x1"].variance = 0.2**2
    # varlist["x2"].variance = 0.2**2

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

    mode = "initial"

    expected = 1

    if mode == "initial":
        time_grid = [0,2,4,6,8,10,12,14,16,18,20]
        varlist["x1"].value = 5.5

        varlist["u1"].value = 0.12

        varlist["u2"].value = 15
        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_asprey2002)
        expected = 2.4e8

    elif mode == "determinant":
        time_grid = [0, 21.2, 22.2, 23.2, 24.2, 25.2, 26.2, 27.2, 28.2, 29.2, 30.2]
        varlist["x1"].value = 8.53

        varlist["u1"].value = 0.2
        varlist["u1"].expand_horizon([5.4], [0.05])

        varlist["u2"].value = 35
        varlist["u2"].expand_horizon([20, 25.2], [22.8, 15.])

        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_asprey2002)
    else:
        raise NotImplementedError
    print(mode)
    print(obj["f"])
    print(f"expected: {expected}")
    print(f"difference obj/expected: {obj['f'] / expected}")

    if plot:
        time_grid = np.linspace(0, 32, 1000)
        sim = par_est.Simulator(m_monod, time_grid, varlist)
        sim.generate_exp_data().plot()

def yeast_oed(mode="initial", estimate=False, plot=False):
    if mode in ["initial", "x_time_zero", "adaptive"]:
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
            res_adaptive = {'u1_t0': 0.0684054281517351, 'u2_t0': 35.0}
        else:
            res = oed.optimize(1e-3)["x_dict"]

        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        sim_res_a = oed.generate_experimental_data(res_adaptive, par_initial)
        breakpoint()
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
        print(obj)
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

        oed_setttings = par_est.AdaptiveSampling(num_control_switches=0, num_sampling_times=len(time_grid), max_time_experiment=time_grid[-1], min_sampling_delay=0.2)

        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid, oed_setttings)
        oed.solver_settings["ipopt"]["linear_solver"] = "ma57"

        if estimate is False:
            res = {'x1': 10.0, 'u1_t0': 0.0684054281517351, 'u2_t0': 35.0, 'time_sp0': 1.13840019411628, 'time_sp1': 1.3384002867398925, 'time_sp2': 1.538400389658182, 'time_sp3': 7.142819465906003, 'time_sp4': 7.342821464722896, 'time_sp5': 7.5428226858854455, 'time_sp6': 7.742823727509143, 'time_sp7': 7.942824797470047, 'time_sp8': 8.142826119543038, 'time_sp9': 8.342828392838715, 'time_sp10': 19.999999629541946}
        else:
            res = oed.optimize(1e-3)["x_dict"]
        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        breakpoint()

        print(res)
        if plot:
            sim_res.plot()
        obj = oed.calculate_objective_and_jacobian(res, "A")
        print(obj)
        # expected = 2.4e8

    elif mode == "optimal":
        varlist["x1"].value = 5.5
        varlist["x1"].fixed = False

        varlist["u1"].guess = 0.2
        varlist["u2"].guess = 5

        varlist["u1"].fixed = False
        varlist["u2"].fixed = False

        varlist_oed = unfix_parameters(varlist)

        extended_time_grid = np.linspace(0, time_grid[-1], 41)
        oed_setttings = par_est.OptimalSampling(end_time_fixed=True, num_sampling_times=len(time_grid))

        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], extended_time_grid, oed_setttings)
        oed.solver_settings["ipopt"]["linear_solver"] = "ma57"

        if estimate is False:
            res = {'x1': 10.0, 'u1': 0.09769358284691101, 'u2': 35.0, 'weight_0': 0.9999999801264828, 'weight_1': 0.9999999939664228, 'weight_2': 0.9999999933918129, 'weight_3': 0.9999999876311378, 'weight_4': 5.795139551782676e-08, 'weight_5': 5.795139294377323e-08, 'weight_6': 5.7951391721197166e-08, 'weight_7': 5.795139141155059e-08, 'weight_8': 0.9999999443618343, 'weight_9': 0.9999999515472652, 'weight_10': 0.9999999585348427, 'weight_11': 5.795139339281108e-08, 'weight_12': 5.795139397256644e-08, 'weight_13': 5.795139444667638e-08, 'weight_14': 5.795139478626125e-08, 'weight_15': 0.9999999739487303, 'weight_16': 0.9999999741735978, 'weight_17': 0.9999999737408243, 'weight_18': 5.795139472822307e-08, 'weight_19': 5.7951394421161226e-08, 'weight_20': 5.795139403671644e-08, 'weight_21': 5.7951393595871046e-08, 'weight_22': 5.795139311815646e-08, 'weight_23': 5.795139262090633e-08, 'weight_24': 5.795139211895278e-08, 'weight_25': 5.795139162449898e-08, 'weight_26': 5.7951391147247304e-08, 'weight_27': 5.7951390694520085e-08, 'weight_28': 5.7951390271529884e-08, 'weight_29': 5.795138988165618e-08, 'weight_30': 5.795138952672992e-08, 'weight_31': 5.795138920729947e-08, 'weight_32': 5.7951388922891324e-08, 'weight_33': 5.7951388672236e-08, 'weight_34': 5.7951388453468e-08, 'weight_35': 5.7951388264299054e-08, 'weight_36': 5.795138810216423e-08, 'weight_37': 5.7951387964342347e-08, 'weight_38': 5.795138784798957e-08, 'weight_39': 0.9999986300353065}
        else:
            res = oed.optimize(1e-3)["x_dict"]
        print(res)
        sim_res = oed.generate_experimental_data(res, par_initial)
        if plot:
            sim_res.plot()
        obj = oed.calculate_objective_and_jacobian(res, "A")
        print(obj)
        # expected = 2.4e8
    breakpoint()

if __name__ == "__main__":
    # asprey2002()

    mode = "initial"
    # mode = "x_time_zero"
    # mode = "adaptive"
    # mode = "optimal"

    yeast_oed(mode, estimate=False, plot=True)
