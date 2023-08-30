import copy

import matplotlib.pyplot as plt
import numpy as np

import par_est

def get_varlist_paper(data_set=1, normalized=True):
    varlist, model, _ = par_est.examples.yeast_growth("monod", True, ode=True, normalize=normalized)

    if normalized:
        p_true = {
            "theta1": 1,
            "theta2": 1,
            "theta3": 1,
            "theta4": 1,
        }
    else:
        p_true = {
            "theta1": 0.310,
            "theta2": 0.18,
            "theta3": 0.55,
            "theta4": 0.05,
        }

    for var_name, var_value in p_true.items():
        varlist[var_name].value = var_value
        varlist[var_name].fixed = False

    time_grid = np.linspace(0, 48, 90)

    varlist["x1"].fixed = True
    varlist["u1"].fixed = True
    varlist["u2"].fixed = True
    varlist["x2"].value = 0


    oed_settings = par_est.AdaptiveSampling(True, num_control_switches=4, 20, 0.1, 4, 4)

    if data_set == 1:
        varlist["x1"].value = 2.5
        u1_paper = [
            0.190,
            0.197,
            0.122,
            0.103,
            0.130,
            0.189,
        ]

        u2_paper = [
            34.86,
            20.470,
            8.880,
            8.992,
            6.976,
            25.454,
        ]
        u2time = [
            0,
            7.589,
            15.119,
            23.675,
            29.288,
            39.144,
        ]
    elif data_set == 2:
        varlist["x1"].value = 5.24
        u1_paper = [
            0.129,
            0.082,
            0.075,
            0.098,
            0.073,
            0.057,
        ]

        u2_paper = [
            19.960,
            22.613,
            5.106,
            4.947,
            4.947,
            4.947,
        ]
        u2time = [
            0,
            5.988,
            14.880,
            24.032,
            30.978,
            45.907,
        ]

    elif data_set == 3:
        varlist["x1"].value = 5.5
        u1_paper = [
            0.125,
            0.125,
            0.125,
            0.125,
            0.2,
            0.061,
        ]

        u2_paper = [
            5.000,
            10.084,
            5.000,
            20.093,
            35.026,
            5.053,
        ]
        u2time = [
            0,
            7.385,
            15.093,
            22.478,
            29.927,
            42.105,
        ]

    varlist["u1"].expand_horizon(u2time[1:], u1_paper[1:])
    varlist["u2"].expand_horizon(u2time[1:], u2_paper[1:])
    varlist["u1"].value = u1_paper[0]
    varlist["u2"].value = u2_paper[0]
    # breakpoint()

    oed = par_est.OptimalExperimentalDesign(model, [varlist], time_grid, oed_settings)
    exp_varlist = oed.generate_experimental_data({}, p_true)
    # plot_controls(exp_varlist)
    # exp_varlist.plot()
    # plt.show()
    parameter_accuracy(model, exp_varlist)
    return varlist, time_grid, oed

def plot_controls(varlist):
    fig, ax = plt.subplots(1, 1)
    df = varlist["u1"].dataframe
    ax.step(varlist["u1"].time_relative, df["u1"], where="post")
    ax.set_xlim([0, 48])
    ax.set_ylim([0, 0.28])
    ax2 = ax.twinx()
    ax2.set_ylim([0, 40])
    df = varlist["u2"].dataframe
    ax2.step(varlist["u2"].time_relative, df["u2"], where="post", c="r")

    plt.show()

def plot_from_parameters(oed, controls):
    a = oed.simulate(controls)
    data = a["y"]
    data = np.insert(data, 0, [2.5, 0], 0)

    fig, ax = plt.subplots(1, 1)
    ax.set_ylim([0, 12])
    ax.plot(oed.time_grid_measurements, data[:, 0])
    ax2 = ax.twinx()
    ax2.set_ylim([0, 20])
    ax.set_xlim([0, 48])
    ax2.plot(oed.time_grid_measurements, data[:, 1], c="r")

def get_varlist_solution(normalized=True):
    varlist, model, _ = par_est.examples.yeast_growth("monod", True, ode=True, normalize=normalized)

    if normalized:
        p_true = {
            "theta1": 1,
            "theta2": 1,
            "theta3": 1,
            "theta4": 1,
        }
    else:
        p_true = {
            "theta1": 0.310,
            "theta2": 0.18,
            "theta3": 0.55,
            "theta4": 0.05,
        }

    solution = {
        "x1": 0.4291013149958769,
        "u1_t0": 0.050002814999123055,
        "u1_t1": 0.050000899485218596,
        "u1_t2": 0.07106221453517607,
        "u1_t3": 0.19999892211383605,
        "u1_t4": 0.19999803160181653,
        "u1_t5": 0.05000346420855041,
        "u1_t6": 0.050002519268820884,
        "u1_t7": 0.1999929673670143,
        "u2_t0": 6.977099975970994,
        "u2_t1": 20.234389745507723,
        "u2_t2": 34.99850802462028,
        "u2_t3": 34.999603418956916,
        "u2_t4": 34.99959234643442,
        "u2_t5": 16.743775923515333,
        "u2_t6": 5.001284867366844,
        "u2_t7": 34.999240688357965,
    }

    for var_name, var_value in p_true.items():
        varlist[var_name].value = var_value
        varlist[var_name].fixed = False

    time_grid = np.linspace(0, 48, 9)

    varlist["x1"].fixed = False
    varlist["u1"].fixed = False
    varlist["u2"].fixed = False
    varlist["x2"].value = 0

    varlist["u1"].ignore_plotting = False
    varlist["u2"].ignore_plotting = False

    oed_settings = par_est.OEDsettings(20, 0.1, 4, 1)

    oed = par_est.OptimalExperimentalDesign(model, [varlist], time_grid, oed_settings)
    exp_varlist = oed.generate_experimental_data(solution, p_true)
    exp_varlist.plot()
    parameter_accuracy(model, exp_varlist)

    return oed, exp_varlist

def parameter_accuracy(model, exp_varlist):
    perturbate = True
    if perturbate:
        rng = np.random.default_rng()
        for var_name in ["x1", "x2"]:
            var = exp_varlist[var_name]
            df = var.dataframe
            df[var_name] = rng.normal(df, np.sqrt(var.variance))
            df[var_name][df[var_name] < 0] = 0

    pe = par_est.ParameterEstimation(model, [exp_varlist])
    pe.solver_settings["ipopt"]["max_iter"] = 20
    res = pe.optimize()
    print(res)
    pe.parameter_analysis(res["x_dict"])
    breakpoint()
    # exp_varlist.plot()


if __name__ == "__main__":
    # plot_y()
    get_varlist_paper(3)
    # get_varlist_solution()
    varlist, model, _ = par_est.examples.yeast_growth_ode("monod", True)
    solution = {
        "x1": 0.4291013149958769,
        "u1_t0": 0.050002814999123055,
        "u1_t1": 0.050000899485218596,
        "u1_t2": 0.07106221453517607,
        "u1_t3": 0.19999892211383605,
        "u1_t4": 0.19999803160181653,
        "u1_t5": 0.05000346420855041,
        "u1_t6": 0.050002519268820884,
        "u1_t7": 0.1999929673670143,
        "u2_t0": 6.977099975970994,
        "u2_t1": 20.234389745507723,
        "u2_t2": 34.99850802462028,
        "u2_t3": 34.999603418956916,
        "u2_t4": 34.99959234643442,
        "u2_t5": 16.743775923515333,
        "u2_t6": 5.001284867366844,
        "u2_t7": 34.999240688357965,
    }

    p_true = {
        "theta1": 0.310,
        "theta2": 0.18,
        "theta3": 0.55,
        "theta4": 0.05,
    }

    p_preliminary = {
        "theta1": 0.310 * 1.016,
        "theta2": 0.18 * 0.544,
        "theta3": 0.55 * 1.046,
        "theta4": 0.05 * 1.188,
    }

    for var_name, var_value in p_true.items():
        varlist[var_name].value = var_value
        varlist[var_name].fixed = False

    time_grid = np.linspace(0, 48, 9)

    varlist["x1"].fixed = False
    varlist["u1"].fixed = False
    varlist["u2"].fixed = False
    varlist["x2"].value = 1e-5

    controls = {"u1": 0.12, "u2": 35, "x1": 5}
    oed_settings = par_est.OEDsettings(20, 0.1, 4, 1)

    oed = par_est.OptimalExperimentalDesign(model, [varlist], time_grid, oed_settings)
    print(oed.calculate_objective_and_jacobian(solution))
    # plot_from_parameters(oed, solution)
    plt.show()
    ss = oed.list_simulators[0]

    oed.solver_settings["ipopt"]["max_iter"] = 500
    # oed.solver_settings["ipopt"]["linear_solver"] = "ma57"
    a = oed.optimize(1e-3, objective_function="A_fd")
    # breakpoint()
    # print(a)

    # oed.max_time_experiment = 20
    # oed.min_sampling_delay = 0.5

    # time_grid_dict = {}
    # for i, time in enumerate(time_grid[1:]):
    #     time_grid_dict["time_sp" + str(i)] = time

    # oed.solver_settings["ipopt"]["linear_solver"] = "ma57"
    # a = oed.calculate_objective_and_jacobian(controls | time_grid_dict)
    print(a)
    breakpoint()
    # plot_fig3(oed, "A")
    # plot_fig3(oed, "D")

    # a = oed.calculate_objective_and_jacobian({"u1": 0.2, "u2": 35})
    # j = a["jac"]
    # 1.75e14
    # print(np.linalg.det(j.T @ j) / 1e14)

    print(oed.optimize(1e-5, objective_function="A_fd"))
