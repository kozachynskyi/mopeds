import copy

import matplotlib.pyplot as plt
import numpy as np

import par_est


def plot_res(pe, parameters):
    import matplotlib.pyplot as plt

    data = pe.calculate_objective_and_residual(p_wls, "wls")
    y_all = data["y"]
    y_m = pe.array_data

    y_0 = pe.list_simulators[0]._initial_state
    x = pe.list_simulators[0].time_grid_relative

    y_all = np.insert(y_all, 0, y_0.T, 0)
    y_m = np.insert(y_m, 0, y_0.T, 0)
    fig, ax = plt.subplots(1, 2)
    for y, y_m, axis in zip(y_all.T, y_m.T, ax):
        axis.plot(x, y)
        axis.scatter(x, y_m)
    fig.show()


def plot_fig3(oed, obj_f_name):
    import matplotlib.pyplot as plt
    from matplotlib import cm, ticker

    varlist = oed.list_input_varlist[0]
    control_bounds = {}
    for var_name in ["u1", "u2"]:
        lb = varlist[var_name].lower_bound
        ub = varlist[var_name].upper_bound
        num_p = 7
        control_bounds[var_name] = [lb, ub, num_p]

    grid, meshgrid = par_est.tools.create_grid(control_bounds.values())

    obj_f = []
    for grid_i in grid:
        controls = {"u1": grid_i[0], "u2": grid_i[1]}
        res = oed.calculate_objective_and_jacobian(controls, obj_f_name)
        obj_f.append(res["f"])
    z = np.reshape(obj_f, meshgrid[0].shape)
    plt.contourf(meshgrid[0], meshgrid[1], z, locator=ticker.LogLocator())
    plt.colorbar()
    plt.show()
    print(obj_f)


if __name__ == "__main__":
    varlist, model, _ = par_est.examples.yeast_growth("monod", True)

    p_preliminary = {
        "theta1": 0.310 * 1.016,
        "theta2": 0.18 * 0.544,
        "theta3": 0.55 * 1.046,
        "theta4": 0.05 * 1.188,
    }

    for var_name, var_value in p_preliminary.items():
        varlist[var_name].value = var_value
        varlist[var_name].fixed = False

    time_grid = np.linspace(0, 48, 8)

    varlist["x1"].fixed = False
    varlist["u1"].fixed = False
    varlist["u2"].fixed = False
    varlist["u2"].guess = 20
    varlist["x1"].value = 0#1e-4

    controls = {"u1": 0.12, "u2": 35, "x1": 5}
    oed_settings = par_est.OEDsettings(20, 0.1, 4, 4)

    oed = par_est.OptimalExperimentalDesign(model, [varlist], time_grid, oed_settings)

    u1_paper = [
        0.192,
        0.196,
        0.122,
        0.103,
        0.129,
        0.189,
    ]

    u2_paper = [
        34.694,
        20.561,
        8.878,
        6.888,
        25.51,
    ]
    u2time = [
        -0.031,
        8.031,
        15.066,
        24,
        29.385,
        39.128,
    ]
    np.array(
        [
            0.0,
            1.71428571,
            3.42857143,
            5.14285714,
            6.85714286,
            8.57142857,
            10.28571429,
            12.0,
            13.71428571,
            15.42857143,
            17.14285714,
            18.85714286,
            20.57142857,
            22.28571429,
            24.0,
            25.71428571,
            27.42857143,
            29.14285714,
            30.85714286,
            32.57142857,
            34.28571429,
            36.0,
            37.71428571,
            39.42857143,
            41.14285714,
            42.85714286,
            44.57142857,
            46.28571429,
        ]
    )
    # DM(1.92023e-05)
    paper = {
        "x1": 2.5,
        "u1_t0": u1_paper[0],
        "u1_t1": u1_paper[0],
        "u1_t2": u1_paper[0],
        "u1_t3": u1_paper[0],
        "u1_t4": u1_paper[0],
        "u1_t5": u1_paper[1],
        "u1_t6": u1_paper[1],
        "u1_t7": u1_paper[1],
        "u1_t8": u1_paper[1],
        "u1_t9": u1_paper[2],
        "u1_t10": u1_paper[2],
        "u1_t11": u1_paper[2],
        "u1_t12": u1_paper[2],
        "u1_t13": u1_paper[2],
        "u1_t14": u1_paper[3],
        "u1_t15": u1_paper[3],
        "u1_t16": u1_paper[3],
        "u1_t17": u1_paper[4],
        "u1_t18": u1_paper[4],
        "u1_t19": u1_paper[4],
        "u1_t20": u1_paper[4],
        "u1_t21": u1_paper[4],
        "u1_t22": u1_paper[4],
        "u1_t23": u1_paper[5],
        "u1_t24": u1_paper[5],
        "u1_t25": u1_paper[5],
        "u1_t26": u1_paper[5],
        "u1_t27": u1_paper[5],
        "u2_t0": u2_paper[0],
        "u2_t1": u2_paper[0],
        "u2_t2": u2_paper[0],
        "u2_t3": u2_paper[0],
        "u2_t4": u2_paper[0],
        "u2_t5": u2_paper[1],
        "u2_t6": u2_paper[1],
        "u2_t7": u2_paper[1],
        "u2_t8": u2_paper[1],
        "u2_t9": u2_paper[2],
        "u2_t10": u2_paper[2],
        "u2_t11": u2_paper[2],
        "u2_t12": u2_paper[2],
        "u2_t13": u2_paper[2],
        "u2_t14": u2_paper[2],
        "u2_t15": u2_paper[2],
        "u2_t16": u2_paper[2],
        "u2_t17": u2_paper[3],
        "u2_t18": u2_paper[3],
        "u2_t19": u2_paper[3],
        "u2_t20": u2_paper[3],
        "u2_t21": u2_paper[3],
        "u2_t22": u2_paper[3],
        "u2_t23": u2_paper[4],
        "u2_t24": u2_paper[4],
        "u2_t25": u2_paper[4],
        "u2_t26": u2_paper[4],
        "u2_t27": u2_paper[4],
    }

    all_u1 = []
    all_u2 = []
    for key, value in paper.items():
        if "u1" in key:
            all_u1.append(value)
        if "u2" in key:
            all_u2.append(value)

    solution = {
        "x1": 0.41150921560729115,
        "u1_t0": 0.19980632487862396,
        "u1_t1": 0.15315730223505433,
        "u1_t2": 0.19622560812407971,
        "u1_t3": 0.19789243439345033,
        "u1_t4": 0.19947916738933305,
        "u1_t5": 0.19951190784769676,
        "u1_t6": 0.1998461316260917,
        "u1_t7": 0.19995061395628277,
        "u1_t8": 0.1935809001806323,
        "u1_t9": 0.19536072578495678,
        "u1_t10": 0.1974824163692348,
        "u1_t11": 0.13921551419993458,
        "u1_t12": 0.19176339013652724,
        "u1_t13": 0.1952277185372248,
        "u1_t14": 0.19714409580064332,
        "u1_t15": 0.19016000198925248,
        "u1_t16": 0.19009369280565921,
        "u1_t17": 0.193698125841117,
        "u1_t18": 0.19594009209969443,
        "u1_t19": 0.06336273570602648,
        "u1_t20": 0.053860954281680445,
        "u1_t21": 0.052610443132993866,
        "u1_t22": 0.05201253973674068,
        "u1_t23": 0.05193174334467016,
        "u1_t24": 0.052792071006538155,
        "u1_t25": 0.05274061373258198,
        "u1_t26": 0.052540149973576956,
        "u1_t27": 0.05289674678446349,
        "u2_t0": 34.9694550838278,
        "u2_t1": 33.34183641077189,
        "u2_t2": 5.101539505955396,
        "u2_t3": 5.022611065684428,
        "u2_t4": 5.008039312933507,
        "u2_t5": 5.005790030704473,
        "u2_t6": 22.613127521092963,
        "u2_t7": 34.99250982142333,
        "u2_t8": 34.09793413323992,
        "u2_t9": 34.44497820011084,
        "u2_t10": 34.73481922832279,
        "u2_t11": 34.53806720274013,
        "u2_t12": 33.81839262092707,
        "u2_t13": 34.436600356638955,
        "u2_t14": 34.70526029067036,
        "u2_t15": 34.70343758787842,
        "u2_t16": 34.39329537588233,
        "u2_t17": 34.62763646301945,
        "u2_t18": 34.76932818290804,
        "u2_t19": 34.22416429409089,
        "u2_t20": 32.59345543931457,
        "u2_t21": 29.100249108924835,
        "u2_t22": 15.850306225181232,
        "u2_t23": 7.550460994452381,
        "u2_t24": 16.40867863843042,
        "u2_t25": 7.3078861225972975,
        "u2_t26": 6.531307237624394,
        "u2_t27": 7.295487413621835,
    }
    fig, ax = plt.subplots(1, 1)
    ax.plot(oed.time_grid_control_switch, all_u1)
    ax.set_ylim([0,0.28])
    ax2 = ax.twinx()
    ax2.set_ylim([0,40])
    ax2.plot(oed.time_grid_control_switch, all_u2, c="r")

    a = oed.simulate(solution)
    data = a["y"]
    data = np.insert(data, 0, [2.5 ,0], 0)

    fig, ax = plt.subplots(1,1)
    ax.set_ylim([0,12])
    ax.plot(oed.time_grid_measurements, data[:,0])
    ax2 = ax.twinx()
    ax2.set_ylim([0,20])
    ax2.plot(oed.time_grid_measurements, data[:,1], c="r")
    plt.show()
    # breakpoint()
    a = oed.optimize(1e-3, objective_function="A_fd")
    breakpoint()
    print(a)

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
