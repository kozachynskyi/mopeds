import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm

import par_est

def plot_res(pe, parameters):
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
    varlist = oed.list_input_varlist[0]
    control_bounds = {}
    for var_name in ["u1", "u2"]:
        lb = varlist[var_name].lower_bound
        ub = varlist[var_name].upper_bound
        num_p = 15
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
    varlist, m_cantois, exp_data = par_est.examples.yeast_growth("cantois")

    time_grid = np.array([0, 5, 10, 15, 20])
    sim_cantois = par_est.Simulator(m_cantois, time_grid, varlist)

    varlist, m_monod, _ = par_est.examples.yeast_growth("monod", False)
    varlist["theta1"].value = 0.531
    varlist["theta2"].value = 7.854
    varlist["theta3"].value = 0.474
    varlist["theta4"].value = 0.019

    # sim_monod = par_est.Simulator(m_monod, time_grid, varlist)
    # sim_monod.generate_exp_data().plot(show=False)
    # sim_cantois.generate_exp_data().plot(show=True)

    pe = par_est.ParameterEstimation(m_monod, exp_data)
    p_preliminary = {
        "theta1": 0.531,
        "theta2": 7.854,
        "theta3": 0.474,
        "theta4": 0.019,
    }

    # pe.solver_settings["ipopt"]["linear_solver"] = "ma57"

    # print(pe.optimize(True, objective_function="wls"))
    p_wls = {'theta1': 0.5490962817206284, 'theta2': 8.359345508729714, 'theta3': 0.47235866666073845, 'theta4': 0.018064211744886684}
    # print(pe.optimize(True, objective_function="ols"))
    # p_ols = {'theta1': 0.37625871088136736, 'theta2': 4.273832042202273, 'theta3': 0.46807664082825035, 'theta4': 0.017818492265246618}

    # print(pe.parameter_analysis(p_preliminary))

    # plot_res(pe, p_wls)
    # plot_res(pe, p_preliminary)


    for var_name, var_value in p_preliminary.items():
        varlist[var_name].value = var_value
        varlist[var_name].fixed = False
    varlist["u1"].fixed = False
    varlist["u2"].fixed = False
    varlist["x1"].fixed = True
    # varlist["u2"].value = 35

    # varlist["u2"].expand_horizon([1], [35])


    controls = {"u1": 0.12, "u2": 35, "x1": 5}
    # oed = par_est.OptimalExperimentalDesign(m_monod, [varlist], time_grid)
    # oed.optimize()
    # a = oed.calculate_objective_and_jacobian(controls)
    # print(a)
    # print(a["jac"].shape)

    oed_settings = par_est.OEDsettings(20, 0.1, 4, 2)

    oed = par_est.OptimalExperimentalDesign(m_monod, [varlist], time_grid, oed_settings)
    # a = oed.optimize(1e-3)
    # print(a)

    plot_fig3(oed, "A")
    # plot_fig3(oed, "D")

    # a = oed.calculate_objective_and_jacobian({"u1": 0.2, "u2": 35})
    # j = a["jac"]
    # 1.75e14
    # print(np.linalg.det(j.T @ j) / 1e14)

    print(oed.optimize(1e-5, objective_function="A_fd"))
