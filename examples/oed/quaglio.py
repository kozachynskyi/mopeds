import numpy as np

import par_est

if __name__ == "__main__":
    varlist, m_cantois, exp_data = par_est.examples.yeast_growth("cantois")

    time_grid = np.array([0, 5, 10, 15, 20])
    sim_cantois = par_est.Simulator(m_cantois, time_grid, varlist)

    varlist, m_monod, _ = par_est.examples.yeast_growth("monod")

    sim_monod = par_est.Simulator(m_monod, time_grid, varlist)
    # sim_monod.generate_exp_data().plot(show=False)
    # sim_cantois.generate_exp_data().plot(show=True)

    pe = par_est.ParameterEstimation(m_monod, exp_data)
    p_preliminary = {
        "theta1": 0.531,
        "theta2": 7.854,
        "theta3": 0.474,
        "theta4": 0.019,
    }
    p2 = {
        "theta1": 1.0228490716969507,
        "theta2": 19.999997933391416,
        "theta3": 0.45371925110918304,
        "theta4": 0.010955977763805047,
    }
    print(pe.calculate_objective_and_residual(p_preliminary, "wls"))
    print(pe.calculate_objective_and_residual(p2, "wls"))

    pe.solver_settings["ipopt"]["linear_solver"] = "ma57"
    print(pe.optimize(True, objective_function="ols"))
    # print(pe.optimize(False, objective_function="wls"))
