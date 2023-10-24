import numpy as np

import par_est
import copy
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

def parameter_analysis_monte_carlo(pe, parameters):
    original_solver_settings = copy.deepcopy(pe.solver_settings)
    original_data = copy.deepcopy(pe.array_data)

    original_guess = copy.deepcopy(pe.guess)
    pe.guess = pe.variables_dict_to_list(parameters)

    pe.solver_settings = {"ipopt": {"print_level": 0}}
    pe.solver_settings["print_time"] = False

    all_res = []
    rng = np.random.default_rng()
    num_replications = 40
    for i in tqdm(range(num_replications)):
        artificial_data = pe.calculate_objective_and_residual(parameters)["y"]
        pe.array_data = rng.normal(artificial_data, 1 / pe.array_inverted_std)

        res = pe.optimize()
        all_res.append(res["x"].toarray())

    df = pd.DataFrame(np.squeeze(all_res), columns=pe.varlist_decision.keys())

    parameter_residual = df - df.mean()
    cov_par = (1 / (num_replications - 1)) * parameter_residual.T @ parameter_residual
    cov_lin = pe.calculate_sensitivity_and_fim(parameters)["cov_par"]
    print(cov_par)
    print(cov_lin)
    print(cov_par / cov_lin)
    pe.array_data = original_data
    pe.guess = original_guess
    pe.solver_settings = original_solver_settings
    breakpoint()


if __name__ == "__main__":
    VAR_LIST, MODEL, EXP_DATA = par_est.examples.bod_model()
    pe = par_est.ParameterEstimationNLE(MODEL, EXP_DATA)

    # Example BOD 5 Bates Page 54
    res = pe.optimize()["x_dict"]
    print("expected par 19.143, 0.5311")
    print("Estimated par:\n", res)

    print("S2 expected 6.498")
    print("S^2: ", pe.calculate_objective_and_residual(res, "ols")["f"] / 4)
    a = pe.calculate_sensitivity_and_fim(res)["cov_par"]
    parameter_analysis_monte_carlo(pe, res)
    breakpoint()

    # Plot of convidence 95%-region page 55
    res = pe.parameter_analysis(res, plot=True)
    plt.show()
