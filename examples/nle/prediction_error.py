import mopeds
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt


        # v = pe_grid.calculate_sensitivity_and_fim(
        #     dict_of_params, list(dict_of_params.keys())
        # )

        # jj = jac_grid["f"]
        # hess = pe_artificial.calculate_sensitivity_and_fim(
        #     dict_of_params, list(dict_of_params.keys())
        # )["cov_par"]
        # hess_e = hess.copy()
        # hess_e[0][1] = 0
        # hess_e[1][0] = 0

        # len_exp = len(experimental_data)
        # len_param = len(dict_of_params)
        # fisher95 = scipy.stats.f(len_param, self.dof).ppf(0.95)

        # inference_results = {}
        # for control in dict_of_controls:
        #     inference_results[control] = np.array(sim_data[control])

        # for response in dict_of_responses:
        #     inference_results[response] = {}
        #     s = np.sqrt(OLS[response] / self.dof)
        #     R = np.linalg.qr(jac[response], mode="reduced")[1]
        #     vv =  np.linalg.norm(jac_grid[response] @ np.linalg.inv(R), axis=1)
        #     my = np.sqrt(np.diag(jj @ hess @ jj.T))
def get_model():
    vl_original, model = mopeds.examples.polynomial_1d()
    vl_original["y"].variance = 0.5**2
    return vl_original, model

def unfix_parameters(list_vl):
    for vl in list_vl:
        vl["a"].fixed = False
        vl["b"].fixed = False
    return list_vl

def parameter_covariance_mc():
    vl_original, model = get_model()

    meas_grid = {"u": [0, 1, 5]}
    true_data, true_params = mopeds.tools.generate_artificial_data_from_grid_nle(model, vl_original, meas_grid, perturbate=False)

    true_data = unfix_parameters(true_data)

    pe = mopeds.ParameterEstimationNLE(model, true_data)
    original_data = copy.deepcopy(pe.array_data)

    rng = np.random.default_rng()

    list_parameters = []
    for i in range(300):
        pe.array_data = rng.normal(original_data, (1 / pe.array_inverted_scaled_std))
        res = pe.optimize(None, "wls", True, reuse_solver=True)
        list_parameters.append(list(res["x_dict"].values()))

    df_params = pd.DataFrame(list_parameters, columns=pe.varlist_decision.keys())
    cov_linearized = pe.calculate_sensitivity_and_fim(true_params)["cov_par"]
    cov_mc = df_params.cov()

    return df_params, cov_linearized, pe

def model_prediction_error_mc():
    vl_original, model = get_model()
    df_params, cov_linearized, pe_covariance = parameter_covariance_mc()

    prediction_grid = {"u": [0, 1, 20]}

    prediction_data, true_params = mopeds.tools.generate_artificial_data_from_grid_nle(model, vl_original, prediction_grid, perturbate=False)

    prediction_data = unfix_parameters(prediction_data)
    pe = mopeds.ParameterEstimationNLE(model, prediction_data)

    list_predictions = []
    for index, row in df_params.iterrows():
        prediction_i = pe.calculate_objective_and_residual(row.to_dict())["y"].flatten()
        list_predictions.append(prediction_i)

    df_predictions = pd.DataFrame(list_predictions)

    jac_prediction = pe.calculate_sensitivity_and_fim(true_params)["jac_sorted"]["y"]

    prediction_line = pe.calculate_objective_and_residual(true_params)["y"].flatten()
    cov_linearized = pe_covariance.calculate_sensitivity_and_fim(row.to_dict())["cov_par"]
    prediction_linearized = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T))


    std_mc = df_predictions.std()


    plt.plot(prediction_line, c="b")
    plt.plot(prediction_line + prediction_linearized, label="lin", c="r")
    plt.plot(prediction_line - prediction_linearized, label="lin", c="r")
    plt.plot(prediction_line + std_mc, label="mc", c="g")
    plt.plot(prediction_line - std_mc, label="mc", c="g")
    plt.legend()
    plt.show()

    # df
    # df.plot()
    # print(cov_linearized)
    # print(cov_mc)
    breakpoint()




if __name__ == "__main__":
    # parameter_covariance_mc()
    model_prediction_error_mc()
