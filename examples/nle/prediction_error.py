import mopeds
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt

MODEL_NAME = 2
NUM = 100


def model_selector():
    if MODEL_NAME == 1:
        return get_model_poly1()
    elif MODEL_NAME == 2:
        return get_model_linear_example()
    elif MODEL_NAME == 3:
        return get_model_bod()
    elif MODEL_NAME == 4:
        return get_model_vle_wilson()
    elif MODEL_NAME == 5:
        return get_model_dynamic_example()
    elif MODEL_NAME == 6:
        return get_model_poly1_multivariate()
        # return get()


def scale_df_all(pe, df):
    varlist_alg = pe.model.varlist_algebraic(pe.list_input_varlist[0])
    for i, row in df.iterrows():
        df.iloc[i] = varlist_alg.scale_to_original(row)
    return df


def get_model_vle_wilson():
    vl_original, model = mopeds.examples.vle_wilson()
    vl_original["e0_T"].variance = (0.1 / 2) ** 2
    vl_original["e0_y_c1"].variance = (0.01 / 2) ** 2

    prediction_grid = {"e0_x_c1": [0.01, 0.99, 20], "e0_P": [1, 2, 1]}
    meas_grid = {"e0_x_c1": [0.01, 0.99, 6], "e0_P": [1, 2, 1]}

    unfix_parameters = []
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            unfix_parameters.append(var.name)

    # meas_variables = ["e0_T", "e0_y_c1"]
    meas_variables = ["e0_y_c1"]

    return (
        vl_original,
        model,
        prediction_grid,
        meas_grid,
        unfix_parameters,
        meas_variables,
    )


def get_model_bod():
    vl_original, model, _ = mopeds.examples.bod_model()
    vl_original["f"].variance = 1.5**2
    prediction_grid = {"x": [0, 10, 10]}
    meas_grid = {"x": [0, 10, 6]}

    unfix_parameters = []
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            unfix_parameters.append(var.name)

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, None


def get_model_poly1():
    vl_original, model = mopeds.examples.polynomial_1d()
    vl_original["y"].variance = 0.1**2
    prediction_grid = {"u": [0, 1, 10]}
    meas_grid = {"u": [0, 1, 5]}

    unfix_parameters = []
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            unfix_parameters.append(var.name)
    # unfix_parameters = ["a"]
    # vl_original["b"].value = 2.01

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, None


def get_model_poly1_multivariate():
    vl_original, model = mopeds.examples.polynomial_1d_multivariate()
    vl_original["y1"].variance = 0.1**2
    vl_original["y2"].variance = 0.1**2
    prediction_grid = {"u1": [0, 1, 5], "u2": [0, 1, 5]}
    # meas_grid = {"u": [0, 1, 5]}
    meas_grid = {"u1": [0, 10, 10], "u2": [0, 10, 10]}

    unfix_parameters = []
    meas_name = ["y1", "y2"]
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            if "y2" not in meas_name:
                if "2" in var.name:
                    continue
            unfix_parameters.append(var.name)

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_name


def get_model_linear_example():
    vl_original, model = mopeds.examples.linear_example()
    vl_original["y"].variance = 0.5**2
    prediction_grid = {"u": [0, 1, 10], "v": [3, 4, 1]}
    meas_grid = {"u": [0, 1, 5], "v": [3, 4, 5]}

    unfix_parameters = []
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            unfix_parameters.append(var.name)

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, ["y", "z"]


def get_model_dynamic_example():
    vl_original, model, _ = mopeds.examples.yeast_growth()

    prediction_grid_controls = mopeds.tools.controls_grid_from_dict(
        {"u1": [0.05, 0.2, 10], "u2": [5, 35, 1]}
    )
    list_scenarios_prediction = []
    for controls in prediction_grid_controls:
        list_scenarios_prediction.append(
            {"controls": controls, "time_grid": np.linspace(0, 20, 6), "initials": {}}
        )

    measurement_grid_controls = mopeds.tools.controls_grid_from_dict(
        {"u1": [0.05, 0.2, 3], "u2": [5, 35, 3]}
    )
    list_scenarions_measurement = []
    for inits in [{"x1": 5}, {"x1": 7}]:
        for controls in measurement_grid_controls:
            list_scenarions_measurement.append(
                {
                    "controls": controls,
                    "time_grid": np.linspace(0, 20, 6),
                    "initials": inits,
                }
            )

    unfix_parameters = ["theta1", "theta2", "theta3", "theta4"]

    return (
        vl_original,
        model,
        list_scenarios_prediction,
        list_scenarions_measurement,
        unfix_parameters,
        ["x1", "x2"],
    )


def parameter_covariance_mopeds():
    vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_variables = (
        model_selector()
    )

    # true_parameters = {"a": 20}
    # unfix_parameters.pop(0)
    true_parameters = {}

    if MODEL_NAME == 5:
        pe = mopeds.ParameterEstimation
    else:
        pe = mopeds.ParameterEstimationNLE

    analyzer = mopeds.tools.ErrorAnalyzer(
        vl_original,
        model,
        prediction_grid,
        meas_grid,
        unfix_parameters,
        meas_variables,
        true_parameters=true_parameters,
        pe_class=pe,
    )
    if MODEL_NAME == 5:
        pe = analyzer.pe_main
        pe.solver_settings["ipopt"]["linear_solver"] = "ma57"
        pe.solver_settings["ipopt"]["tol"] = 1e-2
        pe.solver_settings["ipopt"]["max_iter"] = 30
        pe.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"
    analyzer.parameter_covariance_mc(plot=False, num_samples=NUM)
    analyzer.check_linearization_df_params()
    # analyzer.plot_estimation_accuracy()
    # fig, axes = analyzer.plot_parameter_covariance(normalize_parameters=False)
    # analyzer.plot_parameter_covariance_ellipse(normalize_parameters=False)
    analyzer.model_prediction_error_mc(plot=False)
    # print(analyzer.df_s2.mean())
    # print(analyzer.df_s2_true.mean())

    print(res := analyzer.analyze_model_prediction())
    arr = np.array(analyzer.list_prediction_std)
    mean = arr.reshape((arr.shape[0] * arr.shape[1], arr.shape[2])).mean(axis=0)
    mean = arr.reshape((arr.shape[0] * arr.shape[1], arr.shape[2])).std(axis=0)
    print(mean)

    v = analyzer
    breakpoint()

    if MODEL_NAME != 5:
        pass
        # analyzer.model_prediction_error_mc(plot=False)
        # analyzer.plot_model_prediction_MC(without_outliers=True)
        # v = analyzer.analyze_model_prediction()

    plt.show()


if __name__ == "__main__":
    parameter_covariance_mopeds()
