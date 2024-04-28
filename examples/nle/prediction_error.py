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

    def unfix_parameters(list_vl):
        for vl in list_vl:
            vl["e0_greek_lambdaA_c2_j1"].fixed = False
            vl["e0_greek_lambdaA_c1_j2"].fixed = False
        return list_vl

    # meas_variables = ["e0_T", "e0_y_c1"]
    meas_variables = ["e0_y_c1"]

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_variables

def get_model_bod():
    vl_original, model, _ = mopeds.examples.bod_model()
    vl_original["f"].variance = 1.5**2
    prediction_grid = {"x": [0, 10, 10]}
    meas_grid = {"x": [0, 10, 6]}

    def unfix_parameters(list_vl):
        for vl in list_vl:
            vl["theta1"].fixed = False
            vl["theta2"].fixed = False
        return list_vl

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, None

def get_model_poly1():
    vl_original, model = mopeds.examples.polynomial_1d()
    vl_original["y"].variance = 0.1**2
    prediction_grid = {"u": [0, 1, 10]}
    meas_grid = {"u": [0, 1, 5]}

    def unfix_parameters(list_vl):
        for vl in list_vl:
            vl["a"].fixed = False
            vl["b"].fixed = False
        return list_vl

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, None

def get_model_linear_example():
    vl_original, model = mopeds.examples.linear_example()
    vl_original["y"].variance = 0.5**2
    prediction_grid = {"u": [0, 1, 10], "v": [3, 4, 1]}
    meas_grid = {"u": [0, 1, 5], "v": [3, 4, 5]}

    def unfix_parameters(list_vl):
        for vl in list_vl:
            vl["a"].fixed = False
            vl["b"].fixed = False
            vl["c"].fixed = False
            vl["d"].fixed = False
        return list_vl

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, ["y", "z"]

def parameter_covariance_mopeds():
    vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_variables = model_selector()
    unfix_parameters = []
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            unfix_parameters.append(var.name)

    # true_parameters = {"a": 20}
    # unfix_parameters.pop(0)
    true_parameters = {}

    analyzer = mopeds.tools.ErrorAnalyzer(vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_variables, true_parameters=true_parameters)
    analyzer.parameter_covariance_mc(plot=False, num_samples=NUM)

    # analyzer.plot_estimation_accuracy()
    analyzer.model_prediction_error_mc(plot=False)
    analyzer.plot_model_prediction_MC(without_outliers=True)
    # v = analyzer.analyze_model_prediction()

    plt.show()

if __name__ == "__main__":
    parameter_covariance_mopeds()
