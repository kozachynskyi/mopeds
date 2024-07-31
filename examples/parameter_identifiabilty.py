import mopeds
import numpy as np
import pandas as pd
import copy
import matplotlib.pyplot as plt

MODEL_NAME = 5
NUM = 100

def model_selector():
    if MODEL_NAME == 1:
        return get_model_poly1()
    elif MODEL_NAME == 5:
        return get_model_dynamic_example()

def get_model_poly1():
    vl_original, model = mopeds.examples.polynomial_1d()
    vl_original["y"].variance = 0.1**2
    prediction_grid = {"u": [0, 1, 10]}
    meas_grid = {"u": [1, 1.00002, 20]}

    unfix_parameters = []
    for var in vl_original.values():
        if isinstance(var, mopeds.VariableParameter):
            unfix_parameters.append(var.name)
    # unfix_parameters = ["a"]
    # vl_original["b"].value = 2.01

    return vl_original, model, prediction_grid, meas_grid, unfix_parameters, None

def get_model_dynamic_example():
    vl_original, model, _ = mopeds.examples.yeast_growth()

    prediction_grid_controls = mopeds.tools.controls_grid_from_dict({"u1": [0.05, 0.2, 10], "u2": [5, 35, 1]})
    list_scenarios_prediction = []
    for controls in prediction_grid_controls:
        list_scenarios_prediction.append({"controls": controls, "time_grid": np.linspace(0,20,6), "initials": {}})

    measurement_grid_controls = mopeds.tools.controls_grid_from_dict({"u1": [0.05, 0.2, 3], "u2": [5, 35, 3]})
    list_scenarions_measurement = []
    for inits in [{"x1": 5}, {"x1": 7}]:
        for controls in measurement_grid_controls:
            list_scenarions_measurement.append({"controls": controls, "time_grid": np.linspace(0,20,6), "initials": inits})

    unfix_parameters = ["theta1", "theta2", "theta3", "theta4"]

    return vl_original, model, list_scenarios_prediction, list_scenarions_measurement, unfix_parameters, ["x1", "x2"]

def parameter_identifiability():
    vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_variables = model_selector()

    if MODEL_NAME == 5: 
        pe = mopeds.ParameterEstimation
    else:
        pe = mopeds.ParameterEstimationNLE

    analyzer = mopeds.tools.ErrorAnalyzer(vl_original, model, prediction_grid, meas_grid, unfix_parameters, meas_variables, pe_class=pe)

    if MODEL_NAME == 5:
        pe = analyzer.pe_main 
        pe.solver_settings["ipopt"]["linear_solver"] = "ma57"
        # pe.solver_settings["ipopt"]["tol"] = 1e-1
        # pe.solver_settings["ipopt"]["max_iter"] = 3000
        # pe.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"
    # analyzer.parameter_covariance_mc(100)
    analyzer.parameter_identifiability()
    # analyzer.plot_parameter_covariance_ellipse(normalize_parameters=False)
    # plt.show()

if __name__ == "__main__":
    parameter_identifiability()
