""" Here come methods that use mopeds as import.
Separated from utilities to avoid dependency hell"""
from __future__ import annotations

import copy
from warnings import warn

import numpy as np
import pandas as pd

from mopeds import Model, Simulator, SimulatorNLE, VariableList, VariableParameter, VariableControl


class ErrorAnalyzer():
    def __init__(self, variable_list, model, prediction_grid, measurement_grid, selected_parameters, measurement_names, *, rng=None):
        self.variable_list = variable_list
        self._model = model
        self.prediction_grid = prediction_grid
        self.measurement_grid = measurement_grid
        self.selected_parameters = selected_parameters
        self.measurement_names = measurement_names

        self.rng = rng
        if rng is None:
            self.rng = np.random.default_rng()

        self.create_pe_objects()

    def unfix_parameters(self, list_varlist):
        for vl in list_varlist:
            for par_name in self.selected_parameters:
                vl[par_name].fixed = False
        return list_varlist

    def create_pe_objects(self):
        from mopeds import ParameterEstimationNLE

        prediction_data, true_params, measurement_names = generate_artificial_data_from_grid_nle(self._model, self.variable_list, self.prediction_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        # Rewrite measurement names, if it was None
        self.measurement_names = measurement_names

        prediction_data = self.unfix_parameters(prediction_data)

        self.pe_prediction = ParameterEstimationNLE(self._model, prediction_data)

        true_data, true_params, _ = generate_artificial_data_from_grid_nle(self._model, self.variable_list, self.measurement_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        true_data = self.unfix_parameters(true_data)

        self.pe_artificial_data = ParameterEstimationNLE(self._model, true_data)
        self.true_parameters = true_params

    def parameter_covariance_mc(self, num_samples=100, plot=False):
        original_data = copy.deepcopy(self.pe_artificial_data.array_data)
        pe = self.pe_artificial_data

        list_parameters = []
        for i in range(num_samples):
            pe.array_data = self.rng.normal(original_data, (1 / pe.array_inverted_scaled_std))
            res = pe.optimize(None, "wls", direct_optimization=True, reuse_solver=True)
            if pe.solver.stats()["success"]:
                list_parameters.append(list(res["x_dict"].values()))

        self.df_params = pd.DataFrame(list_parameters, columns=pe.varlist_decision.keys())
        if plot:
            self.plot_covariance()

    def plot_covariance(self):
        import matplotlib.pyplot as plt

        cov_linearized = self.pe_artificial_data.calculate_sensitivity_and_fim(self.true_parameters)["cov_par"]
        std_linearized = np.sqrt(np.diag(cov_linearized))
        cov_mc = self.df_params.cov()

        axis = self.df_params.hist(bins=int(self.df_params.shape[0]/20))
        for index, (ax, val) in enumerate(zip(axis.flat, self.true_parameters.values())):
            ax.axvline(val, 0, ax.yaxis.get_data_interval()[1], c="r")
            ax.axvline(val + 2*std_linearized[index], 0, ax.yaxis.get_data_interval()[1], c="r")
            ax.axvline(val - 2*std_linearized[index], 0, ax.yaxis.get_data_interval()[1], c="r")

        plt.figure()
        plt.scatter(self.df_params.iloc[:, 0], self.df_params.iloc[:, 1])

        return axis

    def scale_df_all(self, pe, df):
        varlist_alg = pe.model.varlist_algebraic(pe.list_input_varlist[0])
        for i, row in df.iterrows():
            df.iloc[i] = varlist_alg.scale_to_original(row)
        return df

    def model_prediction_error_mc(self, plot=False, rng=None):
        self.list_predictions = []
        for index, row in self.df_params.iterrows():
            prediction_i = self.pe_prediction.calculate_objective_and_residual(row.to_dict())["df_all"]
            prediction_df = self.scale_df_all(self.pe_prediction, prediction_i)
            self.list_predictions.append(prediction_df)

        if plot:
            self.plot_model_prediction()

    def plot_model_prediction(self):
        import matplotlib.pyplot as plt

        jac_prediction_all = self.pe_prediction.calculate_sensitivity_and_fim(self.true_parameters)["jac_sorted"]
        prediction_df = self.pe_prediction.calculate_objective_and_residual(self.true_parameters)["df_all"]
        prediction_df = self.scale_df_all(self.pe_prediction, prediction_df)

        prediction_original_model = self.pe_artificial_data.calculate_objective_and_residual(self.df_params.iloc[-1].to_dict())["df_all"]
        prediction_original_model = self.scale_df_all(self.pe_artificial_data, prediction_original_model)

        cov_linearized = self.pe_artificial_data.calculate_sensitivity_and_fim(self.df_params.iloc[-1].to_dict())["cov_par"]

        plt.figure()

        for meas_index, meas_name in enumerate(self.measurement_names):
            df_predictions = pd.DataFrame(np.array(self.list_predictions)[:,:, meas_index])

            data = self.pe_artificial_data.array_data_unscaled[:, meas_index]
            jac_prediction = jac_prediction_all[meas_name]
            prediction_line = prediction_df[meas_name]

            prediction_linearized = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T))
            std_mc = df_predictions.std()


            plt.scatter(prediction_original_model[meas_name], data, ls="", c="black")
            plt.plot(prediction_line, prediction_line, c="b")
            plt.plot(prediction_line + 2*prediction_linearized, prediction_line, label="lin", c="r")
            plt.plot(prediction_line - 2*prediction_linearized, prediction_line, label="lin", c="r")
            plt.plot(prediction_line + 2*std_mc, prediction_line, label="mc", c="g")
            plt.plot(prediction_line - 2*std_mc, prediction_line, label="mc", c="g")
            plt.legend()


def create_grid(bounds: list[list[float]]) -> list[list[float]]:
    """Create a grid in a given bounds. Bounds is dictionary, with variable names as keys(),
    and values() as a list with 3 elements: [lower_bound, upper_bound, num_points]
    """
    linspace_list = []
    for bound in bounds:
        linspace_list.append(np.linspace(start=bound[0], stop=bound[1], num=bound[2]))
    meshgrid = np.meshgrid(*linspace_list)

    grid = [n_grid.ravel() for n_grid in meshgrid]
    grid = np.array(grid).transpose().tolist()
    return grid, meshgrid


def generate_varlist_with_data(
    variable_list: VariableList,
    model: Model,
    time_grid: np.ndarray,
    algebraic: bool = False,
    perturbate: bool = False,
    rng: np.random.Generator | None = None
) -> VariableList:
    if rng is None:
        rng = np.random.default_rng()
    # Simulated ODE/DAE and replaces StateVariable values with simulated data
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    sim = Simulator(model, time_grid, var_list_fixed)
    var_list_exp = sim.simulate(algebraic=True)[2]

    # Replace empty state variables with results from simulation
    variable_list_with_data = copy.deepcopy(variable_list)
    for key, var in var_list_exp.items():
        if not isinstance(var, VariableControl):
            df = var.dataframe
            if perturbate:
                std = var_list_fixed[key].variance ** 0.5
                value = rng.normal(var.dataframe, std)
                df[key] = value

            variable_list_with_data[key].dataframe = df

    return variable_list_with_data

def generate_artificial_data_from_grid_nle(
    model: Model,
    variable_list: VariableList,
    control_bounds: dict[list[float]],
    perturbate: bool = True,
    rng: np.random.Generator = None,
    measurement_names: list[str] = None,
    *,
    keep_in_bounds: bool = True,
) -> tuple[list[VariableList], dict[str, float]]:
    """Wrapper around generate_artificial_data_nle, where controls are generated based on uniform grid.
    control_bounds is dictionary, with variable names as keys(),
    and values() as a list with 3 elements: [lower_bound, upper_bound, num_points]
    """
    grid, _ = create_grid(list(control_bounds.values()))
    control_grid = []
    for grid_point in grid:
        control_grid.append(dict(zip(control_bounds.keys(), grid_point)))
    return generate_artificial_data_nle(model, variable_list, control_grid, perturbate, rng, measurement_names, keep_in_bounds=keep_in_bounds)


def generate_artificial_data_nle(
    model: Model,
    variable_list: VariableList,
    controls: list[dict[float]],
    perturbate: bool = True,
    rng: np.random.Generator = None,
    measurement_names: list[str] = None,
    *,
    keep_in_bounds: bool = True,
    ) -> tuple[list[VariableList], dict[str, float]]:
    """Generate artificial data that can immediately be used by Parameter Estimator.
    Returns list of varlists and a dictionary with parameter values that were used
    to generate data.

    Parameter values that are used are taken from variable list.
    controls is a list with dictionary, representing control variables and their respective values.
    For example, [{"P": 1, "T": 80}, {"P": 2, "T": 90}] will create two sets of artificial data, generated
    at pressure 1 and T 80 and pressure 2 and temperature 90.
    perturbate: if True, generated data is perturbated based on variance in variable_list
    rng: is a rng object, user can use it to predefine the randomization of the noise
    measurement_names: list with variable names, for which artificial data should be generated
    keep_in_bounds: if perturbated value is out of variable bounds -> make it equal to the closes bound
    """
    if rng is None:
        rng = np.random.default_rng()

    if measurement_names is None:
        measurement_names = model.varlist_algebraic(variable_list).keys()

    variable_list_original = copy.deepcopy(variable_list)
    true_parameters = {}

    for var in variable_list.values():
        var.fixed = True

    for var in model.varlist_independent(variable_list).values():
        if isinstance(var, VariableParameter):
            var_varlist = variable_list_original[var.name]
            true_parameters[var_varlist.name] = var_varlist.value[0]

    sim_fixed = SimulatorNLE(model, variable_list)

    varlist_list = []
    for grid_point in controls:
        variable_list_optimizer = copy.deepcopy(variable_list_original)
        sim_fixed.change_independent_variables(grid_point)
        varlist_results = sim_fixed.simulate()[2]

        # Set startings values
        for variable_name, variable in varlist_results.items():
            if variable_name in measurement_names:
                value = variable.value[0]
                if perturbate:
                    value = rng.normal(value, np.sqrt(variable.variance))
                    if keep_in_bounds:
                        value = min(variable.upper_bound, max(variable.lower_bound, value))
                variable_list_optimizer[variable_name].guess = value
                variable_list_optimizer[variable_name].value = value
        for var_name, var_value in grid_point.items():
            variable_list_optimizer[var_name].value = var_value
        varlist_list.append(variable_list_optimizer)

    return varlist_list, true_parameters, measurement_names

def generate_varlist_with_data_NLE(
    model,
    variable_list,
    control_bounds,
    perturbate: bool = True,
    rng: np.random.Generator = None,
    measurement_names: list[str] = None,
) -> tuple[list[VariableList], dict[str, float]]:
    warn("Deprecated API, use generate_artificial_data_from_grid_nle", FutureWarning, 2)
    return generate_artificial_data_from_grid_nle(model, variable_list, control_bounds, perturbate, rng, measurement_names)


def analyze_scaling_nle(model, varlist, control_bounds):
    """Change the control variables in a given bounds, calculate all algeraic variables and provide lower and upper bounds for them"""
    all_data, true_parameters = generate_varlist_with_data_NLE(model, varlist, control_bounds=control_bounds, perturbate=False)
    v = all_data[0].dataframe
    vv = all_data[1].dataframe
    g = pd.concat([vl.dataframe for vl in all_data])

    algebraic_names = varlist.get_algebraic().keys()
    selected_data = g[algebraic_names]
    return_bounds = dict(zip(algebraic_names, zip(selected_data.min(), selected_data.max())))
    return return_bounds
