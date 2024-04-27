""" Here come methods that use mopeds as import.
Separated from utilities to avoid dependency hell"""
from __future__ import annotations

import copy
from warnings import warn
from itertools import combinations

import numpy as np
import pandas as pd

from mopeds import Model, Simulator, SimulatorNLE, VariableList, VariableParameter, VariableControl


class ErrorAnalyzer():
    def __init__(self, variable_list, model, prediction_grid, measurement_grid, selected_parameters, measurement_names, *, rng=None, true_parameters=None, pe_class=None):
        if pe_class is None:
            from mopeds import ParameterEstimationNLE
            self.PE_class = ParameterEstimationNLE
        else:
            self.PE_class = pe_class

        self.variable_list = variable_list
        self.variable_list_true = copy.deepcopy(variable_list)

        if true_parameters is not None:
            for par_name, par_value in true_parameters.items():
                self.variable_list_true[par_name].value = par_value

        self.true_parameters = {}
        for var in self.variable_list_true.get_independent().values():
            if isinstance(var, VariableParameter):
                self.true_parameters[var.name] = var.value[0]

        self._model = model
        self.prediction_grid = prediction_grid
        self.measurement_grid = measurement_grid
        self.selected_parameters = selected_parameters
        self.measurement_names = measurement_names

        # Used for parameter estimation during MC
        self.pe_main = None
        # Used for generation of the perturbated for MC
        self.pe_artificial_data = None
        # Used for calculation of predicted values
        self.pe_prediction = None

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

        prediction_data, true_params, measurement_names = generate_artificial_data_from_grid_nle(self._model, self.variable_list, self.prediction_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        # Rewrite measurement names, if it was None
        self.measurement_names = measurement_names

        prediction_data = self.unfix_parameters(prediction_data)

        self.pe_prediction = self.PE_class(self._model, prediction_data)

        true_data, true_params, _ = generate_artificial_data_from_grid_nle(self._model, self.variable_list_true, self.measurement_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        true_data = self.unfix_parameters(true_data)

        self.pe_artificial_data = self.PE_class(self._model, true_data)

        control_grid, true_params, measurement_names = generate_artificial_data_from_grid_nle(self._model, self.variable_list, self.measurement_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        control_grid = self.unfix_parameters(control_grid)
        self.pe_main = self.PE_class(self._model, control_grid)
        
    def get_s2_and_df(self, pe, parameters_dict):
        obj_and_residual = pe.calculate_objective_and_residual(parameters_dict, objective_function="ols")
        estimation_df = self.scale_df_all(pe, obj_and_residual["df_all"])
        estimation_df = estimation_df[pe.names_of_measurements]

        scaled_residuals = pe._unscale_resudials(obj_and_residual["residuals"])
        dof = np.count_nonzero(pe.array_data_mask) - (len(pe.varlist_decision) / len(pe.names_of_measurements))
        measurement_variance_estimate = np.diag(scaled_residuals.T @ scaled_residuals) / dof

        return np.sqrt(measurement_variance_estimate), estimation_df

    def parameter_covariance_mc(self, num_samples=100, plot=False):
        original_data = copy.deepcopy(self.pe_artificial_data.array_data)

        pe = self.pe_main

        list_parameters = []
        self.list_estimation = []
        self.failed_pes = 0
        list_s2 = []
        list_s2_true = []

        for i in range(num_samples):
            perturbated_data = self.rng.normal(original_data, (1 / pe.array_inverted_scaled_std))
            pe.array_data = perturbated_data
            res = pe.optimize(None, "wls", direct_optimization=True, reuse_solver=True)

            if pe.solver.stats()["success"]:
                list_parameters.append(list(res["x_dict"].values()))

                s2_esimated, df_estimated = self.get_s2_and_df(self.pe_main, res["x_dict"])
                list_s2.append(s2_esimated)
                self.list_estimation.append(df_estimated)

                self.pe_artificial_data.array_data = perturbated_data
                s2_true, _ = self.get_s2_and_df(self.pe_artificial_data, self.true_parameters)
                list_s2_true.append(s2_true)

            else:
                self.failed_pes += 1

        self.pe_artificial_data.array_data = copy.deepcopy(original_data)

        self.df_s2 = pd.DataFrame(list_s2, columns=pe.names_of_measurements)
        self.df_s2_true = pd.DataFrame(list_s2_true, columns=pe.names_of_measurements)
        self.df_params = pd.DataFrame(list_parameters, columns=pe.varlist_decision.keys())

        if plot:
            self.plot_parameter_covariance(normalize_parameters=True)
            self.plot_parameter_covariance(normalize_parameters=False)
            self.plot_parameter_variance()
            self.plot_estimation_accuracy()

    @property
    def last_estimated_parameters(self):
        return self.df_params.iloc[0].to_dict()

    @property
    def df_params_normalized(self):
        """Source min-max scaling of sklearn
        https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html#sklearn-preprocessing-minmaxscaler
        """
        X = self.df_params
        min, max = (-1, 1)
        X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
        X_scaled = X_std * (max - min) + min

        X_pure = self.df_params[self.no_outliers]
        min, max = (-1, 1)
        X_std_pure = (X_pure - X_pure.min(axis=0)) / (X_pure.max(axis=0) - X_pure.min(axis=0))
        X_scaled_pure = X_std_pure * (max - min) + min
        return X_scaled, X_scaled_pure

    @property
    def no_outliers(self):
        real_s2 = 1 / self.pe_artificial_data.array_inverted_std
        normalized = self.df_s2 - real_s2[0, :]
        bound = 3 * self.df_s2_true.std()

        selected_indexes = (normalized >= -bound) & ( normalized <= bound)
        return selected_indexes.to_numpy()

    @property
    def bins_number(self):
        return max(1, int(self.df_params.shape[0]* 0.2))

    def get_parameter_df(self, without_outliers):
        if without_outliers:
            no_outliers = self.no_outliers
            count_outliers = no_outliers.shape[0] - no_outliers.sum()
            df = self.df_params[self.no_outliers]
        else:
            count_outliers = 0
            df = self.df_params
        return df, count_outliers

    def plot_parameter_variance(self, *, without_outliers=False, parameters=None):
        df, count_outliers = self.get_parameter_df(without_outliers)

        axis = self.df_params.hist(bins=self.bins_number)
        fig = axis.flat[0].get_figure()

        if parameters is not None:
            try:
                cov_linearized = self.pe_main.calculate_sensitivity_and_fim(self.last_estimated_parameters)["cov_par"]
                std_linearized = np.sqrt(np.diag(cov_linearized))
                for index, (ax, val) in enumerate(zip(axis.flat, parameters)):
                    ax.axvline(val, 0, ax.yaxis.get_data_interval()[1], c="r")
                    ax.axvline(val + 2*std_linearized[index], 0, ax.yaxis.get_data_interval()[1], c="r")
                    ax.axvline(val - 2*std_linearized[index], 0, ax.yaxis.get_data_interval()[1], c="r")

            except Exception:
                pass

        fig_name = "Parameter Variance MC"
        if without_outliers:
            fig_name = fig_name + ", without outliers"
            fig.supxlabel(f"Outliers = {count_outliers}")

        fig.suptitle(fig_name)

    def plot_parameter_covariance(self, *, normalize_parameters=True, without_outliers=False):
        import matplotlib.pyplot as plt
        from matplotlib.axes import Axes

        num_par = len(self.pe_prediction.varlist_decision)
        fig, axes = plt.subplots(ncols=num_par - 1, nrows=num_par - 1, layout="constrained")
        if isinstance(axes, Axes) == 1:
            axes = [axes]

        comb = list(combinations(range(num_par), 2))
        par_names = list(self.pe_main.varlist_decision.keys())

        lb = np.asarray(self.pe_main.varlist_decision.scale_to_original(self.pe_main.lower_bound))
        ub = np.asarray(self.pe_main.varlist_decision.scale_to_original(self.pe_main.upper_bound))

        if normalize_parameters:
            df_normalized_all, df_normalized_no_outliers = self.df_params_normalized
            if without_outliers:
                df = df_normalized_no_outliers
                count_outliers = df_normalized_all.shape[0] - df_normalized_no_outliers.shape[0]
            else:
                df = df_normalized_all
        else:
            df, count_outliers = self.get_parameter_df(without_outliers)

        for (par1_index, par2_index) in comb:
            if len(axes) == 1:
                ax = axes[0]
            else:
                ax = axes[par2_index - 1, par1_index]

            par_1_df = df.iloc[:, par1_index]
            par_2_df = df.iloc[:, par2_index]

            ax.scatter(par_1_df, par_2_df)

            if np.isclose(par_1_df, lb[par1_index]).any():
                ax.axvline(par_1_df.min(), 0, 1, c="g")
            if np.isclose(par_1_df, ub[par1_index]).any():
                ax.axvline(par_1_df.max(), 0, 1, c="g")
            if np.isclose(par_2_df, lb[par2_index]).any():
                ax.axhline(par_2_df.min(), 0, 1, c="g")
            if np.isclose(par_2_df, ub[par2_index]).any():
                ax.axhline(par_2_df.max(), 0, 1, c="g")

            if par1_index == 0:
                ax.set_ylabel(f"{par_names[par2_index]}")

            if par2_index == len(par_names) - 1:
                ax.set_xlabel(f"{par_names[par1_index]}")

        fig_name = "Parameter Covariance MC"
        if normalize_parameters:
            fig_name = fig_name + ", normalized values"
        if without_outliers:
            fig_name = fig_name + ", without outliers"
            fig.supxlabel(f"Outliers = {count_outliers}")

        fig.suptitle(fig_name)

    def calculate_cosine_simmilarity(self):
        raise NotImplementedError
        from sklearn.metrics.pairwise import cosine_similarity
        jac_prediction_all = self.pe_prediction.calculate_sensitivity_and_fim(self.true_parameters)["jac_full"]
        prediction_df = self.pe_prediction.calculate_objective_and_residual(self.true_parameters)["df_all"]
        prediction_df = self.scale_df_all(self.pe_prediction, prediction_df)
        list_similarity = []

        list_std_mc = []
        for meas_index, meas_name in enumerate(self.measurement_names):
            df_index = self.pe_prediction.list_simulators[0].mapping_algebraic_variables[meas_name]
            df_predictions = pd.DataFrame(np.array(self.list_predictions)[:,:, df_index])
            v = df_predictions.std() + df_predictions.mean()
            for index, row in self.df_params.iterrows():
                jac_prediction = jac_prediction_all[meas_name]
                cov_linearized = self.pe_main.calculate_sensitivity_and_fim(row.to_dict())["cov_par"]
                prediction_linearized = np.sqrt(np.diag(jac_prediction_all @ cov_linearized @ jac_prediction_all.T))

            list_similarity.append(cosine_similarity([std_mc], [prediction_linearized]))

        self.cosine_similarity = np.array(list_similarity)

    def scale_df_all(self, pe, df):
        varlist_alg = pe.model.varlist_algebraic(pe.list_input_varlist[0])
        for i, row in df.iterrows():
            df.iloc[i] = varlist_alg.scale_to_original(row)
        return df

    def model_prediction_error_mc(self, plot=False):
        self.list_predictions = []
        for index, row in self.df_params.iterrows():
            prediction_i = self.pe_prediction.calculate_objective_and_residual(row.to_dict())["df_all"]
            prediction_df = self.scale_df_all(self.pe_prediction, prediction_i)
            self.list_predictions.append(prediction_df)

        if plot:
            self.plot_model_prediction()

    def plot_parameter_prediction(self):
        import matplotlib.pyplot as plt
        prediction_original_model = self.pe_main.calculate_objective_and_residual(self.last_estimated_parameters)["df_all"]
        prediction_original_model = self.scale_df_all(self.pe_main, prediction_original_model)

        fig = plt.figure()
        for meas_index, meas_name in enumerate(self.measurement_names):
            data = self.pe_artificial_data.array_data_unscaled[:, meas_index]
            plt.scatter(prediction_original_model[meas_name], data, ls="", c="black")
            plt.plot(*[[data.min(), data.max()]]*2)
        fig.suptitle("Parameter prediction")

    def plot_estimation_accuracy(self, *, without_outliers=False):
        if without_outliers:
            no_outliers = self.no_outliers
            count_outliers = no_outliers.shape[0] - no_outliers.sum()
            df = self.df_s2[no_outliers]

        else:
            df = self.df_s2

        axis = df.hist(bins=self.bins_number)
        fig = axis.flat[0].get_figure()
        std_s2 = df.std()
        mean_s2 = df.mean()
        real_s2 = 1 / self.pe_artificial_data.array_inverted_std

        for ax, val, std, std_real in zip(axis.flat, mean_s2, std_s2, real_s2[0,:]):
            ax.axvline(std_real, 0, 1, c="g")
            ax.axvline(val, 0, 1, c="r")
            ax.axvline(val + 2*std, 0, 1, c="r")
            ax.axvline(val - 2*std, 0, 1, c="r")
            ax.set_title(ax.get_title() + f"\nReal s2 {std_real}, estimated {round(val, 5)}")
        
        fig_name = "Estimation accuracy of parameters"

        if without_outliers:
            fig_name = fig_name + ", without outliers"
        fig.suptitle(fig_name)

        if without_outliers:
            fig.supxlabel(f"Outliers = {count_outliers}")

    def plot_model_prediction(self):
        import matplotlib.pyplot as plt

        jac_prediction_all = self.pe_prediction.calculate_sensitivity_and_fim(self.true_parameters)["jac_sorted"]
        prediction_df = self.pe_prediction.calculate_objective_and_residual(self.true_parameters)["df_all"]
        prediction_df = self.scale_df_all(self.pe_prediction, prediction_df)

        cov_linearized = self.pe_main.calculate_sensitivity_and_fim(self.last_estimated_parameters)["cov_par"]
        # cov_linearized = np.identity(len(self.selected_parameters)) * cov_linearized


        for meas_index, meas_name in enumerate(self.measurement_names):
            fig = plt.figure()
            df_index = self.pe_prediction.list_simulators[0].mapping_algebraic_variables[meas_name]
            df_predictions = pd.DataFrame(np.array(self.list_predictions)[:,:, df_index])

            jac_prediction = jac_prediction_all[meas_name]
            prediction_line = prediction_df[meas_name]

            prediction_linearized = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T))
            std_mc = df_predictions.std()

            plt.plot(prediction_line, prediction_line, c="b")
            plt.plot(prediction_line + 2*prediction_linearized, prediction_line, label="lin", c="r")
            plt.plot(prediction_line - 2*prediction_linearized, prediction_line, label="lin", c="r")
            plt.plot(prediction_line + 2*std_mc, prediction_line, label="mc", c="g")
            plt.plot(prediction_line - 2*std_mc, prediction_line, label="mc", c="g")
            plt.title(meas_name)
            plt.legend()

            fig.suptitle("Parameter prediction accuracy")


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
