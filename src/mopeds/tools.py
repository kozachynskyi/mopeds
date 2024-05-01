""" Here come methods that use mopeds as import.
Separated from utilities to avoid dependency hell"""
from __future__ import annotations

import copy
from warnings import warn
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

import mopeds
from functools import cached_property
import tqdm


class ErrorAnalyzer():
    def __init__(self, variable_list, model, prediction_grid, measurement_grid, selected_parameters, measurement_names, *, rng=None, true_parameters=None, pe_class=None):
        if pe_class is None:
            self.PE_class = mopeds.ParameterEstimationNLE
        else:
            self.PE_class = pe_class

        self.variable_list = variable_list
        self.variable_list_true = copy.deepcopy(variable_list)

        if true_parameters is not None:
            for par_name, par_value in true_parameters.items():
                self.variable_list_true[par_name].value = par_value

        self.true_parameters = {}
        for var in self.variable_list_true.get_independent().values():
            if isinstance(var, mopeds.VariableParameter):
                self.true_parameters[var.name] = var.value[0]

        self._model = model
        self.prediction_grid = prediction_grid
        self.measurement_grid = measurement_grid
        self.selected_parameters = selected_parameters

        if measurement_names is None:
            self.measurement_names = model.varlist_algebraic(variable_list).keys()
        else:
            self.measurement_names = measurement_names

        if rng is None:
            self.rng = np.random.default_rng()
        else:
            self.rng = rng

    def unfix_parameters(self, list_varlist):
        for vl in list_varlist:
            for par_name in self.selected_parameters:
                vl[par_name].fixed = False
        return list_varlist

    @cached_property
    def pe_main(self) -> mopeds.ParameterEstimation:
        control_grid, true_params, measurement_names = generate_artificial_data_from_grid_nle(self._model, self.variable_list, self.measurement_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        control_grid = self.unfix_parameters(control_grid)
        return self.PE_class(self._model, control_grid)

    @cached_property
    def pe_artificial_data(self) -> mopeds.ParameterEstimation:
        true_data, true_params, _ = generate_artificial_data_from_grid_nle(self._model, self.variable_list_true, self.measurement_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        true_data = self.unfix_parameters(true_data)
        return self.PE_class(self._model, true_data)

    @cached_property
    def pe_prediction(self) -> mopeds.ParameterEstimation:
        prediction_data, true_params, measurement_names = generate_artificial_data_from_grid_nle(self._model, self.variable_list, self.prediction_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        prediction_data = self.unfix_parameters(prediction_data)
        return self.PE_class(self._model, prediction_data)

    @cached_property
    def pe_true_prediction(self) -> mopeds.ParameterEstimation:
        prediction_data, true_params, measurement_names = generate_artificial_data_from_grid_nle(self._model, self.variable_list_true, self.prediction_grid, perturbate=False, measurement_names=self.measurement_names, rng=self.rng)

        prediction_data = self.unfix_parameters(prediction_data)
        return self.PE_class(self._model, prediction_data)

        
    def get_s2_and_df(self, pe, parameters_dict):
        obj_and_residual = pe.calculate_objective_and_residual(parameters_dict, objective_function="ols")
        estimation_df = self.scale_df_all(pe, obj_and_residual["df_all"])
        estimation_df = estimation_df[pe.names_of_measurements]

        scaled_residuals = pe._unscale_resudials(obj_and_residual["residuals"])
        measurement_variance_estimate = np.diag(scaled_residuals.T @ scaled_residuals) / pe.dof

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
        return selected_indexes.to_numpy().all(axis=1)

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

    def scale_df_all(self, pe, df):
        varlist_alg = pe.model.varlist_algebraic(pe.list_input_varlist[0])
        for i, row in df.iterrows():
            df.iloc[i] = varlist_alg.scale_to_original(row)
        return df

    def model_prediction_error_mc(self, plot=False):
        self.list_predictions = []
        for index, row in tqdm.tqdm(self.df_params.iterrows()):
            prediction_i = self.pe_prediction.calculate_objective_and_residual(row.to_dict())["df_all"]
            prediction_df = self.scale_df_all(self.pe_prediction, prediction_i)
            self.list_predictions.append(prediction_df)

        if plot:
            self.plot_model_prediction()

    @property
    def threshold(self):
        return stats.t.ppf(0.975, self.pe_main.dof)

    def analyze_model_prediction(self):
        true_prediction_all = self.pe_true_prediction.calculate_objective_and_residual(self.true_parameters)["df_all"]
        true_prediction_all = self.scale_df_all(self.pe_true_prediction, true_prediction_all)
        true_prediction = true_prediction_all[self.measurement_names]

        list_prediction_quality = []
        list_prediction_quality_without_outliers = []
        metrics = {}
        self.list_prediction_std = []
        
        df_all_without_outliers, count_outliers = self.get_parameter_df(True)

        for meas_name in self.measurement_names:
            df_predictions = pd.concat(self.list_predictions, axis=1, keys=range(len(self.list_predictions)))
            df_predictions_without_outliers = df_predictions[df_all_without_outliers.index].swaplevel(axis=1)
            df_predictions = df_predictions.swaplevel(axis=1)

            res = stats.shapiro(df_predictions[meas_name], axis=1, nan_policy="omit")
            metrics[meas_name + "_shapiro_stat"] = res.statistic.mean()
            metrics[meas_name + "_shapiro_p"] = res.pvalue.mean()
            try:
                res = stats.shapiro(df_predictions_without_outliers[meas_name], axis=1, nan_policy="omit")
                metrics[meas_name + "_shapiro_stat_without_outliers"] = res.statistic.mean()
                metrics[meas_name + "_shapiro_p_without_outliers"] = res.pvalue.mean()
            except Exception:
                metrics[meas_name + "_shapiro_stat_without_outliers"] = np.nan
                metrics[meas_name + "_shapiro_p_without_outliers"] = np.nan

        for index, row in self.df_params.iterrows():
            try:
                df_predictions = self.list_predictions[index][self.measurement_names]
                cov_linearized = self.pe_main.calculate_sensitivity_and_fim_fast(row.to_dict())[2]

                jac_prediction = self.pe_prediction.calculate_sensitivity_and_fim_fast(row.to_dict())[0]
                prediction_std = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T)).reshape(self.pe_prediction.array_data.shape, order="F")
                lb = df_predictions - self.threshold * prediction_std
                ub = df_predictions + self.threshold * prediction_std
                in_bounds = ((lb <= true_prediction) & (ub >= true_prediction)).all()
                list_prediction_quality.append(in_bounds.all())
                if index in df_all_without_outliers.index:
                    list_prediction_quality_without_outliers.append(in_bounds.all())
                self.list_prediction_std.append(prediction_std)
            except Exception:
                pass

        array_quality = np.array(list_prediction_quality)
        array_quality_without_outliers = np.array(list_prediction_quality_without_outliers)
        metrics["predictions_in_bounds"] = array_quality.sum() / array_quality.shape[0]
        try:
            metrics["predictions_in_bounds_without_outliers"] = array_quality_without_outliers.sum() / array_quality_without_outliers.shape[0]
        except Exception:
            metrics["predictions_in_bounds_without_outliers"] = np.nan 

        return metrics


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

    def plot_parameter_prediction(self, parameters):
        import matplotlib.pyplot as plt
        true_values = self.pe_true_prediction.calculate_objective_and_residual(self.true_parameters)["df_all"]
        true_values = self.scale_df_all(self.pe_true_prediction, true_values)

        jac_prediction_all = self.pe_prediction.calculate_sensitivity_and_fim(parameters)["jac_sorted"]
        prediction_df = self.pe_prediction.calculate_objective_and_residual(parameters)["df_all"]
        prediction_df = self.scale_df_all(self.pe_prediction, prediction_df)

        cov_linearized = self.pe_main.calculate_sensitivity_and_fim(parameters)["cov_par"]
        
        threshold = stats.t.ppf(0.975, self.pe_main.dof)

        for meas_index, meas_name in enumerate(self.measurement_names):
            fig = plt.figure()

            plt.plot(true_values, true_values, label="true", c="black", ls="dashed")

            jac_prediction = jac_prediction_all[meas_name]

            prediction_linearized = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T))
            prediction_line = prediction_df[meas_name]
            plt.plot(true_values, prediction_line, label="prediction", c="b")
            plt.plot(true_values, prediction_line + threshold*prediction_linearized,  label="95% CI", c="r")
            plt.plot(true_values, prediction_line - threshold*prediction_linearized,label="95% CI", c="r")

            plt.legend()
            plt.title(meas_name)

            fig.suptitle("Predicted model accuracy")

    def plot_model_prediction_MC(self, *, without_outliers=False):
        import matplotlib.pyplot as plt
        true_values_all = self.pe_true_prediction.calculate_objective_and_residual(self.true_parameters)["df_all"]
        true_values_all = self.scale_df_all(self.pe_true_prediction, true_values_all)

        df_predictions = pd.concat(self.list_predictions, axis=1, keys=range(len(self.list_predictions)))
        if without_outliers:
            df_params, count_outliers = self.get_parameter_df(without_outliers)
            if df_params.empty:
                return

            df_predictions = df_predictions[df_params.index]

        df_predictions = df_predictions.swaplevel(axis=1)

        try:
            cov_linearized = self.pe_artificial_data.calculate_sensitivity_and_fim(self.true_parameters)["cov_par"]
            jac_prediction = self.pe_true_prediction.calculate_sensitivity_and_fim(self.true_parameters)["jac_full"]
            prediction_std_all = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T)).reshape(self.pe_prediction.array_data.shape, order="F")
        except Exception:
            pass

        fig, axis = plt.subplots(nrows=1, ncols=len(self.measurement_names), squeeze=False)

        for meas_index, (meas_name, ax) in enumerate(zip(self.measurement_names, axis.flat)):
            true_values = true_values_all[meas_name]
            min = df_predictions[meas_name].min(axis=1)
            max = df_predictions[meas_name].max(axis=1)

            ax.plot(min - true_values, label="MC, min", c="g")
            ax.plot(max - true_values, label="MC, max", c="g")
            ax.plot([0]*true_values.shape[0], label="true", c="black", ls="dashed")

            try:
                prediction_std = prediction_std_all[:, meas_index]
                ax.plot(-self.threshold*prediction_std, label="Lin, CI 95%", c="red")
                ax.plot(self.threshold*prediction_std, label="Lin, CI 95%", c="red")
            except Exception:
                pass

            ax.legend()
            ax.set_title(meas_name)

        fig_name = "Predicted model accuracy, MC" 
        if without_outliers:
            fig_name = fig_name + ", without outliers"
            fig.supxlabel(f"Outliers = {count_outliers}")

        fig.suptitle(fig_name)


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
    variable_list: mopeds.VariableList,
    model: mopeds.Model,
    time_grid: np.ndarray,
    algebraic: bool = False,
    perturbate: bool = False,
    rng: np.random.Generator | None = None
) -> mopeds.VariableList:
    if rng is None:
        rng = np.random.default_rng()
    # Simulated ODE/DAE and replaces StateVariable values with simulated data
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    sim = mopeds.Simulator(model, time_grid, var_list_fixed)
    var_list_exp = sim.simulate(algebraic=True)[2]

    # Replace empty state variables with results from simulation
    variable_list_with_data = copy.deepcopy(variable_list)
    for key, var in var_list_exp.items():
        if not isinstance(var, mopeds.VariableControl):
            df = var.dataframe
            if perturbate:
                std = var_list_fixed[key].variance ** 0.5
                value = rng.normal(var.dataframe, std)
                df[key] = value

            variable_list_with_data[key].dataframe = df

    return variable_list_with_data

def generate_artificial_data_from_grid_nle(
    model: mopeds.Model,
    variable_list: mopeds.VariableList,
    control_bounds: dict[list[float]],
    perturbate: bool = True,
    rng: np.random.Generator = None,
    measurement_names: list[str] = None,
    *,
    keep_in_bounds: bool = True,
) -> tuple[list[mopeds.VariableList], dict[str, float]]:
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
    model: mopeds.Model,
    variable_list: mopeds.VariableList,
    controls: list[dict[float]],
    perturbate: bool = True,
    rng: np.random.Generator = None,
    measurement_names: list[str] = None,
    *,
    keep_in_bounds: bool = True,
    ) -> tuple[list[mopeds.VariableList], dict[str, float]]:
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
        if isinstance(var, mopeds.VariableParameter):
            var_varlist = variable_list_original[var.name]
            true_parameters[var_varlist.name] = var_varlist.value[0]

    sim_fixed = mopeds.SimulatorNLE(model, variable_list)

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
) -> tuple[list[mopeds.VariableList], dict[str, float]]:
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
