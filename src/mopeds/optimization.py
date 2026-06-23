from __future__ import annotations

import copy
import logging
from abc import abstractmethod
from collections.abc import Callable
from itertools import combinations
from typing import Sequence
from warnings import warn
import itertools
from functools import partial, cached_property

import casadi as ca
import numpy as np
import pandas as pd
from scipy import linalg
from tqdm import tqdm

from mopeds import (
    Model,
    Simulator,
    SimulatorNLE,
    VariableAlgebraic,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableList,
    VariableParameter,
    VariableState,
    tools,
    utilities,
    get_options,
    _consistent_scaling_decorator,
    options,
)


def eigsorted(cov):
    vals, vecs = np.linalg.eig(cov)
    order = np.flip(vals.argsort())
    return vals[order], vecs[:, order]


class Optimizer(object):
    def __init__(
        self,
        model: Model,
        variable_lists: list[VariableList],
        simulator_name: str,
        simulator_settings: dict | None,
    ) -> None:
        self.solver_name: str
        self.solver_settings: dict
        self._created_with_options = get_options()

        if not isinstance(variable_lists, list):
            raise (Exception("Variable list should be nested of type list"))
        self.logger: logging.Logger = logging.getLogger(__name__)
        self.model: Model = model

        # Deepcopy is used to avoid manipulating input variable list
        self.list_input_varlist: list[VariableList] = copy.deepcopy(variable_lists)
        for varlist in self.list_input_varlist:
            self.model.subsitute_casadi_symbols(varlist)

        # Each varlist holds respective variables
        self.varlist_decision: VariableList = VariableList()
        self.varlist_decision_direct: VariableList = VariableList()
        self.varlist_parameter: VariableList = VariableList()
        self.varlist_control: VariableList = VariableList()
        self.varlist_state: VariableList = VariableList()
        self.varlist_algebraic: VariableList = VariableList()
        self.nlpsol_g: None | ca.MX = None

        self.simulator_name: str = simulator_name
        self.simulator_settings: dict | None = simulator_settings

        self.list_simulators: Sequence[Simulator | SimulatorNLE]

    @abstractmethod
    def optimize(self, scale):
        # Runs optimization once
        raise (NotImplementedError)

    @_consistent_scaling_decorator
    def variables_dict_to_list(
        self, variables_dict: dict[str, float], *, scaling=None
    ) -> list[float]:
        """Takes a dictionary with {"var_name": var_value} and transforms to list
        corresponding to the order of self.varlist_decision variables"""
        if not isinstance(scaling, bool):
            scaling = self._created_with_options["variable_scaling"]

        with options(variable_scaling=scaling):
            selected_variables: list[float] = []
            for var_name in variables_dict.keys():
                if var_name not in self.varlist_decision.keys():
                    if "time_sp" in var_name or "weight_" in var_name:
                        raise ValueError(
                            f"Variable {var_name} is not a decision variable!"
                        )
                    print(f"Supplied value for variables {var_name} is ignored!")
            for var_name, var in self.varlist_decision.items():
                try:
                    scaled_value = var.scale_from_original(variables_dict[var_name])
                    selected_variables.append(scaled_value)
                except KeyError:
                    raise KeyError(f"Missing value for {var_name}")

        return selected_variables

    @_consistent_scaling_decorator
    def _setup_initialization(self) -> None:
        """Sets initials and bounds for optimizer, and as default no scaling.
        If guess equals 0, 1 is used instead to avoid division by 0 during initialization"""
        guess = []
        lower_bound = []
        upper_bound = []

        for var in self.varlist_decision.values():
            if var.guess == 0:
                guess.append(1.0)
            else:
                guess.append(var.scale_from_original(var.guess))
            lower_bound.append(var.scale_from_original(var.lower_bound))
            upper_bound.append(var.scale_from_original(var.upper_bound))

        self.guess: np.ndarray = np.array(guess)
        self.lower_bound: np.ndarray = np.array(lower_bound)
        self.upper_bound: np.ndarray = np.array(upper_bound)
        self.logger.debug(
            f"Initialized:\nguess {self.guess}\nlower_bound {self.lower_bound}\nupper_bound {self.upper_bound}"
        )

    @abstractmethod
    def _objective(self) -> tuple[ca.MX | ca.DM, ca.MX | ca.DM]:
        """Returns a way to calculate and objective. Dependent on optimization type."""
        raise (NotImplementedError)

    def check_result_bounds(self, result) -> pd.DataFrame:
        """Use output of self.optimize() and check if decision variables are not on the bounds
        Returns dataframe with variables that are not in bounds."""
        if result["x_unscaled"].shape[0] == len(self.varlist_decision):
            varlist_decision = self.varlist_decision
            lower_bound = self.lower_bound
            upper_bound = self.upper_bound
        elif result["x_unscaled"].shape[0] == len(self.varlist_decision_direct):
            varlist_decision = self.varlist_decision_direct
            lower_bound = self.lower_bound_direct
            upper_bound = self.upper_bound_direct
        else:
            raise ValueError

        results_with_bound = {}
        for i, (var, res_i) in enumerate(
            zip(varlist_decision.values(), result["x_unscaled"])
        ):
            value = var.scale_to_original(float(res_i))
            lb = var.scale_to_original(lower_bound[i])
            ub = var.scale_to_original(upper_bound[i])
            results_with_bound[var.name] = [lb, value, ub]
        df = pd.DataFrame.from_dict(results_with_bound).T
        selected_df = df[(df[0] >= df[1]) | (df[2] <= df[1])]
        selected_df.columns = ["lb", "value", "ub"]
        return selected_df

    def _change_guess(self, guess_values):
        """Overwrite the guess correcly both in self.guess and self.guess_direct"""
        arr = np.asarray(guess_values)
        if arr.ndim != 1:
            raise ValueError

        num_guess = arr.shape[0]
        self.guess = guess_values
        if self.model.equations_differential is None:
            self.guess_direct[:num_guess] = guess_values

    @_consistent_scaling_decorator
    def _optimize(
        self,
        scale: bool = None,
        direct_optimization: bool = False,
        *,
        reuse_solver: bool = False,
    ) -> dict[str, ca.DM | ca.MX]:
        """Runs optimizer, uses scaling if needed. Returned values is scaled back."""
        if scale is not None:
            warn("Scale argument is deprecated", FutureWarning, 5)

        if direct_optimization:
            varlist_decision = self.varlist_decision_direct
        else:
            varlist_decision = self.varlist_decision

        self.nlpsol_dict = {
            "x": varlist_decision.get_casadi_variables(),
            "f": self._objective()[0],
            "p": self._nlpsol_p_mx,
        }

        if direct_optimization:
            self.nlpsol_dict["g"] = self.nlpsol_g_direct

        if not (hasattr(self, "solver") and reuse_solver):
            self.solver: ca.Function = ca.nlpsol(
                "solver",
                self.solver_name,
                self.nlpsol_dict,
                self.solver_settings,
            )

        if direct_optimization:
            self.nlpsol_args = {
                "x0": self.guess_direct,
                "lbx": self.lower_bound_direct,
                "ubx": self.upper_bound_direct,
                "lbg": [0] * self.nlpsol_g_direct.shape[0],
                "ubg": [0] * self.nlpsol_g_direct.shape[0],
            }

        else:
            self.nlpsol_args = {
                "x0": self.guess,
                "lbx": self.lower_bound,
                "ubx": self.upper_bound,
            }
        self.nlpsol_args["p"] = self._nlpsol_p_values

        res_solver = self.solver.call(self.nlpsol_args)

        res_solver["x_unscaled"] = res_solver["x"].toarray()
        res_solver["x_all"] = np.asarray(
            varlist_decision.scale_to_original(res_solver["x"])
        )
        if direct_optimization:
            res_solver["x"] = res_solver["x_all"][: len(self.varlist_decision)]
        else:
            res_solver["x"] = res_solver["x_all"]

        res_dict = {}
        for solution, var_name in zip(
            res_solver["x"], list(self.varlist_decision.keys())
        ):
            res_dict[var_name] = float(solution[0])

        if direct_optimization:
            res_dict_all = {}
            for solution, var_name in zip(
                res_solver["x_all"], list(varlist_decision.keys())
            ):
                res_dict_all[var_name] = float(solution[0])
        else:
            res_dict_all = res_dict

        res_solver["x_dict"] = res_dict
        res_solver["x_dict_all"] = res_dict

        return res_solver


    @_consistent_scaling_decorator
    def map_objective(self, plot: bool = True) -> None:
        """Calculate objective function for different values of parameters and plot, if needed.
        Currently support only 3 unfixed decision variables."""
        import matplotlib.pyplot as plt

        decision_variables = self.varlist_decision.get_casadi_variables()
        if decision_variables.shape == (3, 1):
            objective_function = ca.Function(
                "objective",
                [decision_variables],
                [self._objective()[0]],
                ["x"],
                ["f"],
            )

            # Generate steps for each variable between lb and ub
            axis_steps = []
            for lb, ub in zip(self.lower_bound, self.upper_bound):
                axis_steps.append(np.linspace(lb, ub, 6))

            # Create every possible combination of cordinates
            # Code the value of objective function as a color information
            xx = []
            yy = []
            zz = []
            objective = []
            for x in axis_steps[0]:
                for y in axis_steps[1]:
                    for z in axis_steps[2]:
                        xx.append(x)
                        yy.append(y)
                        zz.append(z)
                        objective.append(float(objective_function([x, y, z]).toarray()))

            if plot:
                color = np.array(objective)
                fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
                colormap = plt.get_cmap("tab20b")

                ax.set_xlabel(f"x {decision_variables[0]}")
                ax.set_ylabel(f"y {decision_variables[1]}")
                ax.set_zlabel(f"z {decision_variables[2]}")

                quantile = np.quantile(color, 0.1)
                selected_objective = copy.deepcopy(color)
                selected_objective[selected_objective > quantile] = color.max()
                marker_size = (
                    selected_objective.max() - selected_objective
                ) / selected_objective.max() * 20 + 7

                surf = ax.scatter(xx, yy, zz, c=color, s=marker_size, cmap=colormap)
                fig.colorbar(surf)
                plt.show()
        else:
            raise NotImplementedError

    def generate_multistart_guess(self, num_initials: int):
        if isinstance(self, ParameterEstimationNLE_control):
            hammersley_seeds = np.array(
                list(
                    zip(
                        self.lower_bound[0 : self.num_parameters],
                        self.upper_bound[0 : self.num_parameters],
                    )
                )
            )
        else:
            hammersley_seeds = np.array(list(zip(self.lower_bound, self.upper_bound)))

        list_startpoint = utilities.make_startpoints(
            hammersley_seeds, num_initials, sampling="lhs"
        )
        return list_startpoint

    @_consistent_scaling_decorator
    def optimize_multistart(
        self,
        num_initials: int,
        scale: bool = True,
        max_iterations: int = 20,
    ) -> list[dict[str, ca.DM]]:
        """Runs multiple optimizations with gueses spread between upper and lower bound.
        Helps to find feasible starting point for optimization in a few steps.
        WIP: recalcution of algebraic variables doesn't work.
        """
        list_startpoint = self.generate_multistart_guess(num_initials)
        results_with_f = []

        # Optimizer settings and guess are overwritten for multistart, and then returned back
        initial_guess = copy.deepcopy(self.guess)
        initial_settings = copy.deepcopy(self.solver_settings)

        self.solver_settings = {
            "verbose": False,
            "print_time": False,
            "ipopt": {
                "max_iter": max_iterations,
                "print_level": 0,
            },
        }
        for index, guess in tqdm(enumerate(list_startpoint)):
            if isinstance(self, ParameterEstimationNLE_control):
                for index_guess, current_guess in enumerate(guess):
                    self.guess[index_guess] = current_guess
            else:
                self.guess = guess

            print(f"Optimization number {index} started")
            # for sim in self.list_simulators:
            #     sim.calculate_algebraic_initials(apply_intials=True)

            res = self.optimize(scale)
            print(f"Objective: {res['f']}")
            res_f = float(res["f"])
            results_with_f.append([float("inf") if res_f == 0 else res_f, res])

        results_with_f = sorted(results_with_f, key=lambda res: res[0])
        result_list_sorted = [res[1] for res in results_with_f]

        self.solver_settings = initial_settings
        return result_list_sorted

    @_consistent_scaling_decorator
    def _setup_direct_optimization(self, mode: str):
        """Setup all the variables needed for direct optimization of steady state models.
        mode is needed to differentiate between PE and OED"""
        if mode not in ("PE", "OED"):
            raise NotImplementedError
        self.varlist_algebraic_direct = VariableList()

        lower_bound = self.lower_bound.tolist()
        upper_bound = self.upper_bound.tolist()
        lower_bound_inf = self.lower_bound.tolist()
        upper_bound_inf = self.upper_bound.tolist()
        guess = self.guess.tolist()

        guess_dict = {}
        for var in self.varlist_decision.values():
            guess_dict[var.name] = var.guess
        if mode == "OED":
            for par_name, par_value in zip(
                self.varlist_parameter.keys(), self.parameter_values_unscaled
            ):
                guess_dict[par_name] = par_value

        parameter_symbols = []
        parameter_values = []

        for variable_name, var in self.model.varlist_independent(
            self.list_input_varlist[0]
        ).items():
            if isinstance(var, VariableParameter):
                parameter_symbols.append(var.casadi_var)
                if (var.fixed is True) or (mode == "OED"):
                    parameter_values.append(var.scale_from_original(var.value[0]))
                else:
                    parameter_values.append(var.casadi_var)

        equality_constraints = []
        meas_symbols = []
        jacobian_list = []
        self.varlist_decision_direct = copy.deepcopy(self.varlist_decision)

        for sim_index, input_varlist in enumerate(self.list_input_varlist):
            simulator_i = self.list_simulators[sim_index]

            good_initial_guess = (
                simulator_i.simulate(
                    return_varlist=False, unfixed_variables=guess_dict
                )[0]["x"]
                .toarray()
                .flatten()
            )
            control_symbols = []
            control_values = []
            meas_symbols_sim = []

            for variable_name, var in self.model.varlist_independent(
                input_varlist
            ).items():
                if isinstance(var, VariableControlPiecewiseConstant):
                    raise NotImplementedError
                elif isinstance(var, VariableControl):
                    control_symbols.append(var.casadi_var)
                    if (var.fixed) or (mode == "PE"):
                        control_values.append(var.scale_from_original(var.value[0]))
                    else:
                        control_values.append(var.casadi_var)

            varlist_new_algebraic_i = VariableList()

            for variable_name, var_guess in zip(
                self.model.varlist_algebraic(input_varlist).keys(), good_initial_guess
            ):
                var = input_varlist[variable_name]

                new_var = var._create_copy(f"_sim{sim_index}")

                varlist_new_algebraic_i.add_variable(new_var)
                self.varlist_decision_direct.add_variable(new_var)

                lower_bound.append(new_var.scale_from_original(new_var.lower_bound))
                upper_bound.append(new_var.scale_from_original(new_var.upper_bound))
                lower_bound_inf.append(-ca.inf)
                upper_bound_inf.append(ca.inf)
                guess.append(var_guess)

                if var.name in self.names_of_measurements:
                    meas_symbols_sim.append(new_var.casadi_var)

            meas_symbols.append(ca.hcat(meas_symbols_sim))

            scaled_equations = ca.substitute(
                self.model.equations_algebraic,
                input_varlist.get_casadi_variables(),
                input_varlist.get_scaled_casadi_variables(),
            )

            equations_subs_alg = ca.substitute(
                scaled_equations,
                self.model.varlist_algebraic(input_varlist).get_casadi_variables(),
                varlist_new_algebraic_i.get_casadi_variables(),
            )

            if len(parameter_symbols) == 0:
                equations_subs_par = equations_subs_alg
            else:
                equations_subs_par = ca.substitute(
                    equations_subs_alg,
                    ca.vcat(parameter_symbols),
                    ca.vcat(parameter_values),
                )
            if len(control_symbols) == 0:
                equations_subs_all = equations_subs_par
            else:
                equations_subs_all = ca.substitute(
                    equations_subs_par,
                    ca.vcat(control_symbols),
                    ca.vcat(control_values),
                )

            self.varlist_algebraic_direct.update(**varlist_new_algebraic_i)
            # equality_constraints.append(equations_subs_all.printme(sim_index))
            equality_constraints.append(equations_subs_all)

            if mode == "OED":
                jac_function = simulator_i.function.jacobian()
                independent_variables = ca.substitute(
                    simulator_i._independent_variables,
                    ca.vcat(parameter_symbols),
                    ca.vcat(parameter_values),
                )
                args = {
                    "x": varlist_new_algebraic_i.get_casadi_variables(),
                    "p": independent_variables,
                }
                jac_all = jac_function.call(args)
                jac_parameters = ca.solve(jac_all["jac_rhs_x"], -jac_all["jac_rhs_p"])

                index_selected_parameters = []
                for var_name in self.varlist_parameter.keys():
                    index_selected_parameters.append(
                        simulator_i.mapping_independent_variables[var_name]
                    )

                jac_mx_selected_parameters = jac_parameters.get(
                    False, ca.Slice(), index_selected_parameters
                )
                jacobian_list.append(
                    jac_mx_selected_parameters.get(
                        False, self.index_measurements_in_sim, ca.Slice()
                    )
                )

        self.nlpsol_g_direct = ca.cse(ca.vcat(equality_constraints))
        self.lower_bound_direct = np.array(lower_bound)
        self.lower_bound_direct_inf = np.array(lower_bound_inf)
        self.upper_bound_direct = np.array(upper_bound)
        self.upper_bound_direct_inf = np.array(upper_bound_inf)
        self.guess_direct = np.array(guess)

        if mode == "PE":
            self.simulate_all_direct = ca.vcat(meas_symbols)
        elif mode == "OED":
            jac_mx = ca.vcat(jacobian_list)
            if not self.jacobian_scaled_mx_constant.is_empty():
                jac_mx = ca.vcat([self.jacobian_mx_constant, jac_mx])

            number_of_experiments = len(self._previous_measurements) + len(
                self.list_simulators
            )
            std_scaling = ca.repmat(
                self.array_inverted_scaled_std,
                number_of_experiments,
                len(self.varlist_parameter),
            )

            jac_mx_scaled = jac_mx * std_scaling

            self.jacobian_mx_direct = jac_mx
            self.jacobian_scaled_mx_direct = jac_mx_scaled

    @_consistent_scaling_decorator
    def check_decision_bounds(self, plot: bool = False) -> None:
        """Method is simulating model on upper and lower bounds of decision variables.
        Prints if there were porblems simulation some bounds, meaning that optimizer
        will also have problems, when going near that bounds.
        Only first simulator of optimizers is used"""
        if isinstance(self, ParameterEstimation):
            bound_pairs = []
            for lb, ub in zip(self.lower_bound, self.upper_bound):
                bound_pairs.append([lb, ub])

            list_of_bounds = np.array(np.meshgrid(*bound_pairs)).T.reshape(
                -1, len(bound_pairs)
            )

            simulation = self.list_simulators[0]
            results = []

            for set_of_bounds in list_of_bounds:
                bound_dictionary = {}
                for var, bound in zip(
                    list(self.varlist_decision.values()), set_of_bounds
                ):
                    bound_dictionary[var.name] = bound
                try:
                    result = simulation.simulate(
                        algebraic=True,
                        recalculate_algebraic=True,
                        unfixed_variables=set_of_bounds,
                    )
                    if plot:
                        result._get_varlist_to_plot(True).dataframe.plot(
                            subplots=True, title=str(bound_dictionary)
                        )
                    results.append(result)
                except Exception:
                    print(f"Failed for these bounds: {bound_dictionary}")
        else:
            raise NotImplementedError


class PE_base(Optimizer):
    @property
    def _nlpsol_p_mx(self):
        return self.array_data_mx

    @property
    def _nlpsol_p_values(self):
        return self.array_data

    def _objective_ols(self, direct_optimization: bool):
        """Objective function is a trace(Z.T * Z), where Z is a residual matrix with shape:
        numRows -> amount of supplied experiments, numCol -> amount of variables that have measurements
        If experiments do not supply a measurement for one of the measurements, self.array_data_mask will
        have 0 as the respective element of the martix, otherwise 1"""
        if direct_optimization:
            evaluation = self.simulate_all_direct
        else:
            evaluation = self.simulate_all_mx

        residuals = (
            (evaluation - self.array_data_mx)
            * self.array_data_mask
            * np.sqrt(self.experiments_scale)
        )
        objective = ca.sumsqr(residuals)

        return objective, residuals

    def _objective_wls(self, direct_optimization: bool):
        """Objective function is a trace(Z.T * inv(VarY) * Z), where Z is a same matrix as in _objective_ols
        inv(VarY) is the variance of the respective measurements in Z, and has the same shape.
        Thus, covariance of the measurements is assumed to be zero."""
        if direct_optimization:
            evaluation = self.simulate_all_direct
        else:
            evaluation = self.simulate_all_mx

        residuals = (evaluation - self.array_data_mx) * self.array_data_mask
        scaled_residuals = (
            residuals * self.array_inverted_scaled_std * np.sqrt(self.experiments_scale)
        )
        objective = ca.sumsqr(scaled_residuals)
        return objective, residuals

    def _objective_fair(self, direct_optimization: bool):
        if direct_optimization:
            evaluation = self.simulate_all_direct
        else:
            evaluation = self.simulate_all_mx

        c = 2
        residuals = (evaluation - self.array_data) * self.array_data_mask
        scaled_residuals = (
            residuals * self.array_inverted_scaled_std * np.sqrt(self.experiments_scale)
        )
        res_mod = ca.sqrt(scaled_residuals**2)
        objective = 2 * c**2 * (res_mod / c - ca.log(1 + res_mod / c))
        objective = ca.sum1(objective)
        objective = ca.sum2(objective)
        # objective = ca.sumsqr(scaled_residuals)
        return objective, residuals

    @_consistent_scaling_decorator
    def _unscale_y(self, y):
        return self._unscale_simulator_output(y, self.names_of_measurements)

    @_consistent_scaling_decorator
    def _unscale_simulator_output(self, y, variables_names):
        scaling_constants_measurements_v = []
        scaling_constants_measurements_r = []
        for meas_name in variables_names:
            v, r = self.list_input_varlist[0][meas_name]._get_scaling_constants()
            scaling_constants_measurements_v.append(v)
            scaling_constants_measurements_r.append(r)
        scale_factor_v = np.array([scaling_constants_measurements_v])
        scale_factor_r = np.array([scaling_constants_measurements_r])
        scaled_y = y * scale_factor_v + scale_factor_r
        return scaled_y

    @_consistent_scaling_decorator
    def _unscale_df(self, df_all):
        return self._unscale_simulator_output(df_all, df_all.columns)

    @_consistent_scaling_decorator
    def _unscale_residuals(self, residuals):
        scaling_constants_measurements = []
        for meas_name in self.names_of_measurements:
            scaling_constants_measurements.append(
                self.list_input_varlist[0][meas_name]._get_scaling_constants()[0]
            )
        scale_factor = np.array([scaling_constants_measurements])
        scale_factor = np.tile(scale_factor, (residuals.shape[0], 1))
        scaled_residuals = residuals * scale_factor
        return scaled_residuals

    @_consistent_scaling_decorator
    def _unscale_jacobian(self, jacobian):
        scaled_jacobian = self._unscale_jacobian_parameter_values(jacobian)
        scaled_jacobian = self._unscale_jacobian_measurement_values(scaled_jacobian)
        return scaled_jacobian

    @property
    @_consistent_scaling_decorator
    def array_data_unscaled(self):
        unscaled_list = []
        for i, meas_name in enumerate(self.names_of_measurements):
            scaled_column = self.list_input_varlist[0][meas_name].scale_to_original(
                self.array_data[:, i]
            )
            unscaled_list.append(scaled_column)
        return np.array(unscaled_list).T

    @_consistent_scaling_decorator
    def _unscale_jacobian_parameter_values(self, jacobian):
        scale_parameters = np.tile(
            np.array(self.varlist_decision._get_scaling_constants()[0]),
            (jacobian.shape[0], 1),
        )
        scaled_jacobian = jacobian / scale_parameters
        return scaled_jacobian

    @_consistent_scaling_decorator
    def _unscale_jacobian_measurement_values(self, jacobian):
        scaling_constants_measurements = []
        for meas_name in self.names_of_measurements:
            scaling_constants_measurements.append(
                self.list_input_varlist[0][meas_name]._get_scaling_constants()[0]
            )
        scaling_measurements = np.repeat(
            np.asarray([scaling_constants_measurements]),
            len(self.varlist_decision),
            axis=0,
        ).T

        scaling_all = []

        if isinstance(self.list_simulators[0], Simulator):
            for sim in self.list_simulators:
                len_time_grid = sim.time_grid_relative.shape[0] - 1
                scaling_simulator_i = np.repeat(
                    scaling_measurements, len_time_grid, axis=0
                )
                scaling_all.append(scaling_simulator_i)
            scale_measurements = np.concatenate(scaling_all, axis=0)
        else:
            scale_measurements = np.repeat(
                scaling_measurements, len(self.list_simulators), axis=0
            )

        scaled_jacobian = jacobian * scale_measurements
        return scaled_jacobian

    def setup_regularization(
        self,
        contribution: None | float = None,
        reference_parameters: None | np.ndarray = None,
    ):
        if contribution is None:
            self.regularization_contribution = 0
        else:
            self.regularization_contribution = contribution

        if reference_parameters is None:
            self.reference_parameters = np.zeros((len(self.varlist_decision), 1))
        else:
            if reference_parameters.shape[0] != len(self.varlist_decision):
                raise ValueError("Shape of supplied reference_parameters is incorrect")
            else:
                self.reference_parameters = reference_parameters

    def _objective_tikhonov(self, direct_optimization: bool):
        objective, residuals = self._objective_wls(
            direct_optimization=direct_optimization
        )

        penalty = ca.sqrt(
            ca.sumsqr(
                self.varlist_decision.get_casadi_variables() - self.reference_parameters
            )
        )
        regularization_part = 0.5 * (self.regularization_contribution**2) * penalty

        objective = objective + regularization_part

        return objective, residuals

    def optimize(
        self,
        scale=None,
        objective_function="wls",
        direct_optimization=False,
        *,
        reuse_solver=False,
    ):
        if objective_function == "wls":
            self._objective = partial(
                self._objective_wls, direct_optimization=direct_optimization
            )
        elif objective_function == "ols":
            self._objective = partial(
                self._objective_ols, direct_optimization=direct_optimization
            )
        elif objective_function == "fair":
            self._objective = partial(
                self._objective_fair, direct_optimization=direct_optimization
            )
        elif objective_function == "tikh":
            self._objective = partial(
                self._objective_tikhonov, direct_optimization=direct_optimization
            )
        else:
            raise NotImplementedError(
                f"Objective function '{objective_function}' is not supported"
            )

        return self._optimize(
            scale, direct_optimization=direct_optimization, reuse_solver=reuse_solver
        )

    def _setup_experiments_scale(self, scale_experiments):
        if isinstance(self, ParameterEstimation):
            if scale_experiments:
                experiments_scale: int | np.ndarray = self.experiments_weights
            else:
                experiments_scale = 1

            # This attribute is used while calculating Objective, and is either 1 or self.experiments_weights
            # It's used to make some experiments as valuable as others, even if they have less experimental points
            # So if you supply 2 experiments one with 10 and another with 20 time_stamps, effect of each experimental
            # point of second experiments on objective function is decreased by 2
            self.experiments_scale: int | np.ndarray = experiments_scale
        else:
            self.experiments_scale = 1

    @_consistent_scaling_decorator
    def calculate_objective_and_residual(
        self,
        parameters: dict[str, float],
        objective_function: str = "ols",
        experiment_weigts: bool = False,
    ) -> dict[str, float | np.ndarray]:
        if objective_function == "ols":
            obj_f = self._objective_ols(direct_optimization=False)
        elif objective_function == "wls":
            obj_f = self._objective_wls(direct_optimization=False)

        decision_variables = self.varlist_decision.get_casadi_variables()

        casadi_function = ca.Function(
            "objective",
            [decision_variables, self.array_data_mx],
            [obj_f[0], obj_f[1], self.simulate_all_mx, self._simulate_all_mx],
            ["x", "data_arrays"],
            ["f", "residuals", "y", "y_all"],
        )

        selected_parameters = self.variables_dict_to_list(parameters)
        res = casadi_function(x=selected_parameters, data_arrays=self.array_data)

        if isinstance(self, ParameterEstimationNLE):
            algebraic_names = list(
                self.model.varlist_algebraic(self.list_input_varlist[0]).keys()
            )
            df_all = pd.DataFrame(res["y_all"], columns=algebraic_names)
        else:
            all_names = list(
                self.model.varlist_state(self.list_input_varlist[0]).keys()
            )
            if self._use_algebraic_variables:
                all_names.extend(
                    self.model.varlist_algebraic(self.list_input_varlist[0]).keys()
                )

            index_from = 0
            index_till = 0
            list_df = []
            for sim in self.list_simulators:
                time = sim.time_grid_relative[1:]
                index_till += len(time)
                df_one_simulator = pd.DataFrame(
                    res["y_all"][index_from:index_till, :],
                    columns=all_names,
                    index=time,
                )
                index_from += len(time)
                list_df.append(df_one_simulator)

            df_all = pd.concat(
                list_df, keys=range(len(self.list_simulators)), names=["sim", "time"]
            )

        result_np = {
            "f": float(res["f"]),
            "residuals": res["residuals"].toarray(),
            "y": res["y"].toarray(),
            "df_all": df_all,
        }

        return result_np

    def _setup_varlist_decision(self):
        for var in self.model.varlist_independent(self.list_input_varlist[0]).values():
            if isinstance(var, VariableParameter):
                if var.fixed is False:
                    self.varlist_decision.add_variable(var)
        self.setup_regularization()

    @_consistent_scaling_decorator
    def generate_simulate_all_functions(self) -> None:
        """Combines simulate_fast() functions from simulator, and creates MX structure, that is used
        further in objective_function calculation"""
        if isinstance(self.list_simulators[0], Simulator):
            res_dict_name = "xf"
        elif isinstance(self.list_simulators[0], SimulatorNLE):
            if self.list_simulators[0]._solver_name == "ipopt":
                print(
                    "\nSimulators of PE optimizer use IPOPT nlpsol. Results can be incosistent\n"
                )
            res_dict_name = "x"

        list_simulation_T = []

        for simulator in self.list_simulators:
            if isinstance(self.list_simulators[0], Simulator):
                simulator.calculate_algebraic_initials(apply_intials=True)

            res_simulation = simulator.simulate_fast()

            if getattr(self, "_use_algebraic_variables", False):
                data = ca.vcat([res_simulation[res_dict_name], res_simulation["zf"]])
                list_simulation_T.append(data.T)
            else:
                list_simulation_T.append(res_simulation[res_dict_name].T)

        free_variables = self.varlist_decision.get_casadi_variables()
        all_selected_measurements = ca.vcat(list_simulation_T).get(
            False, ca.Slice(), self.index_measurements_in_sim
        )
        self.simulate_all_function = ca.Function(
            "sim_all", [free_variables], [all_selected_measurements]
        )
        self.simulate_all_mx = self.simulate_all_function(free_variables)

        _simulate_all_function = ca.Function(
            "sim_all", [free_variables], [ca.vcat(list_simulation_T)]
        )
        self._simulate_all_mx = _simulate_all_function(free_variables)

    @cached_property
    def _jacobian_function(self):
        decision_variables = self.varlist_decision.get_casadi_variables()
        jac_meas_mx = ca.jacobian(
            self.simulate_all_mx[:, list(range(len(self.names_of_measurements)))],
            decision_variables,
        )
        jac_meas_function = ca.Function("jac_meas", [decision_variables], [jac_meas_mx])
        return jac_meas_function

    @cached_property
    def _jacobian_scaler(self):
        flattened_std = self.array_inverted_std.flatten(order="F")
        return np.tile(flattened_std, (len(self.varlist_decision), 1)).T

    @_consistent_scaling_decorator
    def calculate_sensitivity_and_fim_fast(
        self, parameters: dict[str, float]
    ) -> dict[str, np.ndarray]:
        all_parameter_values = self.variables_dict_to_list(parameters)
        jac_all_dm = self._unscale_jacobian(
            self._jacobian_function(all_parameter_values)
        )
        jac_all_scaled = jac_all_dm * self._jacobian_scaler

        fim_matrix_scaled = jac_all_scaled.T @ jac_all_scaled
        parameter_covariance_matrix = np.linalg.inv(fim_matrix_scaled)  # type: ignore

        return jac_all_dm, jac_all_scaled, parameter_covariance_matrix

    @_consistent_scaling_decorator
    def calculate_jacobian_yao_fast(
        self, parameters: dict[str, float]
    ) -> dict[str, np.ndarray]:
        all_parameter_values = self.variables_dict_to_list(parameters)
        jac_all_dm = self._unscale_jacobian(
            self._jacobian_function(all_parameter_values)
        )

        decision_variables = self.varlist_decision.get_casadi_variables()
        res_simulation = ca.Function(
            "sim", [decision_variables], [self.simulate_all_mx]
        )(all_parameter_values)
        jacobian_index = [0, self.simulate_all_mx.shape[0]]
        jacs = []

        for index_measurement, meas_name in enumerate(self.names_of_measurements):
            scale_factor = self.list_input_varlist[0][
                meas_name
            ]._get_scaling_constants()
            jacobian_slice = ca.Slice(jacobian_index[0], jacobian_index[1])
            jac_meas_dm = jac_all_dm[jacobian_slice, :]
            jacobian_index[0] += self.simulate_all_mx.shape[0]
            jacobian_index[1] += self.simulate_all_mx.shape[0]

            jac_meas_selected_dm = (
                jac_meas_dm * self.array_data_mask[:, index_measurement]
            )

            jac_meas_selected_yao_dm = jac_meas_selected_dm * (
                1
                / (
                    res_simulation[:, index_measurement] * scale_factor[0]
                    + scale_factor[1]
                )
            )
            jacs.append(jac_meas_selected_yao_dm)

        jac_array_yao = np.concatenate(jacs)
        return jac_array_yao

    @_consistent_scaling_decorator
    def calculate_sensitivity_and_fim(
        self, parameters: dict[str, float], parameter_names: list[str] | None = None
    ) -> dict[str, np.ndarray]:
        """Calculate jacobian, scaled_jacobian, parameter_covariance_matrix only for parameters,
        which names are listed in paramtere_names, if None, for all unfixed paramters.
        Jacobian is calculated only for "measured" algebraic variables, the ones which
        had value, when ParameterEstimationNLE was initialized.

        parameters dictionary holds paramter values for all unfixed parameters, example:

        parameters = {"theta1": 1, "theta2": 2}
        parameter_names = ["theta1"]

        Jacobian of measured variables will be calculated for theta1=1 and theta2=2, but reported
        jacobian will only contain parameter theta2.

        Jacobian dimensions: dY/dp [NumOfMeasurements x NumOfParameters]
        """
        decision_variables = self.varlist_decision.get_casadi_variables()
        if parameter_names is not None:
            list_selected_parameters_index = []
            for par_index, par_name in enumerate(self.varlist_decision.keys()):
                if par_name in parameter_names:
                    list_selected_parameters_index.append(par_index)

        all_parameter_values = self.variables_dict_to_list(parameters)

        residuals = self.calculate_objective_and_residual(
            parameters, objective_function="ols"
        )["residuals"]

        scaled_residuals = self._unscale_residuals(residuals)

        measurement_variance_estimate = (
            np.diag(scaled_residuals.T @ scaled_residuals) / self.dof
        )
        print("OLS std: ", np.sqrt(measurement_variance_estimate))

        estimated_inverted_std = copy.deepcopy(self.array_inverted_std)

        if isinstance(measurement_variance_estimate, float):
            measurement_variance_estimate = [measurement_variance_estimate]

        # Avoid division by 0 later
        measurement_variance_estimate[measurement_variance_estimate == 0] = 1e-24

        for index_meas, meas_std in enumerate(measurement_variance_estimate):
            estimated_inverted_std[:, index_meas] = 1 / np.sqrt(meas_std)

        jacobian = {}
        jacobian_scaled = {}
        jacobian_scaled_estimated = {}
        jacobian_yao = {}
        res_simulation = ca.Function(
            "sim", [decision_variables], [self.simulate_all_mx]
        )(all_parameter_values)

        jac_meas_mx = ca.jacobian(
            self.simulate_all_mx[:, list(range(len(self.names_of_measurements)))],
            decision_variables,
        )
        jac_meas_function = ca.Function("jac_meas", [decision_variables], [jac_meas_mx])
        jac_all_dm = self._unscale_jacobian(jac_meas_function(all_parameter_values))
        jacobian_index = [0, self.simulate_all_mx.shape[0]]

        for index_measurement, meas_name in enumerate(self.names_of_measurements):
            scale_factor = self.list_input_varlist[0][
                meas_name
            ]._get_scaling_constants()
            jacobian_slice = ca.Slice(jacobian_index[0], jacobian_index[1])
            jac_meas_dm = jac_all_dm[jacobian_slice, :]
            jacobian_index[0] += self.simulate_all_mx.shape[0]
            jacobian_index[1] += self.simulate_all_mx.shape[0]

            jac_meas_selected_dm = (
                jac_meas_dm * self.array_data_mask[:, index_measurement]
            )
            jac_meas_selected_scaled_dm = (
                jac_meas_selected_dm * self.array_inverted_std[:, index_measurement]
            )
            jac_meas_selected_scaled_estimated_dm = (
                jac_meas_selected_dm * estimated_inverted_std[:, index_measurement]
            )
            jac_meas_selected_yao_dm = jac_meas_selected_dm * (
                1
                / (
                    res_simulation[:, index_measurement] * scale_factor[0]
                    + scale_factor[1]
                )
            )
            if parameter_names is None:
                jacobian[meas_name] = jac_meas_selected_dm
                jacobian_scaled[meas_name] = jac_meas_selected_scaled_dm
                jacobian_scaled_estimated[meas_name] = (
                    jac_meas_selected_scaled_estimated_dm
                )
                jacobian_yao[meas_name] = jac_meas_selected_yao_dm
            else:
                jacobian[meas_name] = jac_meas_selected_dm[
                    :, list_selected_parameters_index
                ]
                jacobian_scaled[meas_name] = jac_meas_selected_scaled_dm[
                    :, list_selected_parameters_index
                ]
                jacobian_scaled_estimated[meas_name] = (
                    jac_meas_selected_scaled_estimated_dm[
                        :, list_selected_parameters_index
                    ]
                )
                jacobian_yao[meas_name] = jac_meas_selected_yao_dm[
                    :, list_selected_parameters_index
                ]

        jac_array = np.concatenate(list(jacobian.values()))
        jac_array_scaled = np.concatenate(list(jacobian_scaled.values()))
        jac_array_scaled_estimated = np.concatenate(
            list(jacobian_scaled_estimated.values())
        )
        jac_array_yao = np.concatenate(list(jacobian_yao.values()))

        # Generate jacobian and hessian on obj function
        obj_func = ca.substitute(
            self._objective_wls(direct_optimization=False)[0],
            self.array_data_mx,
            self.array_data,
        )
        jac_objective = ca.Function(
            "jf",
            [decision_variables],
            [ca.jacobian(obj_func, decision_variables)],
        )(all_parameter_values)
        # Should be twice as big as fim_matrix_scaled
        try:
            hessian_objective_wls = ca.Function(
                "jf",
                [decision_variables],
                [ca.hessian(obj_func, decision_variables)[0]],
            )
            hessian_objective_wls = hessian_objective_wls(all_parameter_values)
        except RuntimeError:
            print("Failed to calculate hessian")
            hessian_objective_wls = None

        try:
            obj_func_tikhonov = ca.substitute(
                self._objective_tikhonov(direct_optimization=False)[0],
                self.array_data_mx,
                self.array_data,
            )
            hessian_objective_tikhonov = ca.Function(
                "jf",
                [decision_variables],
                [ca.hessian(obj_func_tikhonov)],
            )
            hessian_objective_tikhonov = hessian_objective_tikhonov(
                all_parameter_values
            )
        except RuntimeError:
            print("Failed to calculate hessian")
            hessian_objective_tikhonov = None

        if parameter_names is not None:
            jac_objective = jac_objective[:, list_selected_parameters_index]
            if hessian_objective_wls is not None:
                hessian_objective_wls = hessian_objective_wls[
                    list_selected_parameters_index, list_selected_parameters_index
                ]
            if hessian_objective_tikhonov is not None:
                hessian_objective_tikhonov = hessian_objective_tikhonov[
                    list_selected_parameters_index, list_selected_parameters_index
                ]

        fim_matrix = jac_array.T @ jac_array
        fim_matrix_scaled = jac_array_scaled.T @ jac_array_scaled
        parameter_covariance_matrix = np.linalg.inv(fim_matrix_scaled)  # type: ignore

        result = {}
        result["jac_full"] = jac_array
        result["jac_sorted"] = jacobian
        result["jac_scaled_full"] = jac_array_scaled_estimated
        result["jac_scaled_full_theory"] = jac_array_scaled
        result["jac_scaled_sorted"] = jacobian_scaled
        result["jac_yao_full"] = jac_array_yao
        result["jac_yao_sorted"] = jacobian_yao
        result["fim"] = fim_matrix
        result["fim_scaled"] = fim_matrix_scaled
        result["cov_par"] = parameter_covariance_matrix
        result["jac_wls"] = jac_objective
        result["hess_wls"] = hessian_objective_wls
        result["hess_tikh"] = hessian_objective_tikhonov
        result["s2"] = measurement_variance_estimate

        return result

    @_consistent_scaling_decorator
    def parameter_identifiability_chu2012(
        self,
        parameters: dict[str, float],
        unfixed_params: list[str],
        parameters_identifiable: list[str] | None = None,
        parameters_not_identifiable: list[str] | None = None,
    ):
        if parameters_identifiable is None:
            parameters_identifiable = []

        if parameters_not_identifiable is None:
            parameters_not_identifiable = []

        sorted_unfixed_params = []
        for par_name in self.varlist_decision.keys():
            if par_name in unfixed_params:
                sorted_unfixed_params.append(par_name)

        parameters_index = list(range(len(sorted_unfixed_params)))

        S = self.calculate_sensitivity_and_fim_fast(parameters)[1].toarray()
        S = S * np.array(self.variables_dict_to_list(parameters, scaling=False))

        info = []
        best_set = None
        max_det = 0

        for subset_size in range(1, len(unfixed_params) + 1):
            for subset_index in itertools.combinations(parameters_index, subset_size):
                S_selected = S[:, subset_index]
                FIM_selected = S_selected.T @ S_selected
                if subset_size == 1:
                    max_det_i = FIM_selected.item(0)
                else:
                    max_det_i = np.linalg.det(FIM_selected)
                if max_det_i > max_det:
                    max_det = max_det_i
                    best_set = subset_index
                subset_names = np.array(sorted_unfixed_params)[list(subset_index)]
                info_i = [subset_size, subset_names, max_det_i]
                info.append(info_i)

        df = pd.DataFrame(info, columns=["subset_size", "subset_names", "det"])
        parameters_identifiable = np.array(sorted_unfixed_params)[list(best_set)]

        parameters_not_identifiable = list(
            set(sorted_unfixed_params) - set(parameters_identifiable)
        )

        parameters_identifiable_sorted = []
        parameters_not_identifiable_sorted = []
        for par_name in sorted_unfixed_params:
            if par_name in parameters_identifiable:
                parameters_identifiable_sorted.append(par_name)
            else:
                parameters_not_identifiable_sorted.append(par_name)

        print(f"Estimable parameters: {parameters_identifiable_sorted}")
        print(f"Non identifiable parameters: {parameters_not_identifiable_sorted}")

        result = {}
        result["estimable"] = parameters_identifiable_sorted
        result["fixed"] = parameters_not_identifiable_sorted

        return result

    @_consistent_scaling_decorator
    def parameter_identifiability_brun2001(
        self,
        parameters: dict[str, float],
        unfixed_params: list[str],
        parameters_identifiable: list[str] | None = None,
        parameters_not_identifiable: list[str] | None = None,
        eigenvalue_threshold: float = 10e-4,
    ):
        if (
            self._created_with_options["variable_scaling"]
            or get_options()["variable_scaling"]
        ):
            raise NotImplementedError(
                "Brun identification analysis with scaling is dependent on operating system"
            )
        if parameters_identifiable is None:
            parameters_identifiable = []

        if parameters_not_identifiable is None:
            parameters_not_identifiable = []

        sorted_unfixed_params = []
        for par_name in self.varlist_decision.keys():
            if par_name in unfixed_params:
                sorted_unfixed_params.append(par_name)

        parameters_index = list(range(len(sorted_unfixed_params)))

        S = self.calculate_sensitivity_and_fim_fast(parameters)[1].toarray()
        S = S * np.array(self.variables_dict_to_list(parameters, scaling=False))
        S_norm = S / np.linalg.norm(S, axis=0)

        beta_msqr = np.sqrt(np.sum(S**2, axis=0) / S.shape[0])
        parameters_ranked = list(np.array(sorted_unfixed_params)[beta_msqr.argsort()])

        info = []
        identifiable_subset_size = None

        for subset_size in range(2, len(unfixed_params) + 1):
            min_gamma = 20
            for subset_index in itertools.combinations(parameters_index, subset_size):
                S_norm_subset = S_norm[:, subset_index]
                FIM = S_norm_subset.T @ S_norm_subset
                try:
                    gamma_k = 1 / np.sqrt(eigsorted(FIM)[0][-1])
                    rho_k = np.linalg.det(FIM) ** (1 / (2 * S.shape[1]))
                except np.linalg.LinAlgError:
                    gamma_k = np.nan
                    rho_k = np.nan
                subset_names = np.array(sorted_unfixed_params)[list(subset_index)]
                if min_gamma > gamma_k:
                    min_gamma = gamma_k
                info_i = [subset_size, subset_names, rho_k, gamma_k]
                info.append(info_i)

            if min_gamma > 10:
                break
            else:
                identifiable_subset_size = subset_size

        df = pd.DataFrame(info, columns=["subset_size", "subset_names", "rho", "gamma"])
        df_identifiable = df.groupby("subset_size").get_group(identifiable_subset_size)
        parameters_identifiable = list(
            df.loc[df_identifiable.idxmin(numeric_only=True).gamma].subset_names
        )

        parameters_not_identifiable = list(
            set(parameters_ranked) - set(parameters_identifiable)
        )

        parameters_identifiable_sorted = []
        parameters_not_identifiable_sorted = []
        for par_name in sorted_unfixed_params:
            if par_name in parameters_identifiable:
                parameters_identifiable_sorted.append(par_name)
            else:
                parameters_not_identifiable_sorted.append(par_name)

        print(f"Ranked parameters: {parameters_ranked}")
        print(f"Estimable parameters: {parameters_identifiable_sorted}")
        print(f"Non identifiable parameters: {parameters_not_identifiable_sorted}")

        result = {}
        result["ranked"] = parameters_ranked
        result["estimable"] = parameters_identifiable_sorted
        result["fixed"] = parameters_not_identifiable_sorted

        return result

    @_consistent_scaling_decorator
    def parameter_identifiability_lopez2013(
        self,
        parameters: dict[str, float],
        unfixed_params: list[str],
        parameters_identifiable: list[str] | None = None,
        parameters_not_identifiable: list[str] | None = None,
        eigenvalue_threshold: float = 10e-4,
    ):
        if parameters_identifiable is None:
            parameters_identifiable = []

        if parameters_not_identifiable is None:
            parameters_not_identifiable = []

        sorted_unfixed_params = []
        for par_name in self.varlist_decision.keys():
            if par_name in unfixed_params:
                sorted_unfixed_params.append(par_name)

        S = self.calculate_sensitivity_and_fim_fast(parameters)[1].toarray()
        S = S * np.array(self.variables_dict_to_list(parameters, scaling=False))

        # S = S * np.array(self.variables_dict_to_list(parameters))

        svd = np.linalg.svd(S, full_matrices=True)
        Q, R, P = linalg.qr(S, pivoting=True)

        conditional_number = (svd[1][0] / svd[1]) > 1000
        if (~conditional_number).all():
            num_identifiable = svd[1].shape[0]
        else:
            num_identifiable = list(conditional_number).index(True)
        parameters_ranked = np.array(sorted_unfixed_params)[P]
        parameters_identifiable = parameters_ranked[:num_identifiable]
        parameters_not_identifiable = parameters_ranked[num_identifiable:]

        parameters_identifiable_sorted = []
        parameters_not_identifiable_sorted = []
        for par_name in sorted_unfixed_params:
            if par_name in parameters_identifiable:
                parameters_identifiable_sorted.append(par_name)
            else:
                parameters_not_identifiable_sorted.append(par_name)

        print(f"Ranked parameters: {list(parameters_ranked)}")
        print(f"Estimable parameters: {parameters_identifiable_sorted}")
        print(f"Non identifiable parameters: {parameters_not_identifiable_sorted}")

        result = {}
        result["ranked"] = list(parameters_ranked)
        result["estimable"] = parameters_identifiable_sorted
        result["fixed"] = parameters_not_identifiable_sorted

        return result

    @_consistent_scaling_decorator
    def parameter_identifiability_quaiser2009(
        self,
        parameters: dict[str, float],
        unfixed_params: list[str],
        parameters_identifiable: list[str] | None = None,
        parameters_not_identifiable: list[str] | None = None,
        eigenvalue_threshold: float = 10e-4,
    ):
        """Do parameter ranking based on Quasier 2009, however use scaled sensitivity as in Yao 2003.
        Threshold is taken from Quasier 2009.
        Return ranked parameters in descending order, and divide them in identifiable and not"""

        if parameters_identifiable is None:
            parameters_identifiable = []

        if parameters_not_identifiable is None:
            parameters_not_identifiable = []

        sorted_unfixed_params = []
        for par_name in self.varlist_decision.keys():
            if par_name in unfixed_params:
                sorted_unfixed_params.append(par_name)

        S = self.calculate_sensitivity_and_fim_fast(parameters)[1].toarray()
        S = S * np.array(self.variables_dict_to_list(parameters, scaling=False))
        fim_matrix = S.T @ S

        for i in range(fim_matrix.shape[0]):
            vals, vecs = eigsorted(fim_matrix)

            index_max = np.argmax(np.abs(vecs[:, -1]))
            current_parameter_name = sorted_unfixed_params.pop(index_max)
            if np.abs(vals[-1]) > eigenvalue_threshold:
                parameters_identifiable.insert(0, current_parameter_name)
            else:
                parameters_not_identifiable.insert(0, current_parameter_name)
            fim_matrix = np.delete(fim_matrix, index_max, axis=0)
            fim_matrix = np.delete(fim_matrix, index_max, axis=1)

        parameters_ranked = []
        for parameter_name in parameters_identifiable + parameters_not_identifiable:
            parameters_ranked.append(parameter_name)

        print(f"Ranked parameters: {parameters_ranked}")
        print(f"Estimable parameters: {parameters_identifiable}")
        print(f"Non identifiable parameters: {parameters_not_identifiable}")

        result = {}
        result["ranked"] = parameters_ranked
        result["estimable"] = parameters_identifiable
        result["fixed"] = parameters_not_identifiable

        return result

    @_consistent_scaling_decorator
    def parameter_identifiability_yao2003(
        self,
        parameters: dict[str, float],
        unfixed_params: list[str],
        threshold: float = 4e-2,
    ):
        """Do parameter ranking based on Yao 2003. Cut-off value taken from Yao 2003.
        Return ranked parameters in descending order, and divide them in identifiable and not"""
        parameter_values_all: list[float] = []
        selected_parameters: list[float] = []
        unranked_parameters: list[str] = []

        sorted_unfixed_params = []
        for par_name in self.varlist_decision.keys():
            if par_name in unfixed_params:
                sorted_unfixed_params.append(par_name)

        for var_name in parameters.keys():
            if var_name in self.varlist_decision.keys():
                parameter_values_all.append(parameters[var_name])

        for var_name in parameters.keys():
            if var_name in unfixed_params:
                selected_parameters.append(parameters[var_name])
                unranked_parameters.append(var_name)

        results_sensitivity = self.calculate_jacobian_yao_fast(parameters)
        jacobian_yao = results_sensitivity * np.array(
            self.variables_dict_to_list(parameters, scaling=False)
        )

        XK = np.zeros(jacobian_yao.shape)

        parameters_identifiable = []
        parameters_not_identifiable = []

        for i in range(len(unranked_parameters)):
            if i == 0:
                eucnorm = np.linalg.norm(jacobian_yao, axis=0)
            else:
                eucnorm = np.linalg.norm(R, axis=0)

            index_most_identifiable_par = np.argsort(eucnorm)[-1]
            most_identifiable_parameter = unranked_parameters[
                index_most_identifiable_par
            ]

            if max(eucnorm) < threshold:
                parameters_not_identifiable.append(most_identifiable_parameter)
                break
            else:
                parameters_identifiable.append(most_identifiable_parameter)

            if i == 0:
                XK = jacobian_yao[:, index_most_identifiable_par].reshape(
                    (jacobian_yao.shape[0], 1)
                )
            else:
                XK = np.append(
                    XK,
                    jacobian_yao[:, index_most_identifiable_par].reshape(
                        (jacobian_yao.shape[0], 1)
                    ),
                    axis=1,
                )

            Z_hat = XK.dot(np.linalg.inv(XK.T.dot(XK))).dot(XK.T).dot(jacobian_yao)
            R = jacobian_yao - Z_hat

        parameters_identifiable_sorted = []
        parameters_not_identifiable_sorted = []
        for par_name in sorted_unfixed_params:
            if par_name in parameters_identifiable:
                parameters_identifiable_sorted.append(par_name)
            else:
                parameters_not_identifiable_sorted.append(par_name)

        print(f"Estimable parameters: {parameters_identifiable_sorted}")
        print(f"Non identifiable parameters: {parameters_not_identifiable_sorted}")

        result = {}
        result["estimable"] = parameters_identifiable
        result["fixed"] = parameters_not_identifiable

        return result

    @_consistent_scaling_decorator
    def parameter_analysis(self, parameters: dict[str, float], plot=True):
        import scipy.stats

        num_par = len(self.varlist_decision)

        selected_parameters = self.variables_dict_to_list(parameters, scaling=False)
        result_sens = self.calculate_sensitivity_and_fim(parameters)

        parameter_covariance_matrix = result_sens["cov_par"]

        parameter_variance = np.diag(parameter_covariance_matrix)
        parameter_std = np.sqrt(parameter_variance).flatten()

        students_t_dist_95 = scipy.stats.t.ppf(0.975, self.dof)
        marginal_conf_interval_95 = (parameter_std * students_t_dist_95).T

        print(parameter_std)
        for par, var_value in zip(selected_parameters, marginal_conf_interval_95):
            print(f"{par} +- {var_value} |  ({var_value * 100 / par}%)")

        if plot:
            import matplotlib.pyplot as plt
            from matplotlib.axes import Axes
            from matplotlib.patches import Ellipse

            fig, axes = plt.subplots(
                ncols=num_par - 1, nrows=num_par - 1, layout="constrained"
            )
            if isinstance(axes, Axes) == 1:
                axes = [axes]

            fisher_f_dist_95 = scipy.stats.f.ppf(0.95, num_par, self.dof)

            comb = list(combinations(range(num_par), 2))
            par_names = list(self.varlist_decision.keys())

            title = ""
            for par_value, var_variance_i, name in zip(
                selected_parameters, marginal_conf_interval_95, par_names
            ):
                title = (
                    title
                    + f"{name}: {round(par_value, 5)} ± {round(var_variance_i, 5)} |  ({round((var_variance_i / par_value) * 100, 1)}%)\n"
                )

            fig.suptitle(title)

            for i in comb:
                if len(axes) == 1:
                    ax = axes[0]
                else:
                    ax = axes[i[1] - 1, i[0]]
                index_subarray = np.ix_(i, i)
                parameters_i = []
                parameters_i.append(selected_parameters[i[0]])
                parameters_i.append(selected_parameters[i[1]])

                marginal_conf_interval_95_i = []
                marginal_conf_interval_95_i.append(marginal_conf_interval_95[i[0]])
                marginal_conf_interval_95_i.append(marginal_conf_interval_95[i[1]])
                cov_m = parameter_covariance_matrix[index_subarray]
                vals, vecs = eigsorted(cov_m)
                theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

                # Width and height are "full" widths, not radius
                for fisher in [fisher_f_dist_95]:  # , fisher_f_dist_99]:
                    width, height = 2 * np.sqrt(num_par * fisher * vals)
                    ellip = Ellipse(
                        xy=[parameters_i[0], parameters_i[1]],
                        width=width,
                        height=height,
                        angle=theta,
                        alpha=0.3,
                        lw=2,
                        linestyle="-",
                        color="red",
                    )

                    ax.add_artist(ellip)
                ax.relim()
                ax.autoscale()

                if i[0] == 0:
                    ax.set_ylabel(f"{par_names[i[1]]}")

                if i[1] == len(par_names) - 1:
                    ax.set_xlabel(f"{par_names[i[0]]}")

                # ax.axvline(parameters_i[0] - marginal_conf_interval_95_i[0])
                # ax.axvline(parameters_i[0] + marginal_conf_interval_95_i[0])
                # ax.axhline(parameters_i[1] - marginal_conf_interval_95_i[1])
                # ax.axhline(parameters_i[1] + marginal_conf_interval_95_i[1])
        return marginal_conf_interval_95

    @property
    def dof(self):
        # Eq 7-13-22 Bard 1974
        return self.array_data.shape[0] - (
            len(self.varlist_decision) / len(self.names_of_measurements)
        )


class ParameterEstimation(PE_base):
    def __init__(
        self,
        model: Model,
        variable_list: list[VariableList],
        simulator_name: str = "idas",
        simulator_settings: dict | None = None,
        *,
        use_idas_constraints: bool = None,
        recalculate_algebraic: bool = False,
    ):
        super().__init__(
            model,
            variable_list,
            simulator_name,
            simulator_settings,
        )

        self._setup_algebraic_flag()

        self._objective: Callable[[], tuple[ca.MX | ca.DM, ca.MX | ca.DM]] = (
            self._objective_ols
        )

        self._setup_simulator(
            use_idas_constraints=use_idas_constraints,
            recalculate_algebraic=recalculate_algebraic,
        )

        self.logger.debug(
            "Created Optimizer object: \n Data Shape {} \n Desicion Variables {}".format(
                self.array_data.shape,
                self.varlist_decision.get_variable_name(),  # type: ignore
            )
        )
        self._setup_initialization()

        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {"max_iter": 300},
        }

        self._setup_experiments_scale(False)

    def _setup_algebraic_flag(self):
        self._use_algebraic_variables = False
        for varlist_input in self.list_input_varlist:
            if len(varlist_input.get_algebraic()) == 0:
                continue
            if varlist_input.get_algebraic().dataframe[1:].empty is False:
                self._use_algebraic_variables = True

    @_consistent_scaling_decorator
    def _setup_simulator(
        self,
        *,
        use_idas_constraints: bool,
        recalculate_algebraic: bool,
    ) -> None:
        # It's not checked if all supplied varlist have same states etc.
        self._setup_varlist_decision()

        # Lists used to calculate experiments_weights
        list_timegrid_length = []
        size_simulation_output = []
        list_simulators = []
        experimental_data = []
        experimental_data_mask = []

        list_inverted_variances = []
        list_inverted_scaled_variances = []

        for simulator_index, varlist_input in enumerate(self.list_input_varlist):
            # Create a time_grid, that "stops" at every experimental data, for every state variable
            ordered_varlist_input = self.model.varlist(varlist_input)
            if not varlist_input.get_common_origin(
                strict=True, variable_type=VariableState
            ):
                raise (
                    ValueError(
                        f"Not all State Variables in one experiment have same time0, so simulations cannot be initialized:\n{varlist_input}"
                    )
                )
            data_frame = pd.DataFrame()

            for var in ordered_varlist_input.values():
                if isinstance(var, VariableState) or (
                    isinstance(var, VariableAlgebraic) and self._use_algebraic_variables
                ):
                    data_frame = data_frame.join(
                        var.scale_from_original(var.dataframe), how="outer"
                    )

                elif isinstance(var, VariableControl):
                    var.fixed = True
                    if isinstance(var, VariableControlPiecewiseConstant):
                        var.fixed = True
                        data_frame = data_frame.join(
                            var.scale_from_original(var.dataframe), how="outer"
                        )
                        # Column should be dropped, because it's needed only for unique timestamp
                        data_frame.drop(columns=var.name, inplace=True)

            if self._use_algebraic_variables:
                # TODO in further stepps I always assume that state variables are there, and algebraic are added
                # Without this sorting steps, self.data_array is not correctly sorted. However, the logic has to be fixed
                new_order = list(ordered_varlist_input.get_state().keys()) + list(
                    ordered_varlist_input.get_algebraic().keys()
                )
                data_frame = data_frame[new_order]

            time_grid_unique = (
                (data_frame.index - data_frame.index[0]).total_seconds().tolist()
            )

            list_timegrid_length.append(float(len(time_grid_unique)))

            simulator_settings = self.simulator_settings

            simulator = Simulator(
                self.model,
                np.array(time_grid_unique),
                ordered_varlist_input,
                self.simulator_name,
                simulator_settings,
                use_idas_constraints=use_idas_constraints,
                recalculate_algebraic=recalculate_algebraic,
            )

            list_simulators.append(simulator)

            # Generate an array (experiment_data_varlist) with Experimental data with the same dimensions as simulation results.
            new_experiment_data_varlist = data_frame.iloc[1:].to_numpy()
            experimental_data.append(new_experiment_data_varlist)
            new_experiment_data_mask_varlist = (
                data_frame.iloc[1:].notna().to_numpy().astype(int)
            )
            experimental_data_mask.append(new_experiment_data_mask_varlist)

            # Generate inverted_variances
            variable_name_list = list(ordered_varlist_input.get_state().keys())
            if self._use_algebraic_variables:
                variable_name_list.extend(
                    list(ordered_varlist_input.get_algebraic().keys())
                )
            inverted_variances_varlist = []
            inverted_scaled_variances_varlist = []
            for var_name in variable_name_list:
                var = ordered_varlist_input[var_name]
                scaled_variance = var.variance / var._get_scaling_constants()[0] ** 2
                inverted_variances_varlist.append(
                    1.0 / (np.full(len(time_grid_unique) - 1, var.variance))
                )
                inverted_scaled_variances_varlist.append(
                    1.0 / (np.full(len(time_grid_unique) - 1, scaled_variance))
                )
            inverted_variances_array = np.column_stack(inverted_variances_varlist)
            inverted_scaled_variances_array = np.column_stack(
                inverted_scaled_variances_varlist
            )

            list_inverted_variances.append(inverted_variances_array)
            list_inverted_scaled_variances.append(inverted_scaled_variances_array)

            size_simulation_output.append(inverted_variances_array.shape)

        # Calculate experiments_weights
        self.list_simulators: Sequence[Simulator] = list_simulators

        array_data = np.concatenate(experimental_data)
        all_measurements_names_list = list(ordered_varlist_input.get_state().keys())
        if self._use_algebraic_variables:
            all_measurements_names_list.extend(
                list(ordered_varlist_input.get_algebraic().keys())
            )

        all_measurements_names = np.array(all_measurements_names_list)

        index_columns_with_all_nans = np.isnan(array_data).all(axis=0)

        max_time_grid = max(list_timegrid_length)
        experiments_weights = []
        for time_grid_length, size_simulation in zip(
            list_timegrid_length, size_simulation_output
        ):
            experiments_weights.append(
                np.full(size_simulation, max_time_grid / time_grid_length)
            )

        """
        This list holds nested arrays with all experimental data. If data is not available for the time_stamp,
        it's replaced with 0. It has follwing form:
        [exp1_var1_time1, exp1_var2_time1, exp1_varN_time1, exp1_var1_time2 ... , exp1_varN_timeN, exp2_var1_time1 ...]
        """
        self.array_data = np.nan_to_num(array_data[:, ~index_columns_with_all_nans])
        self.array_data_mask = np.concatenate(experimental_data_mask)[
            :, ~index_columns_with_all_nans
        ]
        self.array_data_mx = ca.MX.sym("array_data", self.array_data.shape)
        self.array_data_mask_mx = ca.MX.sym("array_data_mask", self.array_data.shape)

        self.names_of_measurements: list[str] = all_measurements_names[
            ~index_columns_with_all_nans
        ].tolist()

        self.index_measurements_in_sim = []
        for name in self.names_of_measurements:
            try:
                index = self.list_simulators[0].mapping_state_variables[name]
            except KeyError:
                index = (
                    len(self.list_simulators[0].mapping_state_variables)
                    + self.list_simulators[0].mapping_algebraic_variables[name]
                )
            self.index_measurements_in_sim.append(index)

        # Inverted variances provided weightning matrix for PE problem
        array_inverted_variance: np.ndarray = np.concatenate(list_inverted_variances)[
            :, ~index_columns_with_all_nans
        ]
        self.array_inverted_std = np.sqrt(array_inverted_variance)

        array_inverted_scaled_variance: np.ndarray = np.concatenate(
            list_inverted_scaled_variances
        )[:, ~index_columns_with_all_nans]
        self.array_inverted_scaled_std = np.sqrt(array_inverted_scaled_variance)

        self.experiments_weights: np.ndarray = np.concatenate(experiments_weights)

        self.generate_simulate_all_functions()

    def optimize(
        self,
        scale=None,
        objective_function="wls",
        *,
        scale_experiments=False,
        reuse_solver=False,
    ) -> dict[str, ca.DM]:
        """Solves optimization problem. Scaling decreases amount of iterations,
        and should always almost be used
        """
        self._setup_experiments_scale(scale_experiments)
        return PE_base.optimize(self, scale, objective_function)

    @_consistent_scaling_decorator
    def plot_simulation(
        self,
        supplied_parameters: dict[str, float] | None = None,
        experiment_names: list[str] | None = None,
        savefig: bool = False,
        algebraic: bool = True,
        plot: bool = True,
    ) -> list[VariableList]:  # noqa: E501
        """Plots experimental points against simulated trajectories, first line, initial guess, than supplied values"""

        if experiment_names is None:
            experiment_names = []
            for index, _ in enumerate(self.list_input_varlist):
                experiment_names.append(f"EXP_{index}")

        if not len(experiment_names) == len(self.list_input_varlist):
            raise ValueError(
                "Length of experiment names is not same as ammount of experiments"
            )

        for input_varlist, simulator, exp_name in zip(
            self.list_input_varlist,
            self.list_simulators,
            experiment_names,
        ):
            res_guess = simulator.simulate(
                algebraic=algebraic,
                recalculate_algebraic=True,
                unfixed_variables=dict(zip(self.varlist_decision.keys(), self.guess)),
            )[2]

            if supplied_parameters is not None:
                res_supplied = simulator.simulate(
                    algebraic=algebraic,
                    recalculate_algebraic=True,
                    unfixed_variables=supplied_parameters,
                )[2]

            if plot:
                if pd.get_option("plotting.backend") == "plotly":
                    if supplied_parameters is None:
                        fig = res_guess.dataframe.plot(markers=True)
                    else:
                        fig = res_supplied.dataframe.plot(markers=True)
                    fig.show()
                else:
                    axes = res_guess.plot(
                        prefix="GUESS ", color="blue", algebraic=algebraic, show=False
                    )
                    if supplied_parameters is not None:
                        axes = res_supplied.plot(
                            prefix="FINAL ", ax=axes, color="red", algebraic=algebraic
                        )
                    axes[0].set_title(exp_name)
                    input_varlist.plot(
                        ax=axes,
                        marker="x",
                        color="black",
                        prefix="EXP ",
                        linestyle="None",
                        algebraic=algebraic,
                        show=True,
                    )

        if supplied_parameters is None:
            return [res_guess]
        else:
            return [res_guess, res_supplied]


class ParameterEstimationNLE(PE_base):
    def __init__(
        self,
        model: Model,
        variable_lists: list[VariableList],
        simulator_settings=None,
        simulator_name="rootfinder",
        *,
        use_simulator_bounds=None,
        SimulatorClass=SimulatorNLE,
    ) -> None:
        if use_simulator_bounds is not None:
            warn(
                "use_simulator_bounds is not used anymore and will be ignored",
                FutureWarning,
                2,
            )
        super().__init__(model, variable_lists, simulator_name, simulator_settings)

        self._setup_simulator(SimulatorClass)
        self.logger.debug(
            "Created Optimizer object: \n Data Shape {} \n Desicion Variables {}".format(
                self.array_data.shape, self.varlist_decision.get_variable_name()
            )
        )
        self._setup_initialization()
        self._setup_direct_optimization("PE")

        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {"max_iter": 300},
        }
        # Set default objective
        self._objective: Callable[[], tuple[ca.MX | ca.DM, ca.MX | ca.DM]] = (
            self._objective_ols
        )

        self._setup_experiments_scale(False)
        self.setup_regularization(0, np.zeros((len(self.varlist_decision), 1)))

    @_consistent_scaling_decorator
    def _setup_simulator(self, SimulatorClass: SimulatorNLE) -> None:
        # It's not checked if all supplied varlist have same states etc.
        if not issubclass(SimulatorClass, SimulatorNLE):
            raise NotImplementedError("Provided simulator_class is not supported")

        self._setup_varlist_decision()

        list_data_mask = []
        list_simulators = []
        list_data = []
        list_inverted_variances = []
        list_inverted_scaled_variances = []

        for varlist_input in self.list_input_varlist:
            ordered_varlist_input = self.model.varlist(varlist_input)
            varlist_data = []
            varlist_data_mask = []
            varlist_variance = []
            varlist_scaled_variance = []
            for var in ordered_varlist_input.values():
                if isinstance(var, VariableControl):
                    if not var.fixed:
                        raise NotImplementedError
                if isinstance(var, VariableParameter):
                    # Avoid situations, where parameters unfixed differently
                    # in each self.list_input_varlist
                    fixed = var.name not in self.varlist_decision.keys()
                    var.fixed = fixed

            simulator = SimulatorClass(
                self.model,
                ordered_varlist_input,
                self.simulator_settings,
                self.simulator_name,
            )
            list_simulators.append(simulator)

            for var in ordered_varlist_input.get_algebraic().values():
                if var.value[0] is None or np.isnan(var.value[0]):
                    varlist_data.append(np.nan)
                    varlist_data_mask.append(0.0)
                else:
                    scaled_value = var.scale_from_original(var.value[0])
                    varlist_data.append(scaled_value)
                    varlist_data_mask.append(1.0)

                scaled_variance = var.variance / var._get_scaling_constants()[0] ** 2
                varlist_variance.append(1.0 / var.variance)
                varlist_scaled_variance.append(1.0 / scaled_variance)
            list_data.append(varlist_data)
            list_data_mask.append(varlist_data_mask)
            list_inverted_variances.append(varlist_variance)
            list_inverted_scaled_variances.append(varlist_scaled_variance)

        self.list_simulators: list[SimulatorNLE] = list_simulators

        array_data = np.array(list_data)
        all_measurements_names = np.array(
            list(ordered_varlist_input.get_algebraic().keys())
        )

        index_columns_with_all_nans = np.isnan(array_data).all(axis=0)

        self.array_data = np.nan_to_num(array_data[:, ~index_columns_with_all_nans])
        self.array_data_mask = np.array(list_data_mask)[:, ~index_columns_with_all_nans]

        self.array_data_mx = ca.MX.sym("array_data", self.array_data.shape)
        self.array_data_mask_mx = ca.MX.sym("array_data_mask", self.array_data.shape)

        self.names_of_measurements: list[str] = all_measurements_names[
            ~index_columns_with_all_nans
        ].tolist()
        array_inverted_variance = np.array(list_inverted_variances)[
            :, ~index_columns_with_all_nans
        ]
        self.array_inverted_std = np.sqrt(array_inverted_variance)
        array_inverted_scaled_variance = np.array(list_inverted_scaled_variances)[
            :, ~index_columns_with_all_nans
        ]
        self.array_inverted_scaled_std = np.sqrt(array_inverted_scaled_variance)

        self.index_measurements_in_sim = []
        for name in self.names_of_measurements:
            index = self.list_simulators[0].mapping_algebraic_variables[name]
            self.index_measurements_in_sim.append(index)

        self.generate_simulate_all_functions()

    @_consistent_scaling_decorator
    def calculate_inference_bounds(
        self,
        dict_of_params: dict,
        dict_of_responses: dict,
        dict_of_controls: dict,
        dict_of_artificial_controls: dict = None,
        rng: np.random.Generator = None,
    ):
        """Method to calculate one-dimensional inference bounds for a given dict of responses.
        If dict_of_artificial_controls is supplied, artificial data is generated and used for instead
        of experimental data. If rng is supplied, the given rng is used for articficial data generation.

        Parameters
        ----------
        dict_of_params : dict
            keys: parameter name
            values: corresponding parameter value
        dict_of_responses : dict
            keys: response names
            values: corresponding response variance
            Only used for artificial data generation.
            If values are None, the defualt response variance value is utilized.
        dict_of_controls : dict
            keys: control names
            values: corresponding list of necessary information about controls with
            list = list([lower bound: float, upper bound: float, number of points: int])
        dict_of_artificial_controls : dict, optional
            dict_of_controls used for artificial data generation, by default None
        rng : np.random.Generator, optional
            rng used for artifical data generation, by default None

        Output
        ------
        inference_results : dict
            Contains all inference bounds related results
        exp_data : dict
            Contains experimental data OR generated artifical data
        sim_results : dict
            Contains data of the simulation needed for inference computation
            Content of same dimension as content of inference_results
        """
        import scipy.stats

        def convert_varlist_to_data_dictionary(
            var_list_list: list[VariableList],
            dict_of_responses: dict,
            dict_of_controls: dict,
        ):
            sim_data = {}
            for var_name in [*dict_of_controls.keys(), *dict_of_responses.keys()]:
                var_values = []
                for simulation in var_list_list:
                    var = simulation[var_name]
                    var_values.append(var.scale_to_original(var.value[0]))
                sim_data[var_name] = np.array(var_values)

            return sim_data

        def generate_simulation_data(
            template_varlist: VariableList,
            dict_of_params: dict,
            dict_of_responses: dict,
            dict_of_controls: dict,
            perturbate: bool,
            rng: np.random.Generator = None,
        ):

            for param, param_value in dict_of_params.items():
                template_varlist[param].value = param_value
            for key, variance in dict_of_responses.items():
                if variance is not None:
                    template_varlist[key].variance = variance

            generated_var_lists, true_parameters, _ = (
                tools.generate_artificial_data_from_grid_nle(
                    self.model,
                    template_varlist,
                    dict_of_controls,
                    perturbate=perturbate,
                    rng=rng,
                    measurement_names=dict_of_responses.keys(),
                    keep_in_bounds=True,
                )
            )

            for parameter in dict_of_params:
                for simulation in generated_var_lists:
                    simulation[parameter].fixed = False

            sim_data = convert_varlist_to_data_dictionary(
                generated_var_lists, dict_of_controls, dict_of_responses
            )

            return generated_var_lists, sim_data

        if dict_of_artificial_controls is None:
            artificial_mode = False
            experimental_data = copy.deepcopy(self.list_input_varlist)
            exp_data = convert_varlist_to_data_dictionary(
                experimental_data,
                dict_of_responses,
                dict_of_controls,
            )
        else:
            artificial_mode = True
            experimental_data, exp_data = generate_simulation_data(
                copy.deepcopy(self.list_input_varlist[0]),
                dict_of_params,
                dict_of_responses,
                dict_of_artificial_controls,
                True,
                rng=rng,
            )

        variable_list_real, sim_data = generate_simulation_data(
            copy.deepcopy(self.list_input_varlist[0]),
            dict_of_params,
            dict_of_responses,
            dict_of_controls,
            False,
        )

        OLS = {}
        if artificial_mode:
            pe_artificial = ParameterEstimationNLE(self.model, experimental_data)
        else:
            pe_artificial = self

        residuals = pe_artificial.calculate_objective_and_residual(dict_of_params)[
            "residuals"
        ]
        OLS_values = np.diag(residuals.T @ residuals)
        OLS = dict(zip(self.names_of_measurements, OLS_values))

        jac = pe_artificial.calculate_sensitivity_and_fim(
            dict_of_params, list(dict_of_params.keys())
        )["jac_sorted"]

        pe_grid = ParameterEstimationNLE(self.model, variable_list_real)
        jac_grid = pe_grid.calculate_sensitivity_and_fim(
            dict_of_params, list(dict_of_params.keys())
        )["jac_sorted"]

        len_exp = len(experimental_data)
        len_param = len(dict_of_params)
        fisher95 = scipy.stats.f(len_param, self.dof).ppf(0.95)

        inference_results = {}
        for control in dict_of_controls:
            inference_results[control] = np.array(sim_data[control])

        for response in dict_of_responses:
            inference_results[response] = {}
            s = np.sqrt(OLS[response] / self.dof)
            R = np.linalg.qr(jac[response], mode="reduced")[1]
            bound = (
                s
                * np.linalg.norm(jac_grid[response] @ np.linalg.inv(R), axis=1)
                * np.sqrt(len_param * fisher95)
            )

            inference_results[response]["s"] = s
            inference_results[response]["R"] = R
            inference_results[response]["bound"] = bound
            inference_results[response]["lower bound"] = sim_data[response] - bound
            inference_results[response]["simulation"] = sim_data[response]
            inference_results[response]["upper bound"] = sim_data[response] + bound

        return inference_results, exp_data, sim_data


class ParameterEstimationNLE_control(ParameterEstimationNLE):
    @_consistent_scaling_decorator
    def _setup_simulator(self):
        # It's not checked if all supplied varlist have same states etc.
        for var in self.list_input_varlist[0].values():
            if isinstance(var, VariableAlgebraic):
                self.varlist_algebraic.add_variable(var)
            elif isinstance(var, VariableParameter):
                self.varlist_parameter.add_variable(var)
                if var.fixed is False:
                    self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableControl):
                self.varlist_control.add_variable(var)

        self.num_parameters = len(self.varlist_decision)

        list_simulators = []
        self.array_data = []
        self.array_data_mask = []

        self.array_controls = []
        self.array_controls_casadi = []

        for index, varlist_input in enumerate(self.list_input_varlist):
            new_varlist = copy.deepcopy(varlist_input)
            for var in varlist_input.values():
                if isinstance(var, VariableControl):
                    if var.fixed is False:
                        new_varlist.pop(var.name)
                        new_var = VariableControl(
                            f"{var.name}_exp{index}",
                            var.value[0],
                            var.lower_bound,
                            var.upper_bound,
                            var.opc_ua_id,
                        )
                        new_var.fixed = False
                        self.varlist_decision.add_variable(new_var)
                        self.array_controls_casadi.append(new_var.casadi_var)
                        new_varlist[var.name] = new_var
                        self.array_controls.append(var.value[0])

            # for var in varlist_input.values():
            #     if isinstance(var, VariableControl):
            #         var.fixed = True

            simulator = SimulatorNLE(
                self.model,
                new_varlist,
                self.simulator_settings,
                self.simulator_name,
            )
            list_simulators.append(simulator)

            for var in varlist_input.values():
                if isinstance(var, VariableAlgebraic):
                    if var.value[0] is None or np.isnan(var.value[0]):
                        self.array_data.append(0)
                        self.array_data_mask.append(0)
                    else:
                        self.array_data.append(var.value[0])
                        self.array_data_mask.append(1)

        self.list_simulators: list[SimulatorNLE] = list_simulators
        self.array_data = ca.DM(self.array_data)
        self.array_controls = ca.DM(self.array_controls)
        self.array_data_mask = np.array(self.array_data_mask)
        self.array_controls_casadi = ca.vcat(self.array_controls_casadi)

    def _objective__(self):
        array_simulation = None

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate_fast()

            if array_simulation is None:
                array_simulation = res_simulation["x"]
            else:
                array_simulation = ca.vertcat(array_simulation, res_simulation["x"])

        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (array_simulation - self.array_data) * self.array_data_mask
        error_controls = self.array_controls_casadi - self.array_controls
        objective = ca.sum1(error**2) + ca.sum1(error_controls**2)

        return objective

    @_consistent_scaling_decorator
    def optimize(self, scale=None, objective_function="ols", direct_optimization=False):
        if objective_function == "wls":
            self._objective = partial(
                self._objective_wls, direct_optimization=direct_optimization
            )
        elif objective_function == "ols":
            self._objective = partial(
                self._objective_ols, direct_optimization=direct_optimization
            )

        res = self._optimize(scale, direct_optimization=direct_optimization)
        res["all"] = res["x"]
        res["x"] = res["all"][0 : self.num_parameters]
        res["p"] = res["all"][self.num_parameters :]
        return res
