from __future__ import annotations

import copy
import logging
from abc import abstractmethod
from collections.abc import Callable
from typing import Sequence
from itertools import combinations

import casadi as ca
import numpy as np
import pandas as pd
from scipy import linalg
from tqdm import tqdm

from par_est import (
    Model,
    Simulator,
    SimulatorNLE,
    VariableAlgebraic,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableList,
    VariableParameter,
    VariableState,
    utilities,
)


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

        if not isinstance(variable_lists, list):
            raise (Exception("Variable list should be nested of type list"))
        self.logger: logging.Logger = logging.getLogger(__name__)
        self.model: Model = model

        # Deepcopy is used to avoid manipulating input variable list
        self.list_input_varlist: list[VariableList] = copy.deepcopy(variable_lists)

        # Each varlist holds respective variables
        self.varlist_decision: VariableList = VariableList()
        self.varlist_parameter: VariableList = VariableList()
        self.varlist_control: VariableList = VariableList()
        self.varlist_state: VariableList = VariableList()
        self.varlist_algebraic: VariableList = VariableList()

        self.simulator_name: str = simulator_name
        self.simulator_settings: dict | None = simulator_settings

        self.list_simulators: Sequence[Simulator | SimulatorNLE]
        self.mapping_simulator_decisions: list = []

    @abstractmethod
    def optimize(self, scale):
        # Runs optimization once
        raise (NotImplementedError)

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
                guess.append(var.guess)
            lower_bound.append(var.lower_bound)
            upper_bound.append(var.upper_bound)

        self.guess: np.ndarray = np.array(guess)
        self.lower_bound: np.ndarray = np.array(lower_bound)
        self.upper_bound: np.ndarray = np.array(upper_bound)
        self.logger.debug(
            f"Initialized:\nguess {self.guess}\nlower_bound {self.lower_bound}\nupper_bound {self.upper_bound}"
        )

    def _setup_scaling(self, scale: bool = False) -> None:
        """Scaling of decision variables should be done before setting a solver and solver settings.
        Sets scaling variables in optimizer and simulator.
        Scaling is used to correctly spread the weight of decision variables in optimizaiton function.
        By default, variable guess is used as a scaling value.
        TODO: Whole loop can be replaced by simple np.where, isn't it?
        """
        if scale:
            self.scaling: np.ndarray | int = np.where(self.guess == 0, 1, self.guess)

            for simulator, mapping in zip(
                self.list_simulators, self.mapping_simulator_decisions
            ):
                simulator._reset_scaling()

                if isinstance(self, ParameterEstimationNLE):
                    for index_simulator, index_decision in mapping.items():
                        current_guess = self.guess[index_decision]
                        if current_guess == 0:
                            simulator.scaling[index_simulator] = 1
                        else:
                            simulator.scaling[index_simulator] = current_guess

                else:
                    # this loop is used because simple "for loop" fails for ca.MX vector
                    if isinstance(simulator, SimulatorNLE):
                        independent_variables = simulator._independent_variables
                    elif isinstance(simulator, Simulator):
                        independent_variables = simulator._independent_variables[0]
                    else:
                        raise NotImplementedError

                    for count in range(independent_variables.size()[0]):
                        var = independent_variables[count]
                        if var.is_symbolic():
                            if var.name() in self.varlist_decision:
                                current_guess = self.varlist_decision[var.name()].guess
                                if current_guess == 0:
                                    simulator.scaling[count] = 1
                                else:
                                    simulator.scaling[count] = current_guess
                            else:
                                simulator.scaling[count] = 1
        else:
            self.scaling = 1
            for simulator in self.list_simulators:
                simulator._reset_scaling()

    @abstractmethod
    def _objective(self) -> tuple[ca.MX | ca.DM, ca.MX | ca.DM]:
        """Returns a way to calculate and objective. Dependent on optimization type."""
        raise (NotImplementedError)

    def _optimize(self, scale: bool) -> dict[str, ca.DM | ca.MX]:
        """Runs optimizer, uses scaling if needed. Returned values is scaled back.
        Scaling should be done before setting a solver and solver settings."""
        self._setup_scaling(scale)

        self.solver: ca.Function = ca.nlpsol(
            "solver",
            self.solver_name,
            {
                "x": self.varlist_decision.get_casadi_variables(),
                "f": self._objective()[0],
            },
            self.solver_settings,
        )

        lb_scaled = self.lower_bound / self.scaling
        ub_scaled = self.upper_bound / self.scaling

        # Scaling of negative numbers requires a switch bounds
        for index, (lb, ub) in enumerate(zip(lb_scaled, ub_scaled)):
            if lb > ub:
                lb_scaled[index] = ub
                ub_scaled[index] = lb

        res_solver = self.solver(
            x0=self.guess / self.scaling,
            lbx=lb_scaled,
            ubx=ub_scaled,
        )

        res_solver["x"] = res_solver["x"] * self.scaling

        res_dict = {}
        for solution, var_name in zip(
            res_solver["x"].toarray(), list(self.varlist_decision.keys())
        ):
            res_dict[var_name] = float(solution[0])

        res_solver["x_dict"] = res_dict

        return res_solver

    def map_objective(self, plot: bool = True) -> None:
        """Calculate objective function for different values of parameters and plot, if needed.
        Currently support only 3 unfixed decision variables."""
        import matplotlib.pyplot as plt

        decision_variables = self.varlist_decision.get_casadi_variables()
        if decision_variables.shape == (3, 1):
            self._setup_scaling(False)
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

        list_startpoint = utilities.make_startpoints(hammersley_seeds, num_initials)
        results_with_f = []

        # Optimizer settings and guess are overwritten for multistart, and then returned back
        initial_guess = copy.deepcopy(self.guess)
        initial_settings = copy.deepcopy(self.solver_settings)

        self.solver_settings = {
            "verbose": False,
            "print_time": False,
            "ipopt": {
                "hessian_approximation": "limited-memory",
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
        self.scaling = initial_guess
        return result_list_sorted

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
                    result = simulation.generate_exp_data(True, True, set_of_bounds)
                    if plot:
                        result._get_varlist_to_plot(True).dataframe.plot(
                            subplots=True, title=str(bound_dictionary)
                        )
                    results.append(result)
                except Exception:
                    print(f"Failed for these bounds: {bound_dictionary}")
        else:
            raise NotImplementedError


class ParameterEstimation(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_list: list[VariableList],
        simulator_name: str = "idas",
        simulator_settings: dict | None = None,
        *,
        use_idas_constraints: bool = False,
        use_algebraic_vars: bool = False,
        recalculate_algebraic: bool = False,
    ):
        super().__init__(
            model,
            variable_list,
            simulator_name,
            simulator_settings,
        )

        if use_algebraic_vars:
            objective = self._objective_alg
            raise NotImplementedError()
        else:
            objective = self._objective_state
        self._objective: Callable[[], tuple[ca.MX | ca.DM, ca.MX | ca.DM]] = objective

        self._setup_simulator(
            use_idas_constraints=use_idas_constraints,
            use_algebraic_vars=use_algebraic_vars,
            recalculate_algebraic=recalculate_algebraic,
        )

        self.logger.debug(
            "Created Optimizer object: \n Data Shape {} \n Desicion Variables {}".format(
                self.experimental_data.shape, self.varlist_decision.get_variable_name()  # type: ignore
            )
        )
        self._setup_initialization()

        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {"hessian_approximation": "limited-memory", "max_iter": 300},
        }

        for sim in self.list_simulators:
            sim.calculate_algebraic_initials(apply_intials=True)

    def _setup_simulator(
        self,
        *,
        use_idas_constraints: bool,
        use_algebraic_vars: bool,
        recalculate_algebraic: bool,
    ) -> None:
        # It's not checked if all supplied varlist have same states etc.
        for var in self.list_input_varlist[0].values():
            if isinstance(var, VariableState):
                self.varlist_state.add_variable(var)
            elif isinstance(var, VariableParameter):
                self.varlist_parameter.add_variable(var)
                if var.fixed is False:
                    self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableControl):
                self.varlist_control.add_variable(var)

        # Lists used to calculate experiments_weights
        list_timegrid_length = []
        size_simulation_output = []
        list_simulators = []
        experimental_data = []
        experimental_data_mask = []
        inverted_variances: list[np.ndarray] = []

        for varlist_input in self.list_input_varlist:
            # Create a time_grid, that "stops" at every experimental data, for every state variable
            if not varlist_input.get_common_origin(
                strict=True, variable_type=VariableState
            ):
                raise (
                    ValueError(
                        f"Not all State Variables in one experiment have same time0, so simulations cannot be initialized:\n{varlist_input}"
                    )
                )
            data_frame = pd.DataFrame()

            for variable_name in self.model.varlist_all.keys():
                try:
                    var = varlist_input[variable_name]
                except KeyError:
                    continue

                if isinstance(var, VariableState) or (
                    isinstance(var, VariableAlgebraic) and use_algebraic_vars
                ):
                    data_frame = data_frame.join(var.dataframe, how="outer")

                elif isinstance(var, VariableControl):
                    var.fixed = True
                    if isinstance(var, VariableControlPiecewiseConstant):
                        var.fixed = True
                        data_frame = data_frame.join(var.dataframe, how="outer")
                        # Column should be dropped, because it's needed only for unique timestamp
                        data_frame.drop(columns=var.name, inplace=True)

            time_grid_unique = (
                (data_frame.index - data_frame.index[0]).total_seconds().tolist()
            )

            list_timegrid_length.append(float(len(time_grid_unique)))

            list_simulators.append(
                Simulator(
                    self.model,
                    np.array(time_grid_unique),
                    varlist_input,
                    self.simulator_name,
                    self.simulator_settings,
                    use_idas_constraints=use_idas_constraints,
                    recalculate_algebraic=recalculate_algebraic,
                )
            )

            # Generate an array (experiment_data_varlist) with Experimental data with the same dimensions as simulation results.
            new_experiment_data_varlist = (
                data_frame.iloc[1:].fillna(0).to_numpy().flatten()
            )
            experimental_data.append(new_experiment_data_varlist)
            new_experiment_data_mask_varlist = (
                data_frame.iloc[1:].notna().to_numpy().flatten().astype(int)
            )
            experimental_data_mask.append(new_experiment_data_mask_varlist)

            # Generate inverted_variances
            if use_algebraic_vars:
                variable_name_list = list(
                    [
                        *self.model.varlist_state.keys(),
                        *self.model.varlist_algebraic.keys(),
                    ]
                )
            else:
                variable_name_list = list(self.model.varlist_state.keys())
            inverted_variances_varlist = []
            for var_name in variable_name_list:
                var = varlist_input[var_name]
                inverted_variances_varlist.append(
                    1.0 / (np.full(len(time_grid_unique) - 1, var.variance))
                )
            inverted_variances_array = np.column_stack(
                inverted_variances_varlist
            ).flatten()
            inverted_variances.append(inverted_variances_array)

            size_simulation_output.append(len(inverted_variances_array))

        # Calculate experiments_weights
        self.list_simulators: Sequence[Simulator] = list_simulators

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
        self.experimental_data: np.ndarray = np.concatenate(experimental_data)
        # This list holds nested arrays with True or False values, specifying whether experimental data should be
        # used in an optimization function or not. Same form as self.experimental_data

        self.experimental_data_mask: np.ndarray = np.concatenate(experimental_data_mask)
        # Inverted variances provided weightning matrix for PE problem
        self.inverted_variances: np.ndarray = np.concatenate(inverted_variances)

        self.experiments_weights: np.ndarray = np.concatenate(experiments_weights)

    def _objective_state(self) -> tuple[ca.MX | ca.DM, ca.MX | ca.DM]:
        array_simulation = None

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate()
            if array_simulation is None:
                array_simulation = res_simulation["xf"][:]
            else:
                array_simulation = ca.vertcat(array_simulation, res_simulation["xf"][:])
        # res_simulation = self.list_simulators[0].simulate()
        # array_simulation = res_simulation["xf"][:]

        # for simulator in self.list_simulators[1:]:
        #     res_simulation = simulator.simulate()
        #     array_simulation = ca.vertcat(array_simulation, res_simulation["xf"][:])

        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (
            array_simulation - self.experimental_data
        ) * self.experimental_data_mask
        objective = ca.sum1(
            self.experiments_scale * self.inverted_variances * (error**2)
        )

        return objective, error

    def _objective_alg(self) -> tuple[ca.MX | ca.DM, ca.MX | ca.DM]:
        array_simulation = None

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate()

            if array_simulation is None:
                res_all = ca.vertcat(res_simulation["xf"], res_simulation["zf"])
                array_simulation = res_all[:]
            else:
                res_all = ca.vertcat(res_simulation["xf"], res_simulation["zf"])
                array_simulation = ca.vertcat(array_simulation, res_all[:])

        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (
            array_simulation - self.experimental_data
        ) * self.experimental_data_mask
        objective = ca.sum1(
            self.experiments_scale * self.inverted_variances * (error**2)
        )

        return objective, error

    def optimize(self, scale=True, *, scale_experiments=False) -> dict[str, ca.DM]:
        """Solves optimization problem. Scaling decreases amount of iterations,
        and should always almost be used
        """
        if scale_experiments:
            experiments_scale: int | np.ndarray = self.experiments_weights
        else:
            experiments_scale = 1

        # This attribute is used while calculating Objective, and is either 1 or self.experiments_weights
        # It's used to make some experiments as valuable as others, even if they have less experimental points
        # So if you supply 2 experiments one with 10 and another with 20 time_stamps, effect of each experimental
        # point of second experiments on objective function is decreased by 2
        self.experiments_scale: int | np.ndarray = experiments_scale
        res = self._optimize(scale)
        return res

    def plot_simulation(
        self,
        supplied_parameters: list[float] | None = None,
        experiment_names: list[str] | None = None,
        savefig: bool = False,
        algebraic: bool = True,
        plot: bool = True,
    ) -> list[VariableList]:  # noqa: E501
        """Plots experimental points against simulated trajectories, first line, initial guess, than supplied values"""
        self._setup_scaling(False)

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
            res_guess = simulator.generate_exp_data(
                algebraic=algebraic,
                recalculate_algebraic=True,
                unfixed_variables=self.guess,
            )

            if supplied_parameters is not None:
                res_supplied = simulator.generate_exp_data(
                    algebraic=algebraic,
                    recalculate_algebraic=True,
                    unfixed_variables=supplied_parameters,
                )

            if plot:
                if pd.get_option("plotting.backend") == "plotly":
                    if supplied_parameters is None:
                        fig = res_guess.dataframe.plot(markers=True)
                    else:
                        fig = res_supplied.dataframe.plot(markers=True)
                    fig.show()
                else:
                    axes = res_guess.plot(
                        prefix="GUESS ", color="blue", algebraic=algebraic
                    )
                    if supplied_parameters is not None:
                        axes = res_supplied.plot(
                            prefix="FINAL ", ax=axes, color="red", algebraic=algebraic
                        )
                    input_varlist.plot(
                        ax=axes,
                        marker="x",
                        color="black",
                        prefix="EXP ",
                        linestyle="None",
                        algebraic=algebraic,
                    )
                    axes[0].set_title(exp_name)

        if supplied_parameters is None:
            return [res_guess]
        else:
            return [res_guess, res_supplied]


class OptimalExperimentalDesign(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_list: list[VariableList],
        time_grid_relative: np.ndarray,
        simulator_name: str = "idas",
        simulator_settings: dict = None,
        *,
        reinitialize_algebraic: bool = False,
    ) -> None:
        super().__init__(model, variable_list, simulator_name, simulator_settings)

        # User specified time_grid is used for initilizaiton of Simulators
        self.time_grid_original: np.ndarray = time_grid_relative

        self._setup_simulator()
        self._setup_initialization()

        # User specified time grid might not include every time_stamp of VariableControlPiecewiseConstant
        # So the corrected time_grid of Simulator is used in optimization
        self.time_grid_modified = self.list_simulators[0].time_grid_relative

        self.solver_name: str = "ipopt"
        self.solver_settings: dict = {
            "verbose": False,
            # "monitor": ["nlp_grad_f", "nlp_f"],
            "ipopt": {
                "hessian_approximation": "limited-memory",
                "max_iter": 100,
                # "print_level": 6,
            },
        }
        if reinitialize_algebraic:
            for sim in self.list_simulators:
                sim.calculate_algebraic_initials(apply_intials=True)

    def _setup_simulator(self, *, use_idas_constraints: bool = False) -> None:
        """Initializes simulator class. Parameter variables are fixed, and an index of an unfixed
        parameter is saved in self.select_independent list.
        This list is used during the calculation of the objective, to ignore jacobian of fixed parameters.
        self.index_all_states is used additionaly to self.select_independent list to get required jacobian.
        """
        parameter_values = []
        select_independent = []
        inverted_variances = []

        for variable_name in self.model.varlist_all.keys():
            try:
                var = self.list_input_varlist[0][variable_name]
            except KeyError:
                continue

            if isinstance(var, VariableControl):
                if not var.fixed:
                    if isinstance(var, VariableControlPiecewiseConstant):
                        for var_control in var.variable_list.values():
                            if not var_control.fixed:
                                self.varlist_decision.add_variable(var_control)
                    else:
                        self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableParameter):
                if var.fixed is False:
                    self.varlist_parameter.add_variable(var)
                    parameter_values.append(var.value[0])
                    var.fixed = True
                    # index has + 1 to account for tau variable, that is a parameter of a simulation.
                    index = list(self.model.varlist_independent).index(var.name) + 1
                    select_independent.append(index)

            elif isinstance(var, VariableState):
                inverted_variances.append(1 / var.variance)

        # [par_1, par_2, ... par_N]
        self.parameter_values: list[float] = parameter_values
        # Index of parameters, that are unfixed. It's used to select Jacobian for only this variables
        self.select_independent: list[int] = select_independent

        self.inverted_variances: np.ndarray = np.array(inverted_variances)
        self.index_all_states = list(range(len(self.model.varlist_state)))

        self.list_simulators: list[Simulator] = [
            Simulator(
                self.model,
                self.time_grid_original,
                self.list_input_varlist[0],
                self.simulator_name,
                self.simulator_settings,
                simulate_jac=True,
            )
        ]

    def _objective(
        self, analyze: bool = False, values: list[float] | None = None
    ) -> tuple[ca.MX | ca.DM, ca.MX | ca.DM]:
        """Calculates an A OED criteria, beacuse casadi cannot do other.
        "analyze" Flag and values are used for debugging.

        Args:
            analyze: used for debugging, calculates objective and covariance.
            values: this values are used as desicionb variables for calculating of the objective.
        """
        covariance_full = None
        # -1 ignores time point zero in self.time_grid
        num_time = len(self.time_grid_modified) - 1
        # +1 account for tau variable in Simulator class
        num_param = len(self.model.varlist_independent) + 1
        num_state = len(self.model.varlist_state)

        result_simulation = self.list_simulators[0].simulate_jac()
        result_jacobian = result_simulation["jac_xf_p"]

        # Used only for debugging
        if analyze is True:
            jacobian_full = None
            covariance_all = []
            objective_all = []

            self._setup_scaling(False)
            evaluate = ca.Function(
                "eval_fim",
                [self.varlist_decision.get_casadi_variables()],
                [result_jacobian],
            )
            evaluate_sim = ca.Function(
                "eval_fim",
                [self.varlist_decision.get_casadi_variables()],
                [result_simulation["xf"]],
            )
            if values is None:
                result_jacobian = evaluate(self.guess)
                result_sim = evaluate_sim(self.guess)
            else:
                result_jacobian = evaluate(values)
                result_sim = evaluate_sim(values)  # noqa: F841

        # Simulation returns jacobian that has to be split, to get jac at each time point.
        # list_jacobian_at_timepoint contains a list of that jacobians.
        split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param
        list_jacobian_at_timepoint = ca.horzsplit(result_jacobian, split_vector)

        # Jacobian is also scaled based on parameter values.
        parameter_scaling = ca.repmat(ca.DM(self.parameter_values).T, num_state, 1)

        # For jacobian at each timepoint, select only unfixed parameters.
        # Afterwards scale the jacobiand and calculate a parameter covariabce matrix.
        # Finally, sum covariances at every time point.
        for jacobian in list_jacobian_at_timepoint:
            jacobian_selected = jacobian.get(
                False, self.index_all_states, self.select_independent
            )
            jacobian_selected = jacobian_selected * parameter_scaling

            cov_at_timepoint = (
                jacobian_selected.T
                @ np.diag(self.inverted_variances)
                @ jacobian_selected
            )

            if analyze:
                covariance_all.append(cov_at_timepoint)
                objective_all.append(ca.trace(ca.inv(cov_at_timepoint)))

            if covariance_full is None:
                covariance_full = cov_at_timepoint
            else:
                covariance_full = covariance_full + cov_at_timepoint

        if analyze:
            for jacobian in list_jacobian_at_timepoint:
                jacobian_selected = jacobian.get(
                    False, self.index_all_states, self.select_independent
                )
                jacobian_selected = jacobian_selected * parameter_scaling

                if jacobian_full is None:
                    jacobian_full = jacobian_selected
                else:
                    jacobian_full = ca.vertcat(jacobian_full, jacobian_selected)

        error = ca.trace(ca.inv(covariance_full))

        if analyze:
            return error, covariance_full, jacobian_full, covariance_all, objective_all  # type: ignore

        return error, covariance_full

    def optimize(self, scale: bool = False) -> dict[str, ca.DM | ca.MX]:
        """Run optimization.

        Args:
            scale: scaling should be used as default, allows for faster convergence
        """

        # Scaling works unpredictably. It was shown during creation of VariableControlPiecewiseConstant
        if scale:
            raise NotImplementedError
        # Scaling decreases amount of iterations, but ipopt fails gradient check at big amount of timestamps

        return self._optimize(scale)

    def identifiability_analysis(self):
        """Taken from Erik/Diana Subset0. Many questions arrise about how it works."""
        (
            error,
            covariance_full,
            jacobian,
            covariance_all,
            objective_all,
        ) = self._objective(True)
        cond_threshold = 1000
        colin_threshold = 15

        states, parameters = jacobian.shape
        jacobian_original = jacobian
        if states < parameters:
            jacobian = np.pad(
                jacobian,
                ((0, parameters - states), (0, 0)),
                mode="constant",
                constant_values=0,
            )

        u, s, vh = np.linalg.svd(jacobian, False)

        # 2. rank determination of J
        values = abs(s)
        dimSVal = len(values)
        maxVal = np.max(values)
        # minVal = np.min(values)

        CondN_Sub = maxVal / values
        ColIdx_Sub = 1 / values

        smallval = []
        for i in range(0, dimSVal):
            if (
                np.abs(CondN_Sub[i]) <= cond_threshold
                and np.abs(ColIdx_Sub[i]) <= colin_threshold
            ):
                smallval.append(CondN_Sub[i])

        rank = len(smallval)

        # SummationSv = np.sum(values)
        # NeglectSv = np.sum(values[rank:] / SummationSv)

        # Determination of permutation matrix P by construction a RRQR of S
        Q, R, P = linalg.qr(jacobian, pivoting=True)

        # Use this for ranking
        IdentifOrd = P
        # Why
        # IdentifOrd.append(P)

        # Condition Number of J (Golub, 1996 & Hansen, 1998)
        # CondN = maxVal / minVal

        # Still to understand!!!
        # collinIndex.Colind
        # collinIndex.Sub

        Jsqr = jacobian**2

        ParNorm = np.empty(jacobian.shape[1])

        for i in range(0, jacobian.shape[1]):
            ParNorm[i] = np.sqrt(np.sum(Jsqr[:, i] / jacobian.shape[0]))

        # SensitivityOrder = np.argsort(ParNorm)[::-1]
        # B = np.sort(ParNorm)[::-1]

        # SensityOrd = B, SensitivityOrder

        # TotalVariance = np.sum(values ** 2)
        # NeglectVariance = np.sum(values[rank:] ** 2 / TotalVariance)

        unfix_parameters = []
        for index in range(rank):
            parameter_name = list(self.varlist_parameter.values())[
                IdentifOrd[index]
            ].name
            unfix_parameters.append(parameter_name)

        return unfix_parameters, error, covariance_full, jacobian_original


class ParameterEstimationNLE(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_lists: list[VariableList],
        simulator_settings=None,
        simulator_name="rootfinder",
        *,
        use_simulator_bounds=True,
        SimulatorClass=SimulatorNLE,
    ) -> None:
        super().__init__(model, variable_lists, simulator_name, simulator_settings)

        self._setup_simulator(use_simulator_bounds, SimulatorClass)
        self.logger.debug(
            "Created Optimizer object: \n Data Shape {} \n Desicion Variables {}".format(
                self.array_data.shape, self.varlist_decision.get_variable_name()
            )
        )
        self._setup_initialization()

        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {"max_iter": 300},
        }
        # Set default objective
        self._objective: Callable[
            [], tuple[ca.MX | ca.DM, ca.MX | ca.DM]
        ] = self._objective_ols

    def _setup_simulator(
        self, use_simulator_bounds: bool, SimulatorClass: SimulatorNLE
    ) -> None:
        # It's not checked if all supplied varlist have same states etc.
        if not issubclass(SimulatorClass, SimulatorNLE):
            raise NotImplementedError("Provided simulator_class is not supported")

        for var in self.list_input_varlist[0].values():
            if isinstance(var, VariableAlgebraic):
                self.varlist_algebraic.add_variable(var)
            elif isinstance(var, VariableParameter):
                self.varlist_parameter.add_variable(var)
                if var.fixed is False:
                    self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableControl):
                self.varlist_control.add_variable(var)

        array_data = []
        array_data_mask = []
        array_data_mask_new = []
        list_simulators = []
        list_simulator_mappings = []
        inverted_variances = []
        array_data_new = []

        for varlist_input in self.list_input_varlist:
            varlist_data = []
            varlist_data_mask = []
            for var in varlist_input.values():
                if isinstance(var, VariableControl):
                    if not var.fixed:
                        raise NotImplementedError

            simulator = SimulatorClass(
                self.model,
                varlist_input,
                self.simulator_settings,
                self.simulator_name,
                use_bounds=use_simulator_bounds,
            )
            list_simulators.append(simulator)

            list_simulator_mappings.append(self._setup_simulator_mapping(simulator))

            for var in varlist_input.values():
                if isinstance(var, VariableAlgebraic):
                    if var.value[0] is None or np.isnan(var.value[0]):
                        varlist_data.append(np.nan)
                        array_data.append(0.0)
                        array_data_mask.append(0.0)
                        varlist_data_mask.append(0.0)
                    else:
                        varlist_data.append(var.value[0])
                        array_data.append(var.value[0])
                        array_data_mask.append(1.0)
                        varlist_data_mask.append(1.0)

                    inverted_variances.append(1.0 / var.variance)
            array_data_new.append(varlist_data)
            array_data_mask_new.append(varlist_data_mask)

        self.list_simulators: list[SimulatorNLE] = list_simulators

        array_data_new = np.array(array_data_new)
        array_data_mask_new = np.array(array_data_mask_new)

        all_measurements_names = np.array(list(self.model.varlist_algebraic.keys()))
        data_columns_with_all_nans = np.isnan(array_data_new).all(axis=0)
        self.array_data_new = np.nan_to_num(array_data_new[:,~data_columns_with_all_nans])
        self.array_data_mask_new = array_data_mask_new[:,~data_columns_with_all_nans]
        self.names_of_measurements = all_measurements_names[~data_columns_with_all_nans].tolist()

        # List of dicts for each Simulation that shows, which index coresponds to each
        # simulator._independent_variables in self.varlist_ variable
        # For example [{1: 2}], {1: 3}]: in first simulator._independent_variables[1]
        # is the same variable as self.varlist_decision[2]
        self.mapping_simulator_decisions: list[dict[int, int]] = list_simulator_mappings

        self.inverted_variances: np.ndarray = np.array(inverted_variances)

        self.array_data: ca.DM = ca.DM(array_data)
        self.array_data_mask: np.ndarray = np.array(array_data_mask)

    def _setup_simulator_mapping(self, simulator: SimulatorNLE) -> dict[int, int]:
        names_variables_decision = list(self.varlist_decision.keys())

        independent_variables = simulator._independent_variables
        mapping_simulator_decisions = {}

        for count in range(independent_variables.size()[0]):
            var = independent_variables[count]
            if var.is_symbolic():
                if var.name() in names_variables_decision:
                    index = list(self.varlist_decision.keys()).index(var.name())
                    mapping_simulator_decisions[count] = index
                else:
                    raise NotImplementedError

        return mapping_simulator_decisions

    def _simulate_all(self):
        array_simulation = None

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate_sym()

            if array_simulation is None:
                array_simulation = res_simulation["x"]
            else:
                array_simulation = ca.vertcat(array_simulation, res_simulation["x"])

        return array_simulation

    def _objective_ols(self):
        # Ordinary Least Squares
        array_simulation = self._simulate_all()
        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (array_simulation - self.array_data) * self.array_data_mask
        objective = ca.sum1(error**2)

        return objective, error

    def _objective_wls(self):
        # Weighted Least Squares
        array_simulation = self._simulate_all()
        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (array_simulation - self.array_data) * self.array_data_mask
        objective = ca.sum1(self.inverted_variances * error**2)

        return objective, error

    def optimize(self, scale=False, objective_function="wls"):
        if objective_function == "wls":
            self._objective = self._objective_wls
        elif objective_function == "ols":
            self._objective = self._objective_ols

        variance = [[] for x in range(len(self.varlist_algebraic))]

        for sim in self.list_simulators:
            variance_i = sim.get_variance_array()
            variance.append(variance_i)

        return self._optimize(scale)

    def calculate_ols_value(self, parameters: dict[str, float]):
        decision_variables = self.varlist_decision.get_casadi_variables()

        self._setup_scaling(False)
        residuals_function = ca.Function(
            "objective",
            [decision_variables],
            [self._objective_ols()[0], self._objective_ols()[1]],
            ["x"],
            ["f", "residuals"],
        )

        selected_parameters = self.parameters_dict_to_list(parameters)
        result = residuals_function(x=selected_parameters)
        return result

    def parameters_dict_to_list(self, parameters: dict[str, float]) -> list[float]:
        """Takes a dictionary with {"par_name": par_value} and transforms to list
        corresponding to the order of self.varlist_decision variables"""
        selected_parameters: list[float] = []
        for var_name in parameters.keys():
            if var_name in self.varlist_decision.keys():
                selected_parameters.append(parameters[var_name])

        return selected_parameters

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
        self._setup_scaling()
        decision_variables = self.varlist_decision.get_casadi_variables()
        len_alg = len(self.model.varlist_algebraic)

        full_jacobian = []
        scaled_yao_jacobian = []
        variance = []
        select_independent: list[int] = []

        all_parameter_values = self.parameters_dict_to_list(parameters)
        selected_parameters_values = []

        for parameter_index, var_name in enumerate(self.varlist_decision.keys()):
            if parameter_names is not None:
                if var_name not in parameter_names:
                    continue
            index = list(self.model.varlist_independent).index(var_name) + len(
                self.model.varlist_algebraic
            )
            select_independent.append(index)
            selected_parameters_values.append(all_parameter_values[parameter_index])

        for sim_index, sim in enumerate(self.list_simulators):
            select_algebraic = np.reshape(
                self.array_data_mask, (len(self.list_simulators), len_alg)
            )[sim_index]
            select_algebraic = np.asarray(select_algebraic == 1).nonzero()[0]

            res_simulation = sim.simulate_sym_unfixed(all_parameter_values).toarray()

            variance_i = sim.get_variance_array()

            jacobian = ca.Function(
                "jacobian",
                [decision_variables],
                [sim.calculate_jac()["jac"]],
                ["x"],
                ["jac"],
            )

            jac = jacobian(all_parameter_values)

            jacobian_selected = jac.get(
                False, select_algebraic, select_independent
            ).toarray()

            full_jacobian.append(jacobian_selected)

            simulation_vector = res_simulation.take(select_algebraic)
            scaled_jac_i = (jacobian_selected.T * (1 / simulation_vector)).T
            scaled_yao_jacobian.append(scaled_jac_i)

            for index_algebraic in select_algebraic:
                variance.append(variance_i[index_algebraic])

        jac_array = np.concatenate(full_jacobian)
        measurement_variance_matrix = np.diagflat(variance)
        measurement_variance_array_inverted = np.linalg.inv(measurement_variance_matrix)

        fim_matrix = jac_array.T.dot(measurement_variance_array_inverted).dot(jac_array)
        parameter_covariance_matrix = np.linalg.inv(fim_matrix)  # type: ignore

        # jac_array_yao = np.concatenate(scaled_yao_jacobian) * np.sqrt(np.diag(parameter_covariance_matrix))
        jac_array_yao = np.concatenate(scaled_yao_jacobian) * selected_parameters_values

        result = {}
        result["jac"] = jac_array
        result["jac_yao"] = jac_array_yao
        result["cov_par"] = parameter_covariance_matrix
        result["cov_meas"] = measurement_variance_matrix
        result["fim"] = fim_matrix

        return result

    def parameter_identifiability_yao(
        self,
        parameters: dict[str, float],
        unfixed_params: list[str],
        threshold: float = 4e-2,
    ):
        """Do parameter ranking based on Yao 2003. Cut-off value taken from Yao 2003.
        Return ranked parameters in descending order, and divide them in identifiable and not"""
        self._setup_scaling(False)
        parameter_values_all: list[float] = []
        selected_parameters: list[float] = []
        unranked_parameters: list[str] = []

        for var_name in parameters.keys():
            if var_name in self.varlist_decision.keys():
                parameter_values_all.append(parameters[var_name])

        for var_name in parameters.keys():
            if var_name in unfixed_params:
                selected_parameters.append(parameters[var_name])
                unranked_parameters.append(var_name)

        results_sensitivity = self.calculate_sensitivity_and_fim(parameters)

        jacobian_yao = results_sensitivity["jac_yao"]

        XK = np.zeros(jacobian_yao.shape)

        parameters_ranked = []
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
            parameters_ranked.append(most_identifiable_parameter)

            # print(eucnorm[-1])
            if eucnorm[-1] < threshold:
                parameters_not_identifiable.append(most_identifiable_parameter)
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

        print(f"Ranked parameters: {parameters_ranked}")
        print(f"Estimable parameters: {parameters_identifiable}")
        print(f"Non identifiable parameters: {parameters_not_identifiable}")

        result = {}
        result["ranked"] = parameters_ranked
        result["estimable"] = parameters_identifiable
        result["fixed"] = parameters_not_identifiable

        return result

    def parameter_identifiability_eigenvalue(
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

        self._setup_scaling(False)
        if parameters_identifiable is None:
            parameters_identifiable = []

        if parameters_not_identifiable is None:
            parameters_not_identifiable = []

        results_sensitivity = self.calculate_sensitivity_and_fim(
            parameters, unfixed_params
        )

        fim_matrix = results_sensitivity["fim"]

        def eigsorted(cov):
            vals, vecs = np.linalg.eig(cov)
            # vals, vecs = np.linalg.eigh(cov)
            order = np.flip(vals.argsort())
            return vals[order], vecs[:, order]

        # vals, vecs = eigsorted(results_sensitivity["fim"])
        vals, vecs = eigsorted(fim_matrix)

        index_max = np.argmax(np.abs(vecs[:, -1]))
        current_parameter_name = unfixed_params.pop(index_max)
        # print(np.abs(vals[-1]))
        if np.abs(vals[-1]) > eigenvalue_threshold:
            parameters_identifiable.insert(0, current_parameter_name)
        else:
            parameters_not_identifiable.insert(0, current_parameter_name)

        if len(vals) == 1:
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

        else:
            self.parameter_identifiability_eigenvalue(
                parameters,
                unfixed_params,
                parameters_identifiable,
                parameters_not_identifiable,
                eigenvalue_threshold,
            )

    def parameter_analysis(self, parameters: dict[str, float]):
        import scipy.stats

        dof = len(self.list_input_varlist) - len(self.varlist_decision)
        num_par = len(self.varlist_decision)
        len_alg = len(self.model.varlist_algebraic)

        self._setup_scaling(False)

        # jac_ols = ca.jacobian(self._objective_ols()[0], decision_variables)

        # jac_ols_function = ca.Function(
        #     "jac_ols",
        #     [decision_variables],
        #     [jac_ols],
        #     ["x"],
        #     ["f"],
        # )

        selected_parameters = self.parameters_dict_to_list(parameters)
        S_squared, residuals, residuals_all = self.calculate_ols_value(parameters)
        residual_mean_square = S_squared / dof

        result_sens = self.calculate_sensitivity_and_fim(parameters)
        jac_array = result_sens["jac"]
        parameter_covariance_matrix = residual_mean_square * np.linalg.inv(jac_array.T.dot(jac_array))

        measurement_variance_matrix = result_sens["cov_meas"]

        residuals_sorted_in_columns = np.reshape(
            residuals_all, (len(self.list_simulators), len_alg)
        )
        residuals_sorted_in_columns = residuals_sorted_in_columns[:, residuals_sorted_in_columns.all(0)]
        S_squred_sorted_in_columns = np.sum(residuals_sorted_in_columns ** 2, axis=0)
        real_residual_mean_square = S_squred_sorted_in_columns / dof

        data_in_columns = np.reshape(
            self.array_data, (len(self.list_simulators), len_alg)
        )

        parameter_variance = np.diag(parameter_covariance_matrix)
        parameter_std = np.sqrt(parameter_variance).flatten()

        students_t_dist_95 = scipy.stats.t.ppf(0.975, dof)

        for par, var_value in zip(selected_parameters, parameter_std):
            print(f"{par} +- {var_value} |  ({var_value / par})")

        marginal_conf_interval_95 = (parameter_std * students_t_dist_95).T

        import matplotlib.pyplot as plt
        from matplotlib.axes import Axes
        from matplotlib.patches import Ellipse

        def eigsorted(cov):
            vals, vecs = np.linalg.eig(cov)
            # vals, vecs = np.linalg.eigh(cov)
            order = vals.argsort()[::-1]
            return vals[order], vecs[:, order]

        fig, axes = plt.subplots(ncols=num_par - 1, nrows=num_par - 1)
        if isinstance(axes, Axes) == 1:
            axes = [axes]

        fisher_f_dist_95 = scipy.stats.f.ppf(0.95, num_par, dof)

        comb = list(combinations(range(num_par), 2))
        par_names = list(self.varlist_decision.keys())

        title = ""
        for par_value, var_variance_i, name in zip(
            selected_parameters, parameter_std, par_names
        ):
            title = (
                title
                + f"{name}: {round(par_value,5)} ± {round(var_variance_i,5)} |  ({round((var_variance_i / par_value) * 100,1)}%)\n"
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
                width, height = (
                    2 * np.sqrt(num_par * residual_mean_square * fisher * vals)
                )[0]
                ellip = Ellipse(
                    xy=parameters_i, width=width, height=height, angle=theta, alpha=0.3
                )

                ax.add_artist(ellip)
            ax.relim()
            ax.autoscale()

            if i[0] == 0:
                ax.set_ylabel(f"{par_names[i[1]]}")

            if i[1] == len(par_names) - 1:
                ax.set_xlabel(f"{par_names[i[0]]}")

            ax.axvline(parameters_i[0] - marginal_conf_interval_95_i[0])
            ax.axvline(parameters_i[0] + marginal_conf_interval_95_i[0])
            ax.axhline(parameters_i[1] - marginal_conf_interval_95_i[1])
            ax.axhline(parameters_i[1] + marginal_conf_interval_95_i[1])

        plt.show()
        return residuals, S_squared, residual_mean_square


class ParameterEstimationNLE_control(ParameterEstimationNLE):
    def _setup_simulator(self, use_simulator_bounds):
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
                use_bounds=use_simulator_bounds,
            )
            list_simulators.append(simulator)

            self._setup_simulator_mapping(simulator)

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
            res_simulation = simulator.simulate_sym()

            if array_simulation is None:
                array_simulation = res_simulation["x"]
            else:
                array_simulation = ca.vertcat(array_simulation, res_simulation["x"])

        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (array_simulation - self.array_data) * self.array_data_mask
        error_controls = self.array_controls_casadi - self.array_controls
        objective = ca.sum1(error**2) + ca.sum1(error_controls**2)

        return objective

    def optimize(self, scale=False, objective_function="ols"):
        if objective_function == "wls":
            self._objective = self._objective_wls
        elif objective_function == "ols":
            self._objective = self._objective_ols

        res = self._optimize(scale)
        res["all"] = res["x"]
        res["x"] = res["all"][0 : self.num_parameters]
        res["p"] = res["all"][self.num_parameters :]
        return res
