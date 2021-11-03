import copy
import logging
import casadi as ca
import numpy as np
from typing import List, Union, Dict
from scipy import linalg

from par_est import tools
from par_est import (
    VariableControl,
    VariableControlPiecewiseConstant,
    Model,
    VariableParameter,
    Simulator,
    VariableState,
    VariableList,
    VariableAlgebraic,
)


class Optimizer(object):
    def __init__(
        self,
        model: Model,
        variable_lists: List[VariableList],
        simulator_name,
        simulator_settings,
    ):
        if not isinstance(variable_lists, list):
            raise (Exception("Variable list should be nested of type list"))
        self.logger = logging.getLogger(__name__)
        self.model = model
        # Deepcopy is used to avoid manipulating input variable list
        self.list_input_varlist = copy.deepcopy(variable_lists)
        self.varlist_decision = VariableList()
        self.varlist_parameter = VariableList()
        self.varlist_control = VariableList()
        self.varlist_state = VariableList()
        self.varlist_algebraic = VariableList()
        self.simulator_name = simulator_name
        self.simulator_settings = simulator_settings
        self.list_simulators = []  # type: List[Simulator]

        self.guess = None
        self.lower_bound = None
        self.upper_bound = None
        self.scaling = 1

        self.solver = None
        self.solver_settings: Union[None, Dict] = None

    def _setup_simulator(self, *, use_idas_constraints):
        # Creates simulator
        raise (NotImplementedError)

    def optimize(self):
        # Runs optimization once
        raise (NotImplementedError)

    def _setup_initialization(self):
        """ Sets initials and bounds for optimizer, and as default no scaling. """
        guess = []
        lower_bound = []
        upper_bound = []

        for var in self.varlist_decision.values():
            if var.guess == 0:
                guess.append(1)
            else:
                guess.append(var.guess)
            lower_bound.append(var.lower_bound)
            upper_bound.append(var.upper_bound)

        self.guess = np.array(guess)
        self.lower_bound = np.array(lower_bound)
        self.upper_bound = np.array(upper_bound)
        self.logger.debug(
            f"Initialized:\nguess {self.guess}\nlower_bound {self.lower_bound}\nupper_bound {self.upper_bound}"
        )

    def _setup_scaling(self, scale=False):
        """Scaling should be done before setting a solver and solver settings.
        Sets scaling variables in optimizer and simulator.
        TODO: Whole loop can be replaced by simple np.where, isn't it?
        """
        if scale:
            self.scaling = np.where(self.guess == 0, 1, self.guess)
            for simulator in self.list_simulators:
                simulator._reset_scaling()
                # for var in self.simulation._variables fails to iterate
                for count in range(simulator._variables[0].size()[0]):
                    var = simulator._variables[0][count]
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

    def _objective(self):
        """ Returns a way to calculate and objective. Dependent on optimization type. """
        raise (NotImplementedError)

    def _optimize(self, scale):
        """ Scaling should be done before setting a solver and solver settings. """
        self._setup_scaling(scale)

        self.solver = ca.nlpsol(
            "solver",
            self.solver_name,
            {"x": self.varlist_decision.get_casadi_variables(), "f": self._objective()},
            self.solver_settings,
        )

        # Scaling of negative numbers requires a switch bounds
        lb_scaled = self.lower_bound / self.scaling
        ub_scaled = self.upper_bound / self.scaling

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
        return res_solver

    def optimize_multistart(self, num_initials, scale=True, max_iterations=20):
        hammersley_seeds = np.array(list(zip(self.lower_bound, self.upper_bound)))

        list_startpoint = tools.make_startpoints(hammersley_seeds, num_initials)

        result = []

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
        for index, guess in enumerate(list_startpoint):
            self.guess = guess
            print(f"Optimization number {index} started")
            res = self.optimize(scale)
            print(f"Objective: {res['f']}")
            result.append(res)

        self.solver_settings = initial_settings
        return result


class ParameterEstimation(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_list: List[VariableList],
        simulator_name="idas",
        simulator_settings=None,
        *,
        reinitialize_algebraic=False,
        use_idas_constraints=False,
        use_algebraic_vars=False,
        reinitialize_algebraic_experimental=False,
    ):
        super().__init__(
            model,
            variable_list,
            simulator_name,
            simulator_settings,
        )

        if use_algebraic_vars:
            self._objective = self._objective_alg
        else:
            self._objective = self._objective_state

        # This attribute is used while calculating Objective, and is either 1 or self.experiments_weights
        self.experiments_scale = 1
        self.experiments_weights: List[np.ndarray] = []
        self.array_data: List[np.ndarray] = []
        self.array_data_mask: List[np.ndarray] = []
        self.inverted_variances: List[np.ndarray] = []

        self._setup_simulator(use_idas_constraints=use_idas_constraints, use_algebraic_vars=use_algebraic_vars)
        self.logger.debug(
            "Created Optimizer object: \n Data Shape {} \n Desicion Variables {}".format(
                self.array_data.shape, self.varlist_decision.get_variable_name()  # type: ignore
            )
        )
        self._setup_initialization()

        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {"hessian_approximation": "limited-memory", "max_iter": 300},
        }

        if reinitialize_algebraic:
            for sim in self.list_simulators:
                sim.calculate_algebraic_initials(apply_intials=True)

        if reinitialize_algebraic_experimental:
            for sim in self.list_simulators:
                sim.calculate_algebraic_initials(apply_intials=True, analyze=True, experimental=True)

    def _setup_simulator(self, *, use_idas_constraints, use_algebraic_vars):
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

        list_timegrid_length = []
        size_simulation_output = []

        for varlist_input in self.list_input_varlist:
            # Create a time_grid, that "stops" at every experimental data, for every state variable
            time_grid = np.ndarray((1, 0))
            for var in varlist_input.values():
                if isinstance(var, VariableState) or (isinstance(var, VariableAlgebraic) and use_algebraic_vars):
                    time_grid = np.append(time_grid, var.value.time)
                elif isinstance(var, VariableControl):
                    var.fixed = True
                    if isinstance(var, VariableControlPiecewiseConstant):
                        var.fixed = True
                        time_grid = np.append(time_grid, var.time)

            time_grid = np.unique(time_grid)
            list_timegrid_length.append(float(len(time_grid)))

            self.list_simulators.append(
                Simulator(
                    self.model,
                    time_grid,
                    varlist_input,
                    self.simulator_name,
                    self.simulator_settings,
                    use_idas_constraints=use_idas_constraints,
                )
            )

            # Generate an array (experiment_data_varlist) with Experimental data with the same dimensions as simulation results.
            experiment_data_varlist = []
            experiment_data_mask_varlist = []

            if use_algebraic_vars:
                variable_name_list = list([*self.model.varlist_state.keys(), *self.model.varlist_algebraic.keys()])
            else:
                variable_name_list = list(self.model.varlist_state.keys())

            for var_name in variable_name_list:
                var = varlist_input[var_name]
                time_grid_var = np.array(var.value.time)
                # if simulated point has data - set element to True
                experiment_data_mask_var = (
                    1.0 * np.isin(time_grid, time_grid_var)[1:]
                )
                experiment_data_var_real = np.array(var.value.value)[1:]
                # array that would be filled with Experimental data where data_mask is 1
                experiment_data_var_extended = experiment_data_mask_var.copy()

                # data_var is being redimensioned to the output of simulation
                counter = 0
                for timegrid_index, trigger in enumerate(experiment_data_mask_var):
                    if trigger == 1:
                        experiment_data_var_extended[
                            timegrid_index
                        ] = experiment_data_var_real[counter]
                        counter = counter + 1
                experiment_data_varlist.append(experiment_data_var_extended)
                experiment_data_mask_varlist.append(experiment_data_mask_var)

            # Stack data from separate variables and flatten columnwise
            experiment_data_varlist = np.column_stack(experiment_data_varlist).flatten()
            experiment_data_mask_varlist = np.column_stack(
                experiment_data_mask_varlist
            ).flatten()
            self.array_data.append(experiment_data_varlist)
            self.array_data_mask.append(experiment_data_mask_varlist)

            """ Generate arrays with inverted variances and experiments weightning.
            Varainces are used for generation of weighted least squares optimization
            problem. Experiments weightning is used in order to give same weight to
            separate experiments: if one experiment has twice as many experimental
            points, their error is multiplied by 0.5.
            """
            inverted_variances_varlist = []
            for var_name in variable_name_list:
                var = varlist_input[var_name]
                inverted_variances_varlist.append(
                    1.0 / (np.full(len(time_grid) - 1, var.variance))
                )
            inverted_variances_varlist = np.column_stack(
                inverted_variances_varlist
            ).flatten()
            self.inverted_variances.append(inverted_variances_varlist)
            size_simulation_output.append(len(inverted_variances_varlist))

        max_time_grid = max(list_timegrid_length)
        for time_grid_length, size_simulation in zip(
            list_timegrid_length, size_simulation_output
        ):
            self.experiments_weights.append(
                np.full(size_simulation, max_time_grid / time_grid_length)
            )

        self.array_data = np.concatenate(self.array_data)
        self.array_data_mask = np.concatenate(self.array_data_mask)
        self.inverted_variances = np.concatenate(self.inverted_variances)
        self.experiments_weights = np.concatenate(self.experiments_weights)

    def _objective_state(self):
        array_simulation = None

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate()
            if array_simulation is None:
                array_simulation = res_simulation["xf"][:]
            else:
                array_simulation = ca.vertcat(array_simulation, res_simulation["xf"][:])

        # multiply by self.array_data_mask needed to ignore elements were error experimental data is zero
        error = (array_simulation - self.array_data) * self.array_data_mask
        objective = ca.sum1(
            self.experiments_scale * self.inverted_variances * (error ** 2)
        )

        return objective

    def _objective_alg(self):
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
        error = (array_simulation - self.array_data) * self.array_data_mask
        objective = ca.sum1(
            self.experiments_scale * self.inverted_variances * (error ** 2)
        )

        return objective

    def optimize(self, scale=True, *, scale_experiments=False):
        """Solves optimization problem. Scaling decreases amount of iterations,
        and should be used as a first option.
        """
        if scale_experiments:
            self.experiments_scale = self.experiments_weights
        else:
            self.experiments_scale = 1

        return self._optimize(scale)


class OptimalExperimentalDesign(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_list: List[VariableList],
        time_grid,
        simulator_name="idas",
        simulator_settings=None,
        *,
        reinitialize_algebraic=False,
    ):
        super().__init__(model, variable_list, simulator_name, simulator_settings)
        self.time_grid_original = time_grid
        self.parameter_values: List[float] = []
        self.select_independent: List[int] = []
        self.inverted_variances: List[float] = []

        self._setup_simulator()
        self._setup_initialization()
        # Is not the same as self.original_time_grid if VariableControlPiecewiseConstant adds additional time_stamps
        self.time_grid_modified = self.list_simulators[0].time_grid

        self.solver_name = "ipopt"
        self.solver_settings = {
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

    def _setup_simulator(self, *, use_idas_constraints=False):
        """Initializes simulator class. Parameter variables are fixed, and an index of an unfixed
        parameter is saved in self.select_independent list.
        This list is used during the calculation of the objective, to ignore jacobian of fixed parameters.
        self.index_all_states is used additionaly to self.select_independent list to get required jacobian.
        """

        for var in self.list_input_varlist[0].values():
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
                    self.parameter_values.append(var.value)
                    var.fixed = True
                    # index has + 1 to account for tau variable, that is a parameter of a simulation.
                    index = list(self.model.varlist_independent).index(var.name) + 1
                    self.select_independent.append(index)
            elif isinstance(var, VariableState):
                self.inverted_variances.append(1 / var.variance)

        self.inverted_variances = np.array(self.inverted_variances)
        self.index_all_states = list(range(len(self.model.varlist_state)))

        self.list_simulators.append(
            Simulator(
                self.model,
                self.time_grid_original,
                self.list_input_varlist[0],
                self.simulator_name,
                self.simulator_settings,
            )
        )

    def _objective(self, analyze=False, values=None):
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
            return error, covariance_full, jacobian_full, covariance_all, objective_all

        return error

    def optimize(self, scale=False):
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
        """ Taken from Erik/Diana Subset0. Many questions arrise about how it works. """
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

        Jsqr = jacobian ** 2

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
