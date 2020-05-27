import copy
import logging
import casadi as ca
import numpy as np
from typing import List
from scipy.sparse import csc_matrix, vstack
from scipy import linalg

from par_est import (
    VariableControl,
    Model,
    VariableParameter,
    Simulator,
    VariableState,
    VariableList,
)


class Optimizer(object):
    def __init__(
        self,
        model: Model,
        variable_lists: [VariableList],
        integrator_name,
        integrator_settings,
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
        self.integrator_name = integrator_name
        self.integrator_settings = integrator_settings
        self.list_simulators = []  # type: List[Simulator]

        self.guess = None
        self.lower_bound = None
        self.upper_bound = None
        self.scaling = 1

        self.solver = None
        self.solver_settings = None

    def _setup_simulator(self):
        # Creates simulator
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

    def _setup_scaling(self, scale=False):
        """ Scaling should be done before setting a solver and solver settings.
        Sets scaling variables in optimizer and simulator.
        """
        if scale:
            self.scaling = self.guess
            for simulator in self.list_simulators:
                simulator._reset_scaling()
                # for var in self.simulation._variables fails to iterate
                for count in range(simulator._variables.size()[0]):
                    var = simulator._variables[count]
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
            {"x": self.varlist_decision.get_casadi_var(), "f": self._objective()},
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
            x0=self.guess / self.scaling, lbx=lb_scaled, ubx=ub_scaled,
        )

        print(res_solver["x"])
        print(res_solver["x"] * self.scaling)

        return res_solver


class ParameterEstimation(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_list: VariableList,
        integrator_name="idas",
        integrator_settings=None,
    ):
        super().__init__(model, variable_list, integrator_name, integrator_settings)
        self._setup_simulator()
        self.logger.debug(
            "Created Optimizer object: \n Data Shape {} \n Desicion Variables {}".format(
                self.array_data.shape, self.varlist_decision.get_variable_name()
            )
        )
        self._setup_initialization()

    def _setup_simulator(self):
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

        self.array_data = None

        for varlist_input in self.list_input_varlist:
            time_grid = np.ndarray((1, 0))

            for var in varlist_input.values():
                # Generates time_grid based on available exp data
                if isinstance(var, VariableState):
                    time_grid = np.append(time_grid, var.value.time)
                    var.starting_value = var.value.value[0]
                elif isinstance(var, VariableControl):
                    var.fixed = True

            time_grid = np.unique(time_grid)

            self.list_simulators.append(
                Simulator(
                    self.model,
                    time_grid,
                    varlist_input,
                    self.integrator_name,
                    self.integrator_settings,
                )
            )

            for var in varlist_input.values():
                if isinstance(var, VariableState):
                    time_grid_var = np.array(var.value.time)
                    data_mask_var = np.isin(time_grid, time_grid_var)[1:]
                    data_var = np.array(var.value.value)[1:]

                    sparsity_pattern = csc_matrix(data_mask_var.astype(int))
                    sparsity_pattern.data = data_var

                    if self.array_data is None:
                        self.array_data = sparsity_pattern
                    else:
                        self.array_data = vstack([self.array_data, sparsity_pattern])

            self.array_data_sparcity = ca.DM(self.array_data).sparsity()

    def _objective(self):
        array_simulation = None

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate()

            if array_simulation is None:
                array_simulation = res_simulation["xf"]
            else:
                array_simulation = ca.vertcat(array_simulation, res_simulation["xf"])

        error = ca.sumsqr(
            array_simulation.get(False, self.array_data_sparcity) - self.array_data
        )

        return error

    def optimize(self, scale=True):
        """ Solves optimization problem. Scaling decreases amount of iterations,
        and should be used as a first option.
        """
        self.solver_name = "ipopt"
        if self.solver_settings is None:
            self.solver_settings = {
                "verbose": False,
                "ipopt": {"hessian_approximation": "limited-memory", "max_iter": 300},
            }

        return self._optimize(scale)


class OptimalExperimentalDesign(Optimizer):
    def __init__(
        self,
        model: Model,
        variable_list: [VariableList],
        time_grid,
        integrator_name="idas",
        integrator_settings=None,
    ):
        super().__init__(model, variable_list, integrator_name, integrator_settings)
        self.time_grid = time_grid
        self.parameter_values = []
        self.select_independent = []

        self._setup_simulator()
        self._setup_initialization()

    def _setup_simulator(self):
        """ Initializes simulator class. Parameter variables are fixed, and an index of an unfixed
        parameter is saved in self.select_independent list.
        This list is used during the calculation of the objective, to ignore jacobian of fixed parameters.
        self.index_all_states is used additionaly to self.select_independent list to get required jacobian.
        """

        for var in self.list_input_varlist[0].values():
            if isinstance(var, VariableControl):
                if var.fixed is False:
                    self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableParameter):
                if var.fixed is False:
                    self.varlist_parameter.add_variable(var)
                    self.parameter_values.append(var.value)
                    var.fixed = True
                    # index has + 1 to account for tau variable, that is a parameter of a simulation.
                    index = list(self.model.varlist_independent).index(var.name) + 1
                    self.select_independent.append(index)

        self.index_all_states = list(range(len(self.model.varlist_state)))

        self.list_simulators.append(
            Simulator(
                self.model,
                self.time_grid,
                self.list_input_varlist[0],
                self.integrator_name,
                self.integrator_settings,
            )
        )

    def _objective(self, analyze=False, values=None):
        """ Calculates an A OED criteria, beacuse casadi cannot do other.
        "analyze" Flag and values are used for debugging.

        Args:
            analyze: used for debugging, calculates objective and covariance.
            values: this values are used as desicionb variables for calculating of the objective.
        """
        covariance_full = None
        jacobian_full = None
        # -1 ignores time point zero in self.time_grid
        num_time = len(self.time_grid) - 1
        # +1 account for tau variable in Simulator class
        num_param = len(self.model.varlist_independent) + 1
        num_state = len(self.model.varlist_state)

        result_simulation = self.list_simulators[0].simulate(True)
        result_jacobian = result_simulation["jac_xf_p"]

        # Used only for debugging
        if analyze is True:
            self._setup_scaling(False)
            evaluate = ca.Function(
                "eval_fim", [self.varlist_decision.get_casadi_var()], [result_jacobian],
            )
            if values is None:
                result_jacobian = evaluate(self.guess)
            else:
                result_jacobian = evaluate(values)

        # Simulation returns jacobian that has to be split, to get jac at each time point.
        # list_jacobian_at_timepoint contains a list of that jacobians.
        split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param
        list_jacobian_at_timepoint = ca.horzsplit(result_jacobian, split_vector)

        # For now the measurement covariance is fixed, but it should be taken as input.
        covariance_measurement = np.eye(num_state)

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
                jacobian_selected.T @ covariance_measurement @ jacobian_selected
            )

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
            return error, covariance_full, jacobian_full

        return error

    def optimize(self, scale=True):
        """ Run optimization.

        Args:
            scale: scaling should be used as default, allows for faster convergence
        """
        # Scaling decreases amount of iterations, but ipopt fails gradient check at big amount of timestamps
        self.solver_name = "ipopt"
        if self.solver_settings is None:
            self.solver_settings = {
                "verbose": False,
                # "monitor": ["nlp_grad_f", "nlp_f"],
                "ipopt": {
                    "hessian_approximation": "limited-memory",
                    "max_iter": 100,
                    # "print_level": 6,
                },
            }

        return self._optimize(scale)

    def identifiability_analysis(self, reset_self=False):
        """ Taken from Erik/Diana Subset0. Many questions arrise about how it works. """
        _, _, jacobian = self._objective(True)
        cond_threshold = 1000
        colin_threshold = 15

        states, parameters = jacobian.shape

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
        minVal = np.min(values)

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

        SummationSv = np.sum(values)
        NeglectSv = np.sum(values[rank:] / SummationSv)

        # Determination of permutation matrix P by construction a RRQR of S
        Q, R, P = linalg.qr(jacobian, pivoting=True)

        # Use this for ranking
        IdentifOrd = P
        # Why
        # IdentifOrd.append(P)

        # Condition Number of J (Golub, 1996 & Hansen, 1998)
        CondN = maxVal / minVal

        # Still to understand!!!
        # collinIndex.Colind
        # collinIndex.Sub

        Jsqr = jacobian ** 2

        ParNorm = np.empty(jacobian.shape[1])

        for i in range(0, jacobian.shape[1]):
            ParNorm[i] = np.sqrt(np.sum(Jsqr[:, i] / jacobian.shape[0]))

        SensitivityOrder = np.argsort(ParNorm)[::-1]
        B = np.sort(ParNorm)[::-1]

        SensityOrd = B, SensitivityOrder

        TotalVariance = np.sum(values ** 2)
        NeglectVariance = np.sum(values[rank:] ** 2 / TotalVariance)

        unfix_parameters = []
        for index in range(rank):
            parameter_name = list(self.varlist_parameter.values())[
                IdentifOrd[index]
            ].name
            unfix_parameters.append(parameter_name)

        if reset_self:
            new_varlist = copy.deepcopy(self.list_input_varlist)

            for var in new_varlist[0].values():
                if isinstance(var, VariableParameter):
                    if var.name in unfix_parameters:
                        var.fixed = False
                    else:
                        var.fixed = True

            self.__init__(self.model, new_varlist, self.time_grid)

        return unfix_parameters
