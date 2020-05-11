import copy
import logging
import casadi as ca
import numpy as np
from typing import List
from scipy.sparse import csc_matrix, vstack

from par_est import (
    VariableControl,
    Model,
    VariableParameter,
    Simulator,
    VariableState,
    VariableList,
)


class Optimizer(object):
    def __init__(self, model: Model, variable_lists: [VariableList]):
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
        # Sets initials and bounds for optimizer, and as default no scaling
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
        # Scaling should be done before setting a solver and solver settings
        # Sets scaling variables in optimizer and simulator
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
        # Returns a way to calculate and objective. Dependent on optimization type
        raise (NotImplementedError)

    def _optimize(self, scale):
        # Scaling should be done before setting a solver and solver settings
        self._setup_scaling(scale)

        self.solver = ca.nlpsol(
            "solver",
            self.solver_name,
            {"x": self.varlist_decision.get_casadi_var(), "f": self._objective()},
            self.solver_settings,
        )

        # Scaling of negative numbers requires switch bounds
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
    def __init__(self, model: Model, variable_list: VariableList):
        super().__init__(model, variable_list)
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

            self.list_simulators.append(Simulator(self.model, time_grid, varlist_input))

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
        # Scaling decreases amount of iterations, but ipopt fails gradient check at big amount of timestamps
        self.solver_name = "ipopt"
        if self.solver_settings is None:
            self.solver_settings = {
                "verbose": False,
                "ipopt": {"hessian_approximation": "limited-memory", "max_iter": 300},
            }

        return self._optimize(scale)


class OptimalExperimentalDesign(Optimizer):

    """Docstring for OptimalExperimentalDesign. """

    def __init__(self, model: Model, variable_list: [VariableList], time_grid):
        """TODO: to be defined. """
        super().__init__(model, variable_list)
        self.time_grid = time_grid
        self.parameter_values = []
        self.ignore_independent = []
        self.select_independent = []

        self._setup_simulator()
        self._setup_initialization()

    def _setup_simulator(self):

        for var in self.list_input_varlist[0].values():
            if isinstance(var, VariableControl):
                index = list(self.model.varlist_independent).index(var.name) + 1
                self.ignore_independent.append(index)
                if var.fixed is False:
                    self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableParameter):
                index = list(self.model.varlist_independent).index(var.name) + 1
                if var.fixed is False:
                    self.varlist_parameter.add_variable(var)
                    self.parameter_values.append(var.value)
                    var.fixed = True
                    print(var.name)
                    self.select_independent.append(index)
                else:
                    index = list(self.model.varlist_independent).index(var.name) + 1
                    self.ignore_independent.append(index)

        # First parameter is always tau
        self.ignore_independent.append(0)
        self.ignore_independent.sort()

        self.index_all_parameters = list(range(len(self.model.varlist_independent) + 1))
        self.index_all_states = list(range(len(self.model.varlist_state)))

        self.list_simulators.append(
            Simulator(self.model, self.time_grid, self.list_input_varlist[0])
        )

    def _objective(self, analyze=False, values=None):
        # Change of basis https://www.youtube.com/watch?v=P2LTAUO1TdA&list=PLZHQObOWTQDPD3MizzM2xVFitgF8hE_ab&index=13
        # Only trace criterium is programmed in Casadi
        covariance_full = None
        num_time = len(self.time_grid) - 1
        num_param = len(self.model.varlist_independent) + 1
        num_state = len(self.model.varlist_state)
        result_simulation = self.list_simulators[0].simulate(True)
        res_jacobian = result_simulation["jac_xf_p"]

        if analyze is True:
            self._setup_scaling(False)
            evaluate = ca.Function(
                "eval_fim", [self.varlist_decision.get_casadi_var()], [res_jacobian]
            )
            if values is None:
                res_jacobian = evaluate(self.guess)
            else:
                res_jacobian = evaluate(values)

        split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param
        list_jacobian_at_timepoint = ca.horzsplit(res_jacobian, split_vector)

        covariance_measurement = np.eye(num_state)

        parameter_scaling = ca.repmat(ca.DM(self.parameter_values).T, num_state, 1)

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

        error = ca.trace(ca.inv(covariance_full))

        if analyze:
            return error, covariance_full

        return error

    def optimize(self, scale=True):
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
