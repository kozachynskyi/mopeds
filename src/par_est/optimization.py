import copy
import casadi as ca
import numpy as np

from par_est import (
    Control_variable,
    Model,
    Parameter_variable,
    Simulator,
    State_variable,
    VariableList
)


class Optimizer(object):
    def __init__(self, model: Model, variable_lists: [VariableList]):
        if not isinstance(variable_lists, list):
            raise (Exception("Variable list should be nested of type list"))
        self.model = model
        # Deepcopy is used to avoid manipulating input variable list
        self.variable_lists = copy.deepcopy(variable_lists)
        self.decision_var = VariableList()
        self.parameters = VariableList()
        self.controls = VariableList()
        self.states = VariableList()
        self.simulators = []  # type: [Simulator]

        self.guess = None
        self.lower_bound = None
        self.upper_bound = None
        self.scaling = 1

        self.solver = None
        self.solver_settings = None

    def _setup_simulator(self):
        # Sets all variable lists of the class from __input_var_list. Optimizer dependent
        raise (NotImplementedError)

    def _setup_initialization(self):
        # Sets initials and bounds for optimizer, and as default no scaling
        guess = []
        lower_bound = []
        upper_bound = []

        for var in self.decision_var.values():
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
            for simulator in self.simulators:
                simulator._reset_scaling()
                # for var in self.simulation._variables fails to iterate
                for count in range(simulator._variables.size()[0]):
                    var = simulator._variables[count]
                    if var.is_symbolic():
                        if var.name() in self.decision_var:
                            current_guess = self.decision_var[var.name()].guess
                            if current_guess == 0:
                                simulator.scaling[count] = 1
                            else:
                                simulator.scaling[count] = current_guess
                        else:
                            simulator.scaling[count] = 1
        else:
            self.scaling = 1
            for simulator in self.simulators:
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
            {"x": self.decision_var.get_casadi_var(), "f": self._objective()},
            self.solver_settings,
        )

        # Scaling of negative numbers requires switch bounds
        lb_scaled = (self.lower_bound / self.scaling)
        ub_scaled = (self.upper_bound / self.scaling)
        for index, (lb, ub) in enumerate(zip(lb_scaled, ub_scaled)):
            if lb > ub:
                lb_scaled[index] = ub
                ub_scaled[index] = lb


        res_solver = self.solver(
            x0=self.guess / self.scaling,
            lbx=lb_scaled,
            ubx=ub_scaled,
        )

        print(res_solver["x"])
        print(res_solver["x"] * self.scaling)


class ParameterEstimation(Optimizer):
    def __init__(self, model: Model, variable_list: VariableList):
        super().__init__(model, variable_list)
        self._setup_simulator()
        self._setup_initialization()

    def _setup_simulator(self):
        for var in self.variable_lists[0].values():
            if isinstance(var, State_variable):
                self.states.add_variable(var)
            elif isinstance(var, Parameter_variable):
                self.parameters.add_variable(var)
                if var.fixed is False:
                    self.decision_var.add_variable(var)
            elif isinstance(var, Control_variable):
                self.controls.add_variable(var)

        for variable_list in self.variable_lists:
            time_grid = np.ndarray((1, 0))

            for var in variable_list.values():
                # Generates time_grid based on available exp data
                if isinstance(var, State_variable):
                    if var.value.is_correct():
                        time_grid = np.append(time_grid, var.value.time)
                elif isinstance(var, Control_variable):
                    var.fixed = True

            time_grid = np.unique(time_grid)
            self.simulators.append(Simulator(self.model, time_grid, variable_list))

    def _objective(self):
        error = 0

        for simulator in self.simulators:
            res_simulation = simulator.generate_exp_data()

            for var in self.states.values():
                if var.value.is_correct():
                    for count_exp_point, time_point in enumerate(var.value.time[1:]):
                        # Looks up an index in a time_grid that has given time_point
                        res_index = np.nonzero(simulator.time_grid == time_point)
                        res_index = res_index[0][0]

                        calculated_value = res_simulation[var.name].value.value[
                            res_index - 1
                        ]
                        experimental_value = var.value.value[count_exp_point + 1]
                        error_at_timepoint = (
                            0.5 * (calculated_value - experimental_value)
                        ) ** 2
                        error = error + error_at_timepoint

        return error

    def optimize(self, scale=True):
        # Scaling decreases amount of iterations, but ipopt fails gradient check at big amount of timestamps
        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {"max_iter": 300, "derivative_test": "first-order"},
        }

        self._optimize(scale)


class OptimalExperimentalDesign(Optimizer):

    """Docstring for OptimalExperimentalDesign. """

    def __init__(self, model: Model, variable_list: [VariableList], time_grid):
        """TODO: to be defined. """
        super().__init__(model, variable_list)
        self.time_grid = time_grid
        self.parameter_values = []

        self._setup_simulator()
        self._setup_initialization()

    def _setup_simulator(self):
        for var in self.variable_lists[0].values():
            if isinstance(var, Control_variable):
                if var.fixed is False:
                    self.decision_var.add_variable(var)
            elif isinstance(var, Parameter_variable):
                if var.fixed is False:
                    self.parameters.add_variable(var)
                    self.parameter_values.append(var.value)

        self.simulators.append(
            Simulator(self.model, self.time_grid, self.variable_lists[0])
        )

    def _sensitivity_matrix(self):
        res, res_jacobian = self.simulators[0].simulate(True)

        eval_jacobian = ca.Function(
            "eval_jacobian",
            [self.parameters.get_casadi_var(), self.decision_var.get_casadi_var()],
            [res_jacobian],
        )

        sensitivity_matrix = eval_jacobian(
            self.parameter_values, self.decision_var.get_casadi_var()
        )

        parameter_sensitivity_matrix = None

        for count in range(self.simulators[0]._variables.size()[0]):
            var = self.simulators[0]._variables[count]
            if var.is_symbolic():
                if var.name() in self.parameters:
                    # Count + 1 because first variable in sensitivity matrix is tau
                    if parameter_sensitivity_matrix is None:
                        parameter_sensitivity_matrix = sensitivity_matrix[:, count + 1]
                    else:
                        parameter_sensitivity_matrix = ca.horzcat(
                            parameter_sensitivity_matrix,
                            sensitivity_matrix[:, count + 1],
                        )

        # TODO use variabnces
        sc_states = np.ones(len(self.simulators[0]._state_variables)).tolist()
        scale_states = np.diagflat(np.tile(sc_states, len(self.time_grid) - 1))
        scale_parameters = np.diagflat(self.parameter_values)

        sensitivity_scaled = scale_states @ (
            parameter_sensitivity_matrix @ scale_parameters
        )

        return sensitivity_scaled

    def get_fim_matrix(self):
        self._setup_scaling(False)
        sensitivity_matrix = self._sensitivity_matrix()
        evaluate = ca.Function(
            "eval_fim", [self.decision_var.get_casadi_var()], [sensitivity_matrix]
        )
        sensitivity_matrix = evaluate(self.guess)
        fim_matrix = sensitivity_matrix.T @ sensitivity_matrix
        return sensitivity_matrix, fim_matrix

    def _objective(self):
        # Only trace is programmed in Casadi
        sensitivity_matrix = self._sensitivity_matrix()

        fim_matrix = sensitivity_matrix.T @ sensitivity_matrix
        error = ca.trace(ca.inv(fim_matrix))
        # error = ca.eig_symbolic(ca.inv(fim_matrix))

        return error

    def optimize(self, scale=True):
        # Scaling decreases amount of iterations, but ipopt fails gradient check at big amount of timestamps
        self.solver_name = "ipopt"
        self.solver_settings = {
            "verbose": False,
            "ipopt": {
                "hessian_approximation": "limited-memory",
                "max_iter": 20,
                "derivative_test": "first-order",
            },
        }

        self._optimize(scale)
