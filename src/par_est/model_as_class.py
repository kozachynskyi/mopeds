import copy
from collections import OrderedDict
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt
from opcua import ua
from opcua.ua import NumericNodeId
from optipal.client import OptiPALClient


class Variable(object):

    """Docstring for Variable. """

    def __init__(self, name):
        """TODO: to be defined. """
        self.name = name
        self.casadi_var = ca.MX.sym(self.name)
        self.fixed = False
        self.opc_ua_id = None
        self.starting_value = None
        self.value = None
        self.guess = None
        self.lower_bound = None
        self.upper_bound = None


class State_variable(Variable):
    def __init__(self, name, starting_value=None, opc_ua_id=None):
        super().__init__(name)
        self.starting_value = starting_value
        self.value = Experimental_Data()
        self.opc_ua_id = opc_ua_id


class Parameter_variable(Variable):
    def __init__(self, name, value=None, lb=None, ub=None):
        super().__init__(name)
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub


class Control_variable(Variable):
    def __init__(self, name, value=None, lb=None, ub=None, opc_ua_id=None):
        super().__init__(name)
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub
        self.opc_ua_id = opc_ua_id


class VariableList(OrderedDict):
    def __init__(self):
        super().__init__()

    def add_variable(self, variable: Variable):
        """ TODO: add error handling if variable exists"""
        self.update({variable.name: variable})

    def get_casadi_var(self):
        casadi_vars = []
        for var in self.values():
            casadi_vars.append(var.casadi_var)
        return ca.vcat(casadi_vars)

    def get_data_opcua(self, time_start: datetime, time_stop: datetime):
        client = OptiPALClient("opc.tcp://admin@localhost:4840")  # type: OptiPALClient
        client.connect()
        try:
            ns_working = client.get_working_ns_idx()
            for var in self.values():
                values_opcua = []
                time_opcua = []
                if isinstance(var, State_variable):
                    sensor = client.get_node(NumericNodeId(var.opc_ua_id, ns_working))
                    process_value = client.get_child_simple(sensor, ["d:ProcessValue"])
                    results = process_value.read_raw_history(
                        time_start, time_stop, 1000
                    )
                    var.value = Experimental_Data()

                    for result in results:
                        if not time_opcua:
                            time_opcua.append(0.0)
                            time_zero = result.SourceTimestamp
                            var.starting_value = result.Value.Value
                        else:
                            time_from_ref = (
                                result.SourceTimestamp - time_zero
                            ).total_seconds()
                            time_opcua.append(time_from_ref)

                        values_opcua.append(result.Value.Value)

                    var.value.value = np.array(values_opcua)
                    var.value.time = np.array(time_opcua)
        finally:
            client.disconnect()

    def write_data_opcua(self, time_start: datetime):
        client = OptiPALClient("opc.tcp://admin@localhost:4840")  # type: OptiPALClient
        client.connect()
        try:
            time_zero = time_start
            ns_working = client.get_working_ns_idx()
            for var in self.values():
                if isinstance(var, State_variable):
                    sensor = client.get_node(NumericNodeId(var.opc_ua_id, ns_working))
                    process_value = client.get_child_simple(sensor, ["d:ProcessValue"])
                    for value, time in zip(var.value.value, var.value.time):
                        datavalue = ua.DataValue(value)
                        datavalue.SourceTimestamp = time_zero + timedelta(seconds=time)
                        process_value.set_attribute(ua.AttributeIds.Value, datavalue)
        finally:
            client.disconnect()


# Model Generation:
class Model(object):

    """Docstring for model. """

    def __init__(self, variable_list):
        """TODO: to be defined. """
        self.states = VariableList()
        self.variables = VariableList()
        self._all_variables = VariableList()
        self.equations = None

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, State_variable):
                    self.states.add_variable(State_variable(var.name))
                else:
                    self.variables.add_variable(Parameter_variable(var.name))
            else:
                raise (ValueError)

        self._all_variables.update(self.states)
        self._all_variables.update(self.variables)

    def add_equations(self, equations):
        if self.equations is None:
            self.equations = ca.vcat(equations)
        else:
            raise (NotImplementedError)


class Simulator(object):

    """Docstring for Simulator. """

    def __init__(self, model: Model, time_grid, variable_list: VariableList):
        """TODO: to be defined. """
        self.__input_variable_list = variable_list
        self.model = model
        self.tau = ca.MX.sym("tau")
        self.scaling = None
        self.time_grid = time_grid

        self.ode_system = {
            "x": self.model.states.get_casadi_var(),
            "p": ca.vertcat(self.model.variables.get_casadi_var()),
            "ode": self.model.equations,
        }
        self.ode_system_tau = {
            "x": self.model.states.get_casadi_var(),
            "p": ca.vertcat(self.tau, self.model.variables.get_casadi_var()),
            "ode": self.model.equations * self.tau,
        }

        self.integrator = ca.integrator(
            "integrator",
            "cvodes",
            self.ode_system,
            {"grid": self.time_grid, "output_t0": False, "print_stats": False},
        )

        # This integrator uses tau variable and is used in PE and OED
        self.integrator_tau = ca.integrator(
            "integrator",
            "cvodes",
            self.ode_system_tau,
            {
                "tf": 1,
                "output_t0": False,
                "print_stats": False,
                "linear_multistep_method": "adams",
            },
        )

        self._variables = []
        self._state_variables = VariableList()
        self._initial_states = []

        for var in self.__input_variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, State_variable):
                    self._initial_states.append(var.starting_value)
                    self._state_variables.add_variable(var)
                else:
                    if var.fixed:
                        self._variables.append(var.value)
                    else:
                        self._variables.append(var.casadi_var)

        self._variables = ca.vcat(self._variables)
        self._reset_scaling()

    def _reset_scaling(self):
        self.scaling = ca.DM.ones(self._variables.size())

    def simulate(self, derivatives=False):
        prev_time_step = 0
        res_states = []
        res_jacobian = []
        x_init = self._initial_states

        for time_step in self.time_grid[1:]:
            res_integration = self.integrator_tau(
                x0=x_init,
                p=ca.vertcat(
                    time_step - prev_time_step, self._variables * self.scaling
                ),
            )
            if derivatives:
                integrator_jac = self.integrator_tau.factory(
                    "I_fwd", ["x0", "p"], ["jac:xf:p"]
                )
                res_integration_jac = integrator_jac(
                    x0=x_init,
                    p=ca.vertcat(
                        time_step - prev_time_step, self._variables * self.scaling
                    ),
                )

            prev_time_step = time_step
            x_init = res_integration["xf"]

            if time_step == self.time_grid[1]:
                res_states = res_integration["xf"]
                if derivatives:
                    res_jacobian = res_integration_jac["jac_xf_p"]
            else:
                res_states = ca.horzcat(res_states, res_integration["xf"])
                if derivatives:
                    res_jacobian = ca.vertcat(
                        res_jacobian, res_integration_jac["jac_xf_p"]
                    )

        if derivatives:
            return res_states, res_jacobian
        else:
            return res_states

    def generate_exp_data(self):
        res_array = self.simulate()
        variables = VariableList()
        convert_to_numpy = False
        if isinstance(res_array, ca.DM):
            convert_to_numpy = True

        for count, var in enumerate(self._state_variables.values()):
            new_var = copy.deepcopy(var)
            new_var.value = Experimental_Data()
            if convert_to_numpy:
                new_var.value.time = self.time_grid
                new_var.value.value = res_array[count, :].toarray()
                new_var.value.value = np.insert(
                    new_var.value.value, 0, var.starting_value
                )
            else:
                new_var.value.time = self.time_grid[1:]
                new_var.value.value = res_array[count, :]

            variables.add_variable(new_var)

        return variables


class Experimental_Data(object):
    def __init__(self):
        self.time = None
        self.value = None

    def is_correct(self):
        if self.time is None or self.value is None:
            return False
        if self.time.size == self.value.size:
            return True
        else:
            False


class Optimizer(object):
    def __init__(self, model: Model, variable_list: VariableList):
        self.model = model
        self._input_variable_list = variable_list
        self.decision_var = VariableList()
        self.parameters = VariableList()
        self.controls = VariableList()
        self.states = VariableList()
        self.time_grid = None
        self.simulator = None  # type: Simulator

        self.guess = None
        self.lower_bound = None
        self.upper_bound = None
        self.scaling = None

        self.solver = None
        self.solver_settings = None

    def _setup_simulator(self):
        # Sets all variable lists of the class from __input_var_list
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

        self.scaling = 1
        self.simulator.scaling = None

    def _setup_scaling(self, scale=False):
        # Scaling should be done before setting a solver and solver settings
        # Sets scaling variables in optimizer and simulator
        if scale:
            self.simulator._reset_scaling()
            self.scaling = self.guess
            # for var in self.simulation._variables fails to iterate
            for count in range(self.simulator._variables.size()[0]):
                var = self.simulator._variables[count]
                if var.is_symbolic():
                    if var.name() in self.decision_var:
                        current_guess = self.decision_var[var.name()].guess
                        if current_guess == 0:
                            self.simulator.scaling[count] = 1
                        else:
                            self.simulator.scaling[count] = current_guess
                    else:
                        self.simulator.scaling[count] = 1
        else:
            self.scaling = 1
            self.simulator._reset_scaling()

    def _objective(self):
        # Returns a way to calculate and objective
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

        res_solver = self.solver(
            x0=self.guess / self.scaling,
            lbx=self.lower_bound / self.scaling,
            ubx=self.upper_bound / self.scaling,
        )
        print(res_solver["x"])
        print(res_solver["x"] * self.scaling)


class ParameterEstimation(Optimizer):
    def __init__(self, model: Model, variable_list: VariableList):
        super().__init__(model, variable_list)
        self._setup_simulator()
        self._setup_initialization()

    def _setup_simulator(self):
        self.time_grid = np.ndarray((1, 0))

        for var in self._input_variable_list.values():
            # Generates time_grid based on available exp data
            if isinstance(var, State_variable):
                self.states.add_variable(var)
                if var.value.is_correct():
                    self.time_grid = np.append(self.time_grid, var.value.time)
            elif isinstance(var, Parameter_variable):
                self.parameters.add_variable(var)
                if var.fixed is False:
                    self.decision_var.add_variable(var)
            elif isinstance(var, Control_variable):
                var.fixed = True
                self.controls.add_variable(var)

        self.time_grid = np.unique(self.time_grid)

        self.simulator = Simulator(
            self.model, self.time_grid, self._input_variable_list
        )

    def _objective(self):
        res_simulation = self.simulator.generate_exp_data()
        error = 0

        for var in self.states.values():
            if var.value.is_correct():
                for count_exp_point, time_point in enumerate(var.value.time[1:]):
                    # Looks up an index in a time_grid that has given time_point
                    res_index = np.nonzero(self.time_grid == time_point)
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

    def __init__(self, model: Model, variable_list: VariableList, time_grid):
        """TODO: to be defined. """
        super().__init__(model, variable_list)
        self.time_grid = time_grid
        self.parameter_values = []

        self._setup_simulator()
        self._setup_initialization()

    def _setup_simulator(self):
        for var in self._input_variable_list.values():
            if isinstance(var, Control_variable):
                if var.fixed is False:
                    self.decision_var.add_variable(var)
            elif isinstance(var, Parameter_variable):
                if var.fixed is False:
                    self.parameters.add_variable(var)
                    self.parameter_values.append(var.value)

        self.simulator = Simulator(
            self.model, self.time_grid, self._input_variable_list
        )

    def _sensitivity_matrix(self):
        res, res_jacobian = self.simulator.simulate(True)

        eval_jacobian = ca.Function(
            "eval_jacobian",
            [self.parameters.get_casadi_var(), self.decision_var.get_casadi_var()],
            [res_jacobian],
        )

        sensitivity_matrix = eval_jacobian(
            self.parameter_values, self.decision_var.get_casadi_var()
        )

        parameter_sensitivity_matrix = None

        for count in range(self.simulator._variables.size()[0]):
            var = self.simulator._variables[count]
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
        sc_states = [1, 100, 100, 100, 100]
        scale_states = np.diagflat(np.tile(sc_states, len(self.time_grid) - 1))
        scale_parameters = np.diagflat(self.parameter_values)

        sensitivity_scaled = scale_states @ (parameter_sensitivity_matrix @ scale_parameters)

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
        sensitivity_matrix = self._sensitivity_matrix()

        fim_matrix = sensitivity_matrix.T @ sensitivity_matrix
        error = ca.eig_symbolic(ca.inv(fim_matrix))
        # error = ca.trace(ca.inv(fim_matrix))

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


if __name__ == "__main__":

    e0_greek_nu_i1_r1 = -1.0
    e0_greek_nu_i1_r2 = 1.0
    e0_greek_nu_i2_r2 = -1.0
    e0_greek_nu_i3_r1 = 1.0
    e0_greek_nu_i1_r3 = -1.0
    e0_greek_nu_i4_r3 = 1.0
    e0_greek_Deltah_r1 = 0.0045
    e0_greek_Deltah_r2 = -0.0055
    e0_greek_Deltah_r3 = 0.0045
    e0_greek_rho = 800.0
    e0_A = 1.0
    e0_E_r1 = 96000.0
    e0_c_p = 3.5
    e0_E_r2 = 72000.0
    e0_E_r3 = 69000.0
    e0_F = 6.5e-4
    e0_R = 8.314
    e0_V = 1.0

    variable_list = VariableList()

    variable_list.add_variable(State_variable("e0_T", 273.0, 10))
    variable_list.add_variable(State_variable("e0_c_i1", 3.0, 20))
    variable_list.add_variable(State_variable("e0_c_i2", 10.0, 30))
    variable_list.add_variable(State_variable("e0_c_i3", 0.0, 40))
    variable_list.add_variable(State_variable("e0_c_i4", 0.0, 50))

    variable_list.add_variable(Parameter_variable("e0_k_pre_r1", 5000000.0))
    variable_list.add_variable(Parameter_variable("e0_k_pre_r2", 1.0e7))
    variable_list.add_variable(Parameter_variable("e0_k_pre_r3", 500000.0))
    variable_list.add_variable(Parameter_variable("e0_U", 1.4))

    variable_list.add_variable(Control_variable("e0_c_in_i1", 5.0))
    variable_list.add_variable(Control_variable("e0_c_in_i2", 10.0))
    variable_list.add_variable(Control_variable("e0_c_in_i3", 0.0))
    variable_list.add_variable(Control_variable("e0_c_in_i4", 0.0))
    variable_list.add_variable(Control_variable("e0_T_in", 373.0))
    variable_list.add_variable(Control_variable("e0_T_j", 373.0))

    m = Model(variable_list)

    # fmt: off
    tdot = (((((e0_F / e0_V) * ((m._all_variables["e0_T_in"].casadi_var - m._all_variables["e0_T"].casadi_var))) + (((m._all_variables["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (e0_c_p * e0_V))) * ((m._all_variables["e0_T_j"].casadi_var - m._all_variables["e0_T"].casadi_var)))) + (((-e0_greek_Deltah_r1) / (e0_greek_rho * e0_c_p)) * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-e0_greek_Deltah_r2) / (e0_greek_rho * e0_c_p)) * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (((-e0_greek_Deltah_r3) / (e0_greek_rho * e0_c_p)) * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c1dot = ((((e0_F / e0_V) * ((m._all_variables["e0_c_in_i1"].casadi_var - m._all_variables["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m._all_variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c2dot = ((e0_F / e0_V) * ((m._all_variables["e0_c_in_i2"].casadi_var - m._all_variables["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m._all_variables["e0_k_pre_r2"].casadi_var * (m._all_variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c3dot = ((e0_F / e0_V) * ((m._all_variables["e0_c_in_i3"].casadi_var - m._all_variables["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m._all_variables["e0_k_pre_r1"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    c4dot = ((e0_F / e0_V) * ((m._all_variables["e0_c_in_i4"].casadi_var - m._all_variables["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m._all_variables["e0_k_pre_r3"].casadi_var * (m._all_variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m._all_variables["e0_T"].casadi_var))))))
    # fmt: on

    m.add_equations([tdot, c1dot, c2dot, c3dot, c4dot])

    time_grid = np.linspace(10, 1000, 20)
    time_grid = np.insert(time_grid, 0, 0)

    var_list1 = copy.deepcopy(variable_list)
    for var in var_list1.values():
        var.fixed = True
    s = Simulator(m, time_grid, var_list1)
    res = s.simulate()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="\t")

    var_list_exp = s.generate_exp_data()

    var_list2 = copy.deepcopy(var_list1)
    for key, var in var_list_exp.items():
        var_list2[key] = var

    # var_list2["e0_T"].value = Experimental_Data()
    # var_list2["e0_c_i4"].value = Experimental_Data()

    var_list2["e0_k_pre_r1"].fixed = True
    var_list2["e0_k_pre_r1"].guess = 4000000.0
    var_list2["e0_k_pre_r1"].lower_bound = 4000000.0
    var_list2["e0_k_pre_r1"].upper_bound = 6000000.0

    var_list2["e0_k_pre_r2"].fixed = False
    var_list2["e0_k_pre_r2"].guess = 1.0e6
    var_list2["e0_k_pre_r2"].lower_bound = 1.0e6
    var_list2["e0_k_pre_r2"].upper_bound = 1.0e8

    var_list2["e0_k_pre_r3"].fixed = True
    var_list2["e0_k_pre_r3"].guess = 400000.0
    var_list2["e0_k_pre_r3"].lower_bound = 400000.0
    var_list2["e0_k_pre_r3"].upper_bound = 600000.0

    var_list2["e0_U"].fixed = False
    var_list2["e0_U"].guess = 1.1
    var_list2["e0_U"].lower_bound = 1.0
    var_list2["e0_U"].upper_bound = 3.0

    var_list2["e0_c_in_i1"].fixed = False
    var_list2["e0_c_in_i1"].guess = 5.0
    var_list2["e0_c_in_i1"].lower_bound = 4.0
    var_list2["e0_c_in_i1"].upper_bound = 6.0

    var_list2["e0_c_in_i2"].fixed = False
    var_list2["e0_c_in_i2"].guess = 10.0
    var_list2["e0_c_in_i2"].lower_bound = 9.0
    var_list2["e0_c_in_i2"].upper_bound = 11.0

    var_list2["e0_c_in_i3"].fixed = False
    var_list2["e0_c_in_i3"].guess = 0.0
    var_list2["e0_c_in_i3"].lower_bound = 0.0
    var_list2["e0_c_in_i3"].upper_bound = 0.0

    var_list2["e0_c_in_i4"].fixed = False
    var_list2["e0_c_in_i4"].guess = 0.0
    var_list2["e0_c_in_i4"].lower_bound = 0.0
    var_list2["e0_c_in_i4"].upper_bound = 0.0

    var_list2["e0_T_in"].fixed = False
    var_list2["e0_T_in"].guess = 373.0
    var_list2["e0_T_in"].lower_bound = 353.0
    var_list2["e0_T_in"].upper_bound = 393.0

    var_list2["e0_T_j"].fixed = False
    var_list2["e0_T_j"].guess = 373.0
    var_list2["e0_T_j"].lower_bound = 353.0
    var_list2["e0_T_j"].upper_bound = 393.0

    start_time = datetime(2018, 1, 1, 1, 0, 0, 0) + timedelta(days=1)
    end_time = start_time + timedelta(seconds=1000)
    var_list_oed = copy.deepcopy(var_list2)
    # var_list3 = copy.deepcopy(var_list2)
    # var_list_exp.write_data_opcua(start_time)
    # var_list3.get_data_opcua(start_time, end_time)
    pe = ParameterEstimation(m, var_list2)
    # pe.optimize()

    oed = OptimalExperimentalDesign(m, var_list_oed, time_grid)
    a = oed.get_fim_matrix()
    b = a[0].toarray()
    # b = ca.fabs(b)
    c = a[1].toarray()
    # turn_off_states = np.array([1, 1, 1, 1, 1])
    # sc_states = [1, 0.01, 0.01, 0.01, 0.01]
    # sc = np.diagflat(np.tile(turn_off_states, len(time_grid) - 1))
    # sc_full = np.diagflat(np.tile(turn_off_states / sc_states, len(time_grid) - 1))
    # # num_states = 5
    # # sc_full_params = np.tile(sc_params, ((len(time_grid) - 1) * num_states, 1)).T
    # sc_params = [5000000.0, 10000000.0, 500000.0, 1.4]
    # sc_full_params = np.diagflat(sc_params)
    # b_scaled = sc @ b
    # b_scaled_full = sc_full @ (b @ sc_full_params)
    # # sc = np.tile(sc, 2)
    fig = plt.figure()
    fig.add_subplot(151).imshow(b, cmap=cm.Greens_r)
    fig.add_subplot(152).imshow(ca.inv(c), cmap=cm.Greens_r)
    plt.show()
    oed.optimize()
    # pe.optimize(False)
