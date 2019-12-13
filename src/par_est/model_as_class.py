import copy
from collections import OrderedDict

import casadi as ca
import numpy as np


class Variable(object):

    """Docstring for Variable. """

    def __init__(self, name):
        """TODO: to be defined. """
        self.name = name
        self.fixed = False
        self.casadi_var = ca.MX.sym(self.name)
        self.starting_value = None
        self.value = None
        self.guess = None
        self.lower_bound = None
        self.upper_bound = None


class State_variable(Variable):
    def __init__(self, name, starting_value=None):
        super().__init__(name)
        self.starting_value = starting_value


class Parameter_variable(Variable):
    def __init__(self, name, value=None, lb=None, ub=None):
        super().__init__(name)
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub


class Control_variable(Variable):
    def __init__(self, name, value=None, lb=None, ub=None):
        super().__init__(name)
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub


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

    def create_simulation(self):
        variable_list = self.__input_variable_list
        self._variables = []
        self._unfixed_variables = VariableList()
        self._state_variables = VariableList()
        prev_time_step = 0
        res_states = []
        x_init = []

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, State_variable):
                    x_init.append(var.starting_value)
                    self._state_variables.add_variable(var)
                else:
                    if var.fixed:
                        self._variables.append(var.value)
                    else:
                        self._variables.append(var.casadi_var)
                        self._unfixed_variables.add_variable(var)

        self._variables = ca.vcat(self._variables)

        if self.scaling is None:
            self.scaling = ca.DM.ones(self._variables.size())

        for time_step in self.time_grid[1:]:
            res_integration = self.integrator_tau(
                x0=x_init,
                p=ca.vertcat(
                    time_step - prev_time_step, self._variables * self.scaling
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]

            if time_step == self.time_grid[1]:
                res_states = res_integration["xf"]
            else:
                res_states = ca.horzcat(res_states, res_integration["xf"])

        return res_states

    def simulate(self, variable_list=None):
        # I actually do not need this method.
        if isinstance(self.evaluate_states, ca.DM):
            return self.evaluate_states
        else:
            values = []
            for var in self._unfixed_variables.keys():
                values.append(variable_list[var].value)

            function = ca.Function(
                "simulate",
                [self._unfixed_variables.get_casadi_var()],
                [self.evaluate_states],
            )
            return function(values)

    def generate_exp_data(self):
        res_array = self.create_simulation()
        variables = VariableList()
        for count, var in enumerate(self._state_variables.values()):
            new_var = copy.deepcopy(var)
            new_var.value = Experimental_Data()
            new_var.value.time = self.time_grid
            new_var.value.value = res_array[count, :].toarray()
            new_var.value.value = np.insert(new_var.value.value, 0, var.starting_value)
            variables.add_variable(new_var)
        return variables


class Experimental_Data(object):
    def __init__(self):
        self.time = None
        self.value = None

    def is_correct(self):
        if self.time.size == self.value.size:
            return True
        else:
            False


class ParameterEstimation(object):

    """Docstring for ParameterEstimation. """

    def __init__(self, model: Model, variable_list: VariableList):
        """TODO: to be defined. """
        self.model = model
        self.time_grid = np.ndarray((1, 0))
        self.decision_var = VariableList()
        self.states_experiment = VariableList()

        for var in variable_list.values():
            if isinstance(var, State_variable):
                if isinstance(var.value, Experimental_Data):
                    if var.value.is_correct():
                        self.time_grid = np.append(self.time_grid, var.value.time)
                        self.states_experiment.add_variable(var)
            elif isinstance(var, Parameter_variable):
                if var.fixed is False:
                    self.decision_var.add_variable(var)
            elif isinstance(var, Control_variable):
                var.fixed = True

        self.time_grid = np.unique(self.time_grid)

        self.simulation = Simulator(self.model, self.time_grid, variable_list)
        self.simulation.create_simulation()

    def calculate_objective(self):
        # Is done to rescale the simulation
        res = self.simulation.create_simulation()
        error = 0

        # count_state start from 0
        for count_state, var in enumerate(self.states_experiment.values()):
            # count_exp_point starts from 0, but we ignore first experimental column
            # so corresponding value is at count_exp_point + 1
            # res_index looks for index in array with 0 timepoint
            # and sim.evaluate_states doesn't include that time_point, thus - 1
            for count_exp_point, time_point in enumerate(var.value.time[1:]):
                res_index = np.nonzero(self.time_grid == time_point)
                res_index = res_index[0][0]
                error = (
                    error
                    + 0.5
                    * (
                        res[count_state, res_index - 1]
                        - var.value.value[count_exp_point + 1]
                    )
                    ** 2
                )

        return error

    def optimize(self, scale=True):
        # Scaling decreases amount of iterations, but ipopt fails gradient check at big amount of timestamps
        guess = []
        lb = []
        ub = []
        for var in self.decision_var.values():
            guess.append(var.guess)
            lb.append(var.lower_bound)
            ub.append(var.upper_bound)

        guess = np.array(guess)
        lb = np.array(lb)
        ub = np.array(ub)

        if scale:
            self.scaling = guess

            # for var in self.simulation._variables fails to iterate
            for count in range(self.simulation._variables.size()[0]):
                var = self.simulation._variables[count]
                if var.is_symbolic():
                    self.simulation.scaling[count] = self.decision_var[var.name()].guess
        else:
            self.scaling = 1
            self.simulation.scaling = None

        nlp = {"x": self.decision_var.get_casadi_var(), "f": self.calculate_objective()}

        nlp_solver = ca.nlpsol(
            "solver",
            "ipopt",
            nlp,
            {"verbose": False, "ipopt": {"derivative_test": "first-order"}},
        )

        res_solver = nlp_solver(
            x0=guess / self.scaling, lbx=lb / self.scaling, ubx=ub / self.scaling
        )
        print(res_solver["x"])
        print(res_solver["x"] * self.scaling)


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

    variable_list.add_variable(State_variable("e0_T", 273.0))
    variable_list.add_variable(State_variable("e0_c_i1", 3.0))
    variable_list.add_variable(State_variable("e0_c_i2", 10.0))
    variable_list.add_variable(State_variable("e0_c_i3", 0.0))
    variable_list.add_variable(State_variable("e0_c_i4", 0.0))

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

    time_grid = np.linspace(10, 1000, 2)
    time_grid = np.insert(time_grid, 0, 0)

    var_list1 = copy.deepcopy(variable_list)
    for var in var_list1.values():
        var.fixed = True
    s = Simulator(m, time_grid, var_list1)

    var_list_exp = s.generate_exp_data()

    var_list2 = copy.deepcopy(var_list1)
    for key, var in var_list_exp.items():
        var_list2[key] = var

    # TODO Fix the work if this is uncommented
    # var_list2["e0_T"].value = None

    var_list2["e0_U"].fixed = False
    var_list2["e0_U"].guess = 1.1
    var_list2["e0_U"].lower_bound = 1.0
    var_list2["e0_U"].upper_bound = 3.0

    var_list2["e0_k_pre_r1"].fixed = False
    var_list2["e0_k_pre_r1"].guess = 4000000.0
    var_list2["e0_k_pre_r1"].lower_bound = 4000000.0
    var_list2["e0_k_pre_r1"].upper_bound = 6000000.0

    var_list2["e0_k_pre_r2"].fixed = True
    var_list2["e0_k_pre_r2"].guess = 1.0e6
    var_list2["e0_k_pre_r2"].lower_bound = 1.0e6
    var_list2["e0_k_pre_r2"].upper_bound = 1.0e8

    var_list2["e0_k_pre_r3"].fixed = True
    var_list2["e0_k_pre_r3"].guess = 400000.0
    var_list2["e0_k_pre_r3"].lower_bound = 400000.0
    var_list2["e0_k_pre_r3"].upper_bound = 600000.0

    pe = ParameterEstimation(m, var_list2)
    pe.optimize()
    pe.optimize(False)


# # """
# # ODE Routine is currently not yeilding great results
# # """

# # # ODE Routine
# # res_states_ode = []
# # x_init = states_init
# # prev_time_step = 0
# # for time_step in time_grid[1:]:
# # res_integration_tau_ode = integrator_MAN_tau(
# # x0=x_init, p=ca.vertcat(parameters, controls, time_step - prev_time_step),
# # )

# # integrator_jac = integrator_MAN_tau.factory("I_fwd", ["x0", "p"], ["jac:xf:p"])

# # res_integration_jac = integrator_jac(
# # x0=x_init, p=ca.vertcat(parameters, controls, time_step - prev_time_step),
# # )

# # prev_time_step = time_step
# # x_init = res_integration_tau_ode["xf"]
# # if time_step == time_grid[1]:
# # res_jac_ode = res_integration_jac["jac_xf_p"][:, 0:4]
# # else:
# # res_jac_ode = ca.vertcat(res_jac_ode, res_integration_jac["jac_xf_p"][:, 0:4])

# # # res_jacobian = ca.jacobian(res_states_ode, ca.vertcat(parameters, controls))
# # # res_jacobian = ca.jacobian(res_states_ode, parameters)

# # # # Check objective
# # # eval_jacobian = ca.Function("eval_jacobian", [parameters, controls], [res_jac_ode])
# # # sensitivity_matrix = eval_jacobian(parameters_values, controls)[:, 0:4]

# # # # Calculate OBJ TRACE[FIM]
# # # sigma_diag = ca.DM([1, 1, 1, 1, 1]) * 1e20
# # # sigma_full = ca.diag(sigma_diag)

# # # measurement_matrix = ca.repmat(sigma_full, num_steps, num_steps)
# # # sensitivity_matrix = res_jac_ode
# # # fim_matrix = sensitivity_matrix.T @ measurement_matrix @ sensitivity_matrix

# # # eval_fim_matrix = ca.Function("eval_fim", [controls], [fim_matrix])
# # # fim_matrix_inv = ca.inv(fim_matrix)
# # # trace = ca.trace(fim_matrix_inv)

# # # eval_trace = ca.Function("eval_trace", [controls], [trace])

# # # fim = sensitivity_matrix.T @ measurement_matrix @ sensitivity_matrix
# # # trace = ca.trace(ca.inv(fim))
# # #
# # # trace = ca.trace(res_jac_ode@res_jac_ode.T)
# # # nlp_ode = {"x": ca.vertcat(parameters, controls), "f": trace}
# # # nlp_solver_ode = ca.nlpsol(
# # #     "solver",
# # #     "ipopt",
# # #     nlp_ode,
# # #     #     {"verbose": False, "ipopt": {"hessian_approximation": "exact", "max_iter": 200,"derivative_test": 'first-order'}},
# # #     {
# # #         "verbose": True,
# # #         "ipopt": {
# # #             "hessian_approximation": "limited-memory",
# # #             "max_iter": 200,
# # #             "derivative_test": "first-order",
# # #         },
# # #     },
# # # )
# # # #
# # # res_solver_ode = nlp_solver_ode(
# # #     x0=ca.vertcat(parameters_values, controls_guess),
# # #     lbx=ca.vertcat(parameters_values, controls_lb),
# # #     ubx=ca.vertcat(parameters_values, controls_ub),
# # # )
# # # # res_solver_ode = nlp_solver_ode(x0=controls_guess, lbx=controls_lb, ubx=controls_ub)
# # # print(res_solver_ode["x"])
