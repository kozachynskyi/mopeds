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
        self.parameters = VariableList()
        self.variables = VariableList()
        self.equations = None

        for var in variable_list.values():
            if isinstance(var, State_variable):
                self.states.add_variable(State_variable(var.name))
            elif isinstance(var, Parameter_variable):
                self.parameters.add_variable(Parameter_variable(var.name))
            else:
                raise (ValueError)

        self.variables.update(self.states)
        self.variables.update(self.parameters)

    def add_equations(self, equations):
        if self.equations is None:
            self.equations = ca.vcat(equations)
        else:
            raise (NotImplementedError)


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

# parameters_guess = ca.DM([4000000.0, 1.0e6, 400000.0, 2.4])
# parameters_lb = ca.DM([4000000.0, 1.0e6, 400000.0, 0.4])
# parameters_ub = ca.DM([6000000.0, 1.0e8, 600000.0, 3.4])
variable_list.add_variable(Parameter_variable("e0_k_pre_r1", 5000000.0))
variable_list.add_variable(Parameter_variable("e0_k_pre_r2", 1.0e7))
variable_list.add_variable(Parameter_variable("e0_k_pre_r3", 500000.0))
variable_list.add_variable(Parameter_variable("e0_U", 1.4))


variable_list.add_variable(Parameter_variable("e0_c_in_i1", 5.0))
variable_list.add_variable(Parameter_variable("e0_c_in_i2", 10.0))
variable_list.add_variable(Parameter_variable("e0_c_in_i3", 0.0))
variable_list.add_variable(Parameter_variable("e0_c_in_i4", 0.0))
variable_list.add_variable(Parameter_variable("e0_T_in", 373.0))
variable_list.add_variable(Parameter_variable("e0_T_j", 373.0))

m = Model(variable_list)

# fmt: off
tdot = (((((e0_F / e0_V) * ((m.variables["e0_T_in"].casadi_var - m.variables["e0_T"].casadi_var))) + (((m.variables["e0_U"].casadi_var * e0_A) / (e0_greek_rho * (e0_c_p * e0_V))) * ((m.variables["e0_T_j"].casadi_var - m.variables["e0_T"].casadi_var)))) + (((-e0_greek_Deltah_r1) / (e0_greek_rho * e0_c_p)) * (m.variables["e0_k_pre_r1"].casadi_var * (m.variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m.variables["e0_T"].casadi_var))))))) + (((-e0_greek_Deltah_r2) / (e0_greek_rho * e0_c_p)) * (m.variables["e0_k_pre_r2"].casadi_var * (m.variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m.variables["e0_T"].casadi_var))))))) + (((-e0_greek_Deltah_r3) / (e0_greek_rho * e0_c_p)) * (m.variables["e0_k_pre_r3"].casadi_var * (m.variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m.variables["e0_T"].casadi_var))))))
c1dot = ((((e0_F / e0_V) * ((m.variables["e0_c_in_i1"].casadi_var - m.variables["e0_c_i1"].casadi_var))) + (e0_greek_nu_i1_r1 * (m.variables["e0_k_pre_r1"].casadi_var * (m.variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m.variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r2 * (m.variables["e0_k_pre_r2"].casadi_var * (m.variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m.variables["e0_T"].casadi_var))))))) + (e0_greek_nu_i1_r3 * (m.variables["e0_k_pre_r3"].casadi_var * (m.variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m.variables["e0_T"].casadi_var))))))
c2dot = ((e0_F / e0_V) * ((m.variables["e0_c_in_i2"].casadi_var - m.variables["e0_c_i2"].casadi_var))) + (e0_greek_nu_i2_r2 * (m.variables["e0_k_pre_r2"].casadi_var * (m.variables["e0_c_i2"].casadi_var * ca.exp(((-e0_E_r2) / (e0_R * m.variables["e0_T"].casadi_var))))))
c3dot = ((e0_F / e0_V) * ((m.variables["e0_c_in_i3"].casadi_var - m.variables["e0_c_i3"].casadi_var))) + (e0_greek_nu_i3_r1 * (m.variables["e0_k_pre_r1"].casadi_var * (m.variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r1) / (e0_R * m.variables["e0_T"].casadi_var))))))
c4dot = ((e0_F / e0_V) * ((m.variables["e0_c_in_i4"].casadi_var - m.variables["e0_c_i4"].casadi_var))) + (e0_greek_nu_i4_r3 * (m.variables["e0_k_pre_r3"].casadi_var * (m.variables["e0_c_i1"].casadi_var * ca.exp(((-e0_E_r3) / (e0_R * m.variables["e0_T"].casadi_var))))))
# fmt: on

m.add_equations([tdot, c1dot, c2dot, c3dot, c4dot])


class Simulator(object):

    """Docstring for Simulator. """

    def __init__(self, model: Model, time_grid):
        """TODO: to be defined. """
        self.model = model
        self.tau = ca.MX.sym("tau")
        self.time_grid = time_grid

        self.ode_system = {
            "x": self.model.states.get_casadi_var(),
            "p": ca.vertcat(self.model.parameters.get_casadi_var()),
            "ode": self.model.equations,
        }
        self.ode_system_tau = {
            "x": self.model.states.get_casadi_var(),
            "p": ca.vertcat(self.tau, self.model.parameters.get_casadi_var()),
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

    def create_simulation(self, variable_list: VariableList):
        self._parameters = []
        self._unfixed_parameters = VariableList()
        prev_time_step = 0
        res_states = []
        x_init = []

        for var in variable_list.values():
            if isinstance(var, State_variable):
                x_init.append(var.starting_value)
            elif isinstance(var, Parameter_variable):
                if var.fixed:
                    self._parameters.append(var.value)
                else:
                    self._parameters.append(var.casadi_var)
                    self._unfixed_parameters.add_variable(var)

        self._parameters = ca.vcat(self._parameters)

        for time_step in self.time_grid[1:]:
            res_integration = self.integrator_tau(
                x0=x_init, p=ca.vertcat(time_step - prev_time_step, self._parameters)
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]

            if time_step == self.time_grid[1]:
                res_states = res_integration["xf"]
            else:
                res_states = ca.horzcat(res_states, res_integration["xf"])

        self.evaluate_states = res_states

    def simulate(self, variable_list=None):
        # I actually do not need this method.
        if isinstance(self.evaluate_states, ca.DM):
            return self.evaluate_states
        else:
            values = []
            for var in self._unfixed_parameters.keys():
                values.append(variable_list[var].value)

            function = ca.Function(
                "simulate",
                [self._unfixed_parameters.get_casadi_var()],
                [self.evaluate_states],
            )
            return function(values)

    def generate_exp_data(self):
        res_array = self.evaluate_states
        variables = VariableList()
        for count, var in enumerate(self.model.states.values()):
            new_var = copy.deepcopy(var)
            new_var.value = Experimental_Data()
            new_var.value.time = self.time_grid
            new_var.value.value = res_array[count, :].toarray()
            new_var.value.value = np.insert(new_var.value.value,0,var.starting_value)
            variables.add_variable(new_var)
        return variables


time_grid = [0]
time_grid.extend(np.linspace(time_grid[-1] + 10, 1000, 100).tolist())
s = Simulator(m, time_grid)
s2 = Simulator(m, time_grid)

var_list1 = copy.deepcopy(variable_list)
for var in var_list1.values():
    var.fixed = True
var_list1["e0_U"].fixed = False
s.create_simulation(var_list1)
s.simulate(var_list1)

var_list1["e0_U"].fixed = True
s2.create_simulation(var_list1)
s2.simulate()


class Experimental_Data(object):
    def __init__(self):
        self.time = []
        self.value = []


exit()


class ParameterEstimation(object):

    """Docstring for ParameterEstimation. """

    def __init__(self, model: Model, exp_data):
        """TODO: to be defined. """
        self.model = model
        self.exp_data = exp_data

    def set_exp_data(self, exp_data):
        self._exp_data = exp_data

    def calculate_objective(self):
        x_init = self.simulator.model.states_init_values
        prev_time_step = 0
        res_states = []

        # This loop runs iterator for given time_points
        for time_step in self.simulator.time_grid[1:]:
            print(time_step)
            res_integration = self.simulator.integrator_tau(
                x0=x_init,
                p=ca.vertcat(
                    time_step - prev_time_step,
                    self.simulator.model.parameters * self.parameters_scale,
                    self.controls_values,
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]
            if time_step == self.simulator.time_grid[1]:
                res_states = res_integration["xf"]
            else:
                res_states = ca.horzcat(res_states, res_integration["xf"])

        self.res_states = res_states
        self.error = self._exp_data - self.res_states
        self.objective = 0.5 * ca.dot(self.error, self.error)
        return self.objective

    def eval_res_states(self, parameters_values):
        function = ca.Function(
            "eval_states", [self.simulator.model.parameters], [self.error]
        )
        res = function(parameters_values)
        return res

    def eval_error(self, parameters_values):
        function = ca.Function(
            "eval_error", [self.simulator.model.parameters], [self.error]
        )
        res = function(parameters_values)
        return res

    def eval_objective(self, parameters_values):
        function = ca.Function(
            "eval_objective", [self.simulator.model.parameters], [self.objective]
        )
        res = function(parameters_values)
        return res

    def optimize(self, parameters_guess, lower_bound, upper_bound):
        # Settings for optimizer are set here
        self.parameters_scale = parameters_guess
        nlp = {"x": self.simulator.model.parameters, "f": self.calculate_objective()}
        nlp_solver = ca.nlpsol(
            "solver",
            "ipopt",
            nlp,
            {"verbose": False, "ipopt": {"derivative_test": "first-order"}},
        )

        res_solver = nlp_solver(
            x0=parameters_guess / self.parameters_scale,
            lbx=lower_bound / self.parameters_scale,
            ubx=upper_bound / self.parameters_scale,
        )
        print(res_solver["x"])
        print(res_solver["x"] * self.parameters_scale)


if __name__ == "__main__":
    time_grid = [0]
    time_grid.extend(np.linspace(time_grid[-1] + 10, 1000, 100).tolist())

    model = Model([1, None, 1, None], 1)
    sim = Simulator(model, time_grid)

    # Call integrator to generate experimental data for PE
    res_integration = sim.integrator(
        x0=sim.model.states_init_values,
        p=ca.vertcat(sim.model.parameters_values, sim.model.controls_values),
    )
    exp_data = res_integration["xf"]
    print(exp_data)

    pe = ParameterEstimation(sim, model.controls_values)
    pe.set_exp_data(exp_data)
    # obj = pe.calculate_objective()
    # print(obj)

    # print(pe.eval_error(model.parameters_values))
    # print(pe.eval_objective(model.parameters_values))

    parameters_guess = ca.DM([4000000.0, 1.0e6, 400000.0, 2.4])
    parameters_lb = ca.DM([4000000.0, 1.0e6, 400000.0, 0.4])
    parameters_ub = ca.DM([6000000.0, 1.0e8, 600000.0, 3.4])

    pe.optimize(parameters_guess, parameters_lb, parameters_ub)

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
