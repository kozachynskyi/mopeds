import copy

import casadi as ca
import numpy as np

from par_est import VariableList, Model, Experimental_Data, Variable, State_variable


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

        self._state_variables = VariableList()
        # Arrays needed to initialize integrator.
        self._variables = []
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
