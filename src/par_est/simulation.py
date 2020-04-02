import copy

import casadi as ca
import numpy as np

from par_est import (
    VariableList,
    Model,
    ExperimentData,
    Variable,
    VariableState,
    VariableAlgebraic,
)


class Simulator(object):

    """Docstring for Simulator. """

    def __init__(self, model: Model, time_grid, variable_list: VariableList):
        """TODO: to be defined. """
        self.__input_variable_list = copy.deepcopy(variable_list)
        self.model = model
        self.tau = ca.MX.sym("tau")
        self.scaling = None
        self.time_grid = time_grid

        self.ode_system = {
            "x": self.model.varlist_state.get_casadi_var(),
            "p": ca.vertcat(self.model.variables.get_casadi_var()),
            "ode": self.model.differential_equations,
        }
        self.ode_system_tau = {
            "x": self.model.varlist_state.get_casadi_var(),
            "p": ca.vertcat(self.tau, self.model.variables.get_casadi_var()),
            "ode": self.model.differential_equations * self.tau,
        }

        if self.model.algebraic_equations is not None:
            self.ode_system["alg"] = self.model.algebraic_equations
            self.ode_system_tau["alg"] = self.model.algebraic_equations
            self.ode_system["z"] = self.model.varlist_algebraic.get_casadi_var()
            self.ode_system_tau["z"] = self.model.varlist_algebraic.get_casadi_var()

        self.integrator = ca.integrator(
            "integrator",
            "idas",
            self.ode_system,
            {"grid": self.time_grid, "output_t0": False, "print_stats": False},
        )

        # This integrator uses tau variable and is used in PE and OED
        self.integrator_tau = ca.integrator(
            "integrator",
            "idas",
            self.ode_system_tau,
            {
                "tf": 1,
                "output_t0": False,
                "print_stats": False,
                "calc_ic": True,
                # "linear_multistep_method": "adams", # was used for CVODES
            },
        )

        self._state_variables = VariableList()
        self._algebraic_variables = VariableList()
        # Arrays needed to initialize integrator.
        self._variables = []
        self._initial_states = []
        self._initial_alg = []

        for var in self.__input_variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    self._initial_states.append(var.starting_value)
                    self._state_variables.add_variable(var)
                elif isinstance(var, VariableAlgebraic):
                    self._initial_alg.append(var.starting_value)
                    self._algebraic_variables.add_variable(var)
                else:
                    if var.fixed:
                        self._variables.append(var.value)
                    else:
                        self._variables.append(var.casadi_var)

        self._variables = ca.vcat(self._variables)
        self._reset_scaling()

    def _reset_scaling(self):
        self.scaling = ca.DM.ones(self._variables.size())

    def analyze(self):
        # ca.Function("equations", [self._state_variables, self._algebraic_variables, self._variables], [self.model.algebraic_equations, self.model.algebraic_equations], [
        function = ca.Function(
            "eq_sys",
            [self.ode_system["x"], self.ode_system["z"], self.ode_system["p"]],
            [self.ode_system["ode"], self.ode_system["alg"]],
            ["x", "z", "p"],
            ["ode", "alg"],
        )

        algebraic_eqsys = ca.Function(
            "alg_eq_sys",
            [self.ode_system["x"], self.ode_system["z"], self.ode_system["p"]],
            [self.ode_system["alg"]],
            ["x", "z", "p"],
            ["alg"],
        )

        check_initials = function(x=self._initial_states, z=self._initial_alg, p=self._variables)
        jacobian = function.factory("jac_alg", function.name_in(), ["jac:alg:z"])
        check_jacobian = jacobian(x=self._initial_states, z=self._initial_alg, p=self._variables)

        check_alg = algebraic_eqsys(x=self._initial_states, z=self._initial_alg, p=self._variables)

        # should fail by DAE index > 1
        ca.inv(check_jacobian["jac_alg_z"])

        algebraic_eqsys_rootfinder = ca.Function(
            "alg_eq_sys",
            [self.ode_system["z"], ca.vertcat(self.ode_system["x"], self.ode_system["p"])],
            [self.ode_system["alg"]],
            ["x", "p"],
            ["alg"],
        )

        rf = ca.rootfinder("inits", "newton", algebraic_eqsys_rootfinder)
        res = rf(self._initial_alg, ca.vertcat(self._initial_states, self._variables))
        
        check_alg = function(x=self._initial_states, z=res, p=self._variables)
        print(check_alg)
        self._initial_alg = res
        return [res, self._initial_alg]
        # return check_initials, check_jacobian

    def simulate(self, derivatives=False):
        prev_time_step = 0
        res_states = []
        res_algebraic = []
        res_jacobian = []
        x_init = self._initial_states
        alg_init = self._initial_alg

        for time_step in self.time_grid[1:]:
            res_integration = self.integrator_tau(
                x0=x_init,
                z0=alg_init,
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
            alg_init = res_integration["zf"]

            if time_step == self.time_grid[1]:
                res_states = res_integration["xf"]
                res_algebraic = res_integration["zf"]
                if derivatives:
                    res_jacobian = res_integration_jac["jac_xf_p"]
            else:
                res_states = ca.horzcat(res_states, res_integration["xf"])
                res_algebraic = ca.horzcat(res_algebraic, res_integration["zf"])
                if derivatives:
                    res_jacobian = ca.vertcat(
                        res_jacobian, res_integration_jac["jac_xf_p"]
                    )

        if derivatives:
            return res_states, res_jacobian
        else:
            return ca.vertcat(res_states, res_algebraic)

    def generate_exp_data(self):
        res_array = self.simulate()
        variables = VariableList()
        convert_to_numpy = False
        if isinstance(res_array, ca.DM):
            convert_to_numpy = True

        shift_by = 0
        for variable_list in [self._state_variables, self._algebraic_variables]:
            for count, var in enumerate(variable_list.values()):
                new_var = copy.deepcopy(var)
                new_var.value = ExperimentData()
                if convert_to_numpy:
                    new_var.value.time = self.time_grid
                    new_var.value.value = res_array[count + shift_by, :].toarray()
                    new_var.value.value = np.insert(
                        new_var.value.value, 0, var.starting_value
                    )
                else:
                    new_var.value.time = self.time_grid[1:]
                    new_var.value.value = res_array[count + shift_by, :]

                variables.add_variable(new_var)
            shift_by = count + 1

        return variables
