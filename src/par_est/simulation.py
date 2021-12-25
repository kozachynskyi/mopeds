import copy
import logging
from typing import List, Union

import casadi as ca
import numpy as np

from par_est import (
    BadVariableError,
    Model,
    VariableAlgebraic,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableList,
    VariableParameter,
    VariableState,
)


class Simulator(object):
    supported_integrators = ["idas", "cvodes", "collocation"]

    def __init__(  # noqa: C901
        self,
        model: Model,
        input_time_grid: np.ndarray,
        variable_list: VariableList,
        integrator_name="idas",
        integrator_settings=None,
        *,
        use_idas_constraints=False,
        simulate_jac=False,
    ):

        self.logger = logging.getLogger(__name__)
        self.logger.debug(
            "Creating Simulator object: \n timegrid \n {0} \n".format(input_time_grid)
        )

        self.__input_variable_list: VariableList = copy.deepcopy(variable_list)
        self.model = model

        if integrator_name not in self.supported_integrators:
            raise TypeError(f"Provided integrator name {integrator_name} is not supported. Only theese are: {self.supported_integrators}.")
        self.__integrator_name = integrator_name

        self.scaling = None

        self.setup_time_grid(input_time_grid)

        self.ode_system = {
            "x": self.model.varlist_state.get_casadi_variables(),
            "p": ca.vertcat(self.model.varlist_independent.get_casadi_variables()),
            "ode": self.model.equations_differential,
        }

        # Tau variable is used to specify a length of iteration step externally, via tau variable
        self.tau = ca.MX.sym("tau")
        self.ode_system_tau = {
            "x": self.model.varlist_state.get_casadi_variables(),
            "p": ca.vertcat(self.tau, self.model.varlist_independent.get_casadi_variables()),
            "ode": self.model.equations_differential * self.tau,
        }

        if self.model.DAE:
            self.ode_system["alg"] = self.model.equations_algebraic
            self.ode_system_tau["alg"] = self.model.equations_algebraic
            self.ode_system["z"] = self.model.varlist_algebraic.get_casadi_variables()
            self.ode_system_tau["z"] = self.model.varlist_algebraic.get_casadi_variables()

        if integrator_settings is not None:
            self.__integrator_settings = integrator_settings
        else:
            self._set_default_integrator_settings()

        self._setup_constraints_idas(use_idas_constraints)

        # TODO This integrator is not used so far...
        self.integrator = ca.integrator(
            "integrator",
            self.__integrator_name,
            self.ode_system,
            {"grid": self.time_grid_relative, "output_t0": False, "print_stats": True},
        )

        self.integrator_tau = ca.integrator(
            "integrator_tau",
            self.__integrator_name,
            self.ode_system_tau,
            self.__integrator_settings,
        )

        # This integrator is used to output values of algebraic variables at time 0
        # and should be run first to get algebraic variables at time 0 for whole simulation
        integrator_settings_with_output_t0 = copy.deepcopy(self.__integrator_settings)
        integrator_settings_with_output_t0["output_t0"] = True
        self.integrator_tau_with_t0 = ca.integrator(
            "integrator_tau_with_t0",
            self.__integrator_name,
            self.ode_system_tau,
            integrator_settings_with_output_t0,
        )

        """ This nested list holds either a value or a casadi variable of
        each independent variable at every timestamp in time_grid.  First it has form:
        [[var1_t0, var1_t1 ...], [var2_t0, var2_t1 ...], [...]]
        Than it's reformed to:
        [[var1_t0, var2_t0 ...], [var1_t1, var2_t1 ...], [...]]
        And finally nested lists are changed to casadi.MX or DM vectors
        [ca.MX(var1_t0, var2_t0 ...), ca.MX(var1_t1, var2_t1 ...), [...]]
        """
        self._independent_variables: List[Union[float, ca.MX]] = []
        # List of values of State Variables at time 0
        self._initial_state: List[float] = []

        # List of expected values of Algebraic Variables at time 0
        # This list is stored in order to retain information about original guess
        self._initial_algebraic_original: List[float] = []

        # List of expected or recalculated values of Algebraic Variables at time 0
        # This list is further used in calculations
        self._initial_algebraic: List[float] = []

        # This list is used for utility functions, like finding steady state
        self._guess_or_value_of_independent_variables: List[float] = []

        # Here all lists from above are initialized
        self._setup_variables()

        # .factory() method is very expensive so should be requested externally
        if simulate_jac:
            if self.model.DAE is True:
                self.integrator_tau_jac = self.integrator_tau.factory(
                    "integrator_tau_jacobian",
                    self.integrator_tau.name_in(),
                    ["xf", "qf", "zf", "rxf", "rqf", "rzf", "jac:xf:p"],
                )
            else:
                self.integrator_tau_jac = self.integrator_tau.factory(
                    "integrator_tau_jacobian",
                    self.integrator_tau.name_in(),
                    ["xf", "qf", "rxf", "rqf", "jac:xf:p"],
                )

        self._reset_scaling()

        # This code is moved here, so this if statement shouldn't be called every simulation
        if self.model.DAE is True:
            self.simulate = self._simulate_dae
            self.simulate_jac = self._simulate_jac_dae
        else:
            self.simulate = self._simulate_ode
            self.simulate_jac = self._simulate_jac_ode

    def _reset_scaling(self):
        self.scaling = ca.DM.ones(self._independent_variables[0].size())

    def _setup_constraints_idas(self, use_idas_constraints):
        """Holds a list of constraints for state and algebraic variables
        which can be used to constrain the solution of the idas
        to positive or negative numbers """
        self._constraints_idas = []
        variable_names = list(self.model.varlist_state.keys())
        if self.model.DAE:
            variable_names.extend(list(self.model.varlist_algebraic.keys()))
        for var_name in variable_names:
            self._constraints_idas.append(
                self.__input_variable_list[var_name].get_constraint_idas
            )

        if use_idas_constraints:
            if not self.__integrator_name == "idas":
                self.logger.warning(
                    "use_idas_constraints argument is applicable only for idas solver"
                )
            else:
                if all(constraint == 0 for constraint in self._constraints_idas):
                    self.logger.warning(
                        "All idas constraints are 0, so no option is set"
                    )
                else:
                    self.__integrator_settings["constraints"] = self._constraints_idas

    def _setup_variables(self):
        """ Setup all important lists for simulator"""
        num_time_steps = len(self.time_grid_relative) - 1

        for variable_name in self.model.varlist_all.keys():
            try:
                var = self.__input_variable_list[variable_name]
            except KeyError:
                continue

            if isinstance(var, VariableState):
                try:
                    self._initial_state.append(var.value[0])
                except Exception as e:
                    raise (BadVariableError(var)) from e

            elif isinstance(var, VariableAlgebraic):
                self._initial_algebraic.append(var.guess)

            elif isinstance(var, VariableParameter):
                independent_variable = []
                independent_variable.extend(
                    [var.get_value_or_casadi()] * num_time_steps
                )
                self._guess_or_value_of_independent_variables.append(var.get_value_or_guess())
                self._independent_variables.append(independent_variable)
            elif isinstance(var, VariableControl):
                if isinstance(var, VariableControlPiecewiseConstant):
                    var_t0 = var.get_variable_at_time_relative(0)
                    self._guess_or_value_of_independent_variables.append(var_t0.get_value_or_guess())
                    independent_variable = var.get_value_or_casadi(self.time_grid_relative)
                else:
                    independent_variable = []
                    self._guess_or_value_of_independent_variables.append(var.get_value_or_guess())
                    independent_variable.extend([var.get_value_or_casadi()] * num_time_steps)

                self._independent_variables.append(independent_variable)

        # Groups nested lists by time_stamp
        self._independent_variables = list(map(list, zip(*self._independent_variables)))
        self._initial_algebraic_original = copy.deepcopy(self._initial_algebraic)

        # Transforms nested lists in ca.MX or ca.DM array
        for index, column in enumerate(self._independent_variables):
            casadi_mx = ca.vcat(column)
            self._independent_variables[index] = casadi_mx

    def _set_default_integrator_settings(self):
        """ Sane default settings for integrators"""
        if self.__integrator_name == "idas":
            self.__integrator_settings = {
                "tf": 1,
                "expand": True,
                # "calc_ic": False,
                # 'abstol': 1,
                # "reltol": 1,
                # "monitor": "jacF",
                # "print_in": True,
                # "print_out": True,
                # "verbose": True,
                # "print_stats": True,
            }
        elif self.__integrator_name == "cvodes":
            self.__integrator_settings = {
                "tf": 1,
                "expand": True,
                # "linear_multistep_method": "adams",# was used for CVODES
                # "output_t0": False,
                # "use_preconditioner": False,
                # "calc_ic": False,
                # 'abstol': 1e-5,
                # "reltol": 1e-5,
                # "monitor": "jacF",
                # "print_in": True,
                # "print_out": True,
                # "verbose": True,
                # "print_stats": True,
            }
        elif self.__integrator_name == "collocation":
            self.__integrator_settings = {
                "number_of_finite_elements": 3,
                "simplify": True,
                "expand": True,
                "rootfinder": "fast_newton",
                # "monitor": "jacF",
                # "print_in": True,
                # "print_out": True,
                # "verbose": True,
                # "print_stats": True,
            }

    def calculate_steady_state(self):
        if self.model.DAE:
            steady_state_rootfinder = ca.Function(
                "steadystate_eq_sys",
                [
                    ca.vertcat(self.ode_system["x"], self.ode_system["z"]),
                    self.ode_system["p"],
                ],
                [ca.vertcat(self.ode_system["ode"], self.ode_system["alg"])],
                ["x", "p"],
                ["ode_alg"],
            )
        else:
            steady_state_rootfinder = ca.Function(
                "steadystate_eq_sys",
                [
                    self.ode_system["x"],
                    self.ode_system["p"],
                ],
                [self.ode_system["ode"]],
                ["x", "p"],
                ["ode"],
            )

        # rf_settings = {
        #     # "calc_ic": False,
        #     # 'abstol': 1,
        #     # "reltol": 1,
        #     # "monitor": "jacF",
        #     # "print_in": True,
        #     # "nlpsol": "ipopt",
        #     # "print_out": True,
        #     # "verbose": True,
        #     # "print_stats": True,
        #     }
        # rf_settings["nlpsol"] = "ipopt"
        # rf_steadystate = ca.rootfinder("stea_state", "nlpsol", steady_state_rootfinder, rf_settings)
        rf_steadystate = ca.rootfinder("stea_state", "newton", steady_state_rootfinder)

        if self.model.DAE:
            res_steadystate = rf_steadystate(
                ca.vertcat(self._initial_state, self._initial_algebraic_original),
                self._guess_or_value_of_independent_variables,
            )
        else:
            res_steadystate = rf_steadystate(
                self._initial_state,
                self._guess_or_value_of_independent_variables,
            )
        return res_steadystate

    def calculate_algebraic_initials(self, *, apply_intials=False, analyze=False, experimental=False):
        function = ca.Function(
            "eq_sys",
            [self.ode_system["x"], self.ode_system["z"], self.ode_system["p"]],
            [self.ode_system["ode"], self.ode_system["alg"]],
            ["x", "z", "p"],
            ["ode", "alg"],
        )

        algebraic_eqsys_rootfinder = ca.Function(
            "alg_eq_sys",
            [
                self.ode_system["z"],
                ca.vertcat(self.ode_system["x"], self.ode_system["p"]),
            ],
            [self.ode_system["alg"]],
            ["x", "p"],
            ["alg"],
        )

        rf = ca.rootfinder("inits", "newton", algebraic_eqsys_rootfinder)

        if experimental:
            settings = {
                "tf": 1,
                "expand": True,
                "output_t0": True,
            }
            integrator = ca.integrator(
                "integrator_tau",
                "idas",
                self.ode_system_tau,
                settings,
            )

            x_init = self._initial_state
            alg_init = self._initial_algebraic

            res_integration = integrator(
                x0=x_init,
                z0=alg_init,
                p=ca.vertcat(
                    1, self._guess_or_value_of_independent_variables * self.scaling
                ),
            )

            res = res_integration["zf"][:,0]

        else:
            res = rf(
                self._initial_algebraic_original,
                ca.vertcat(self._initial_state, self._guess_or_value_of_independent_variables),
            )

        residual_original = function(
            x=self._initial_state,
            z=self._initial_algebraic_original,
            p=self._guess_or_value_of_independent_variables,
        )
        residual_calculated = function(
            x=self._initial_state, z=res, p=self._guess_or_value_of_independent_variables
        )

        if analyze:
            abs_diff = self._initial_algebraic_original - res
            rel_diff = ca.fabs(abs_diff) / ca.fabs(self._initial_algebraic_original)

            print("Prints Algebraic Variables, that we changed more than 50%")
            for i in range(abs_diff.shape[0]):
                if rel_diff[i] > 0.50:
                    print(self.ode_system["z"][i])
                    print(f"Value After {res[i]}")
                    print(f"Value Before {self._initial_algebraic_original[i]}")

            residual_sum_original = ca.sum1(residual_original["alg"])
            residual_sum_calculated = ca.sum1(residual_calculated["alg"])
            print(
                f"Residual before {residual_sum_original}, after {residual_sum_calculated}."
            )
        if apply_intials:
            self.logger.debug("Fixed algebraic intials")
            self._initial_algebraic = res

    def analyze_WIP(self, state_value=None):
        import par_est.tools as tools  # noqa: F401

        """ This function was working for previous version of the module."""
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

        check_initials = function(  # noqa: F841
            x=self._initial_state, z=self._initial_algebraic, p=self._independent_variables[0]
        )
        jacobian = function.factory(
            "jac_alg",
            function.name_in(),
            ["jac:alg:z", "jac:alg:x", "jac:ode:x", "jac:ode:z"],
        )
        check_jacobian = jacobian(
            x=self._initial_state, z=self._initial_algebraic, p=self._independent_variables[0]
        )

        check_alg = algebraic_eqsys(
            x=self._initial_state, z=self._initial_algebraic, p=self._independent_variables[0]
        )

        # should fail by DAE index > 1
        ca.inv(check_jacobian["jac_alg_z"])
        # tools.plot_array(check_jacobian["jac_alg_z"], self.model.varlist_algebraic)
        # tools.plot_array(check_jacobian["jac_ode_z"], self.model.varlist_algebraic)
        # tools.plot_array(check_jacobian["jac_alg_x"], self.model.varlist_state)
        # tools.plot_array(check_jacobian["jac_ode_x"], self.model.varlist_state)

        algebraic_eqsys_rootfinder = ca.Function(
            "alg_eq_sys",
            [
                self.ode_system["z"],
                ca.vertcat(self.ode_system["x"], self.ode_system["p"]),
            ],
            [self.ode_system["alg"]],
            ["x", "p"],
            ["alg"],
        )

        rf = ca.rootfinder("inits", "newton", algebraic_eqsys_rootfinder)
        if state_value is not None:
            res = rf(
                self._initial_algebraic, ca.vertcat(state_value, self._independent_variables[0])
            )
        else:
            res = rf(
                self._initial_algebraic,
                ca.vertcat(self._initial_state, self._independent_variables[0]),
            )

            check_alg = function(  # noqa: 841
                x=self._initial_state, z=res, p=self._independent_variables[0]
            )
        old_initial = self._initial_algebraic
        # self._initial_algebraic = res
        return [res, old_initial]
        # return check_initials, check_jacobian

    def _simulate_jac_dae(self):
        """Return dictionary with results "xf" - state,
        "zf" - algebraic, "jac_xf_p" - derivatives.
        """
        prev_time_step = 0
        res_states = []
        res_algebraic = []
        res_jacobian = []
        x_init = self._initial_state
        alg_init = self._initial_algebraic

        for time_step, independent_variables in zip(
            self.time_grid_relative[1:], self._independent_variables
        ):
            res_integration = self.integrator_tau_jac(
                x0=x_init,
                z0=alg_init,
                p=ca.vertcat(
                    time_step - prev_time_step, independent_variables * self.scaling
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]
            alg_init = res_integration["zf"]

            res_states.append(res_integration["xf"])
            res_algebraic.append(res_integration["zf"])
            res_jacobian.append(res_integration["jac_xf_p"])

        res_states = ca.hcat(res_states)
        res_algebraic = ca.hcat(res_algebraic)
        res_jacobian = ca.hcat(res_jacobian)

        res = {"xf": res_states, "zf": res_algebraic, "jac_xf_p": res_jacobian}
        return res

    def _simulate_jac_ode(self):
        """Return dictionary with results "xf" - state,
        "zf" - algebraic, "jac_xf_p" - derivatives.
        """
        prev_time_step = 0
        res_states = []
        res_jacobian = []
        x_init = self._initial_state

        for time_step, independent_variables in zip(
            self.time_grid_relative[1:], self._independent_variables
        ):
            res_integration = self.integrator_tau_jac(
                x0=x_init,
                p=ca.vertcat(
                    time_step - prev_time_step, independent_variables * self.scaling
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]

            res_states.append(res_integration["xf"])
            res_jacobian.append(res_integration["jac_xf_p"])

        res_states = ca.hcat(res_states)
        res_jacobian = ca.hcat(res_jacobian)

        res = {"xf": res_states, "jac_xf_p": res_jacobian}
        return res

    def _simulate_dae(self):
        """Return dictionary with results "xf" - state,
        "zf" - algebraic
        """
        prev_time_step = 0
        res_states = []
        res_algebraic = []
        x_init = self._initial_state
        alg_init = self._initial_algebraic

        for time_step, independent_variables in zip(
            self.time_grid_relative[1:], self._independent_variables
        ):
            res_integration = self.integrator_tau(
                x0=x_init,
                z0=alg_init,
                p=ca.vertcat(
                    time_step - prev_time_step, independent_variables * self.scaling
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]
            alg_init = res_integration["zf"]

            res_states.append(res_integration["xf"])
            res_algebraic.append(res_integration["zf"])

        res_states = ca.hcat(res_states)
        res_algebraic = ca.hcat(res_algebraic)

        res = {"xf": res_states, "zf": res_algebraic}
        return res

    def _simulate_ode(self):
        """Return dictionary with results "xf" - state,
        "zf" - algebraic
        """
        prev_time_step = 0
        res_states = []
        x_init = self._initial_state

        for time_step, independent_variables in zip(
            self.time_grid_relative[1:], self._independent_variables
        ):
            res_integration = self.integrator_tau(
                x0=x_init,
                p=ca.vertcat(
                    time_step - prev_time_step, independent_variables * self.scaling
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]

            res_states.append(res_integration["xf"])

        res_states = ca.hcat(res_states)

        res = {"xf": res_states}
        return res

    def generate_exp_data(self, algebraic=False, recalculate_algebraic=True):
        """ Runs simulation and returns results in VariableList class."""
        variables = VariableList()

        if recalculate_algebraic and self.model.DAE:
            self.calculate_algebraic_initials(apply_intials=True)

        result_simulation = self.simulate()
        if not algebraic or not self.model.DAE:
            result_varlist = [copy.deepcopy(self.model.varlist_state)]
            res_array = result_simulation["xf"]
        else:
            result_varlist = [
                copy.deepcopy(self.model.varlist_state),
                copy.deepcopy(self.model.varlist_algebraic),
            ]
            res_array = ca.vertcat(result_simulation["xf"], result_simulation["zf"])

        if not isinstance(res_array, ca.DM):
            raise NotImplementedError("Generation of experimental data is possible only for simulations with fixed independent variables")

        shift_by = 0
        for variable_list in result_varlist:
            for count, var in enumerate(variable_list.values()):
                var.casadi_var = None
                new_var = copy.deepcopy(var)

                if isinstance(var, VariableAlgebraic):
                    value_time_zero = var.guess
                elif isinstance(var, VariableState):
                    value_time_zero = (
                        self.__input_variable_list[var.name].value[0],
                    )
                else:
                    raise (NotImplementedError)

                value = res_array[count + shift_by, :]
                value = np.insert(value, 0, value_time_zero)

                new_var.set_dataframe_from_value_and_time(value, self.time_grid_relative, self.origin_ts)
                new_var.ignore_plotting = self.__input_variable_list[var.name].ignore_plotting

                variables.add_variable(new_var)
            shift_by = count + 1

        return variables

    def setup_time_grid(self, time_grid):
        """Time_grid provided by user may not take into account piecewise controls.
        Thus it might be needed to expand a time grid."""
        for var in self.__input_variable_list.values():
            if isinstance(var, VariableControlPiecewiseConstant):
                time_grid = np.append(time_grid, var.time_relative)

        # Values of provided time_grid are rounded to milisecconds
        # in order to avoid timestamps that are very close to each other
        self.time_grid_relative: np.ndarray = np.unique(time_grid.round(decimals=1))
        self.origin_ts = self.__input_variable_list.get_common_origin()
        self.logger.debug(
            "Timegrid modified: \n self.timegrid \n {0} \n".format(self.time_grid_relative)
        )
