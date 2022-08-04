from __future__ import annotations

import copy
import logging
from typing import cast, Callable

import casadi as ca
import numpy as np

from par_est import (
    BadVariableError,
    Model,
    VariableAlgebraic,
    VariableConstant,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableList,
    VariableParameter,
    VariableState,
)


class Simulator(object):
    """recalculate_algebraic - set to "True" in order to recalculate initial algebraic values "z0"
    Use this option, if optimizer runs into IDAS_CALC_IC problems and is not able to find "z0" by itself
    """

    supported_integrators = ["idas", "cvodes", "collocation"]

    def __init__(  # noqa: C901
        self,
        model: Model,
        input_time_grid: np.ndarray,
        variable_list: VariableList,
        integrator_name: str = "idas",
        integrator_settings: dict | None = None,
        *,
        use_idas_constraints: bool = False,
        simulate_jac: bool = False,
        recalculate_algebraic: bool = False,
    ) -> None:

        self.logger: logging.Logger = logging.getLogger(__name__)
        self.logger.debug(
            "Creating Simulator object: \n timegrid \n {0} \n".format(input_time_grid)
        )

        self.__input_variable_list: VariableList = copy.deepcopy(variable_list)
        self.model: Model = model

        if integrator_name not in self.supported_integrators:
            raise TypeError(
                f"Provided integrator name {integrator_name} is not supported. Only theese are: {self.supported_integrators}."
            )
        self.__integrator_name: str = integrator_name

        self.setup_time_grid(input_time_grid)

        self.ode_system: dict[str, ca.MX] = {
            "x": self.model.varlist_state.get_casadi_variables(),
            "p": ca.vertcat(self.model.varlist_independent.get_casadi_variables()),
            "ode": self.model.equations_differential,
        }

        # Tau variable is used to specify a length of iteration step externally, via tau variable
        self.tau: ca.MX = ca.MX.sym("tau")
        self.ode_system_tau: dict[str, ca.MX] = {
            "x": self.model.varlist_state.get_casadi_variables(),
            "p": ca.vertcat(
                self.tau, self.model.varlist_independent.get_casadi_variables()
            ),
            "ode": self.model.equations_differential * self.tau,
        }

        if self.model.DAE:
            self.ode_system["alg"] = self.model.equations_algebraic
            self.ode_system_tau["alg"] = self.model.equations_algebraic
            self.ode_system["z"] = self.model.varlist_algebraic.get_casadi_variables()
            self.ode_system_tau[
                "z"
            ] = self.model.varlist_algebraic.get_casadi_variables()

        if self.model.DAE:
            self.function_algebraic_equations = ca.Function(
                "alg_eq_sys",
                [
                    self.ode_system["z"],
                    ca.vertcat(self.ode_system["x"], self.ode_system["p"]),
                ],
                [self.ode_system["alg"]],
                ["x", "p"],
                ["alg"],
            )
            self.rootfinder = ca.rootfinder(
                "DAE_zf_init", "newton", self.function_algebraic_equations
            )

        self._set_integrator_settings(integrator_settings)

        self._setup_constraints_idas(use_idas_constraints)

        # TODO This integrator is not used so far...
        self.integrator: ca.Function = ca.integrator(
            "integrator",
            self.__integrator_name,
            self.ode_system,
            {"grid": self.time_grid_relative, "output_t0": False, "print_stats": True},
        )

        self.integrator_tau: ca.Function = ca.integrator(
            "integrator_tau",
            self.__integrator_name,
            self.ode_system_tau,
            self.__integrator_settings,
        )

        # This integrator is used to output values of algebraic variables at time 0
        # and should be run first to get algebraic variables at time 0 for whole simulation
        integrator_settings_with_output_t0 = copy.deepcopy(self.__integrator_settings)
        integrator_settings_with_output_t0["output_t0"] = True
        self.integrator_tau_with_t0: ca.Function = ca.integrator(
            "integrator_tau_with_t0",
            self.__integrator_name,
            self.ode_system_tau,
            integrator_settings_with_output_t0,
        )

        # This list is used for utility functions, like finding steady state
        self._guess_or_value_of_independent_variables: list[float] = []

        # Here all lists from above are initialized
        self._setup_variables()

        # .factory() method is very expensive so should be requested externally
        if simulate_jac:
            if self.model.DAE is True:
                integrator_tau_jac = self.integrator_tau.factory(
                    "integrator_tau_jacobian",
                    self.integrator_tau.name_in(),
                    ["xf", "qf", "zf", "rxf", "rqf", "rzf", "jac:xf:p"],
                )
            else:
                integrator_tau_jac = self.integrator_tau.factory(
                    "integrator_tau_jacobian",
                    self.integrator_tau.name_in(),
                    ["xf", "qf", "rxf", "rqf", "jac:xf:p"],
                )
            self.integrator_tau_jac: ca.Function = integrator_tau_jac

        self._reset_scaling()

        # This code is moved here, so this if statement shouldn't be called every simulation
        if self.model.DAE is True:
            if recalculate_algebraic:
                simulate = self._simulate_dae_calculate_algebraic
            else:
                simulate = self._simulate_dae
            simulate_jac_func = self._simulate_jac_dae
        else:
            simulate = self._simulate_ode
            simulate_jac_func = self._simulate_jac_ode

        self.simulate: Callable[[], dict[str, ca.DM | ca.MX]] = simulate
        self.simulate_jac: Callable[[], dict[str, ca.DM | ca.MX]] = simulate_jac_func

    def _reset_scaling(self) -> None:
        self.scaling: ca.DM = ca.DM.ones(self._independent_variables[0].size())

    def _setup_constraints_idas(self, use_idas_constraints: bool) -> None:
        """Holds a list of constraints for state and algebraic variables
        which can be used to constrain the solution of the idas
        to positive or negative numbers"""
        self._constraints_idas: list[float] = []
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

    def _setup_variables(self) -> None:
        """Setup all important lists for simulator"""
        num_time_steps = len(self.time_grid_relative) - 1
        independent_variables = []
        initial_algebraic = []
        initial_state = []

        for variable_name in self.model.varlist_all.keys():
            try:
                var = self.__input_variable_list[variable_name]
            except KeyError:
                continue

            if isinstance(var, VariableState):
                try:
                    initial_state.append(var.value[0])
                except Exception as e:
                    raise (BadVariableError(var)) from e

            elif isinstance(var, VariableAlgebraic):
                cast(VariableAlgebraic, var)
                initial_algebraic.append(var.guess)

            elif isinstance(var, VariableParameter):
                independent_variable = []
                independent_variable.extend(
                    [var.get_value_or_casadi()] * num_time_steps
                )
                self._guess_or_value_of_independent_variables.append(
                    var.get_value_or_guess()
                )
                independent_variables.append(independent_variable)
            elif isinstance(var, VariableControl):
                if isinstance(var, VariableControlPiecewiseConstant):
                    var_t0 = var.get_variable_at_time_relative(0)
                    self._guess_or_value_of_independent_variables.append(
                        var_t0.get_value_or_guess()
                    )
                    independent_variable = var.get_value_or_casadi(
                        self.time_grid_relative
                    )
                else:
                    independent_variable = []
                    self._guess_or_value_of_independent_variables.append(
                        var.get_value_or_guess()
                    )
                    independent_variable.extend(
                        [var.get_value_or_casadi()] * num_time_steps
                    )

                independent_variables.append(independent_variable)

        """ This nested list holds either a value or a casadi variable of
        each independent variable at every timestamp in time_grid.  First it has form:
        [[var1_t0, var1_t1 ...], [var2_t0, var2_t1 ...], [...]]
        Than it's reformed to:
        [[var1_t0, var2_t0 ...], [var1_t1, var2_t1 ...], [...]]
        And finally nested lists are changed to casadi.MX or DM vectors
        [ca.MX(var1_t0, var2_t0 ...), ca.MX(var1_t1, var2_t1 ...), [...]]
        """
        self._independent_variables: list[ca.MX | ca.DM] = list(
            map(list, zip(*independent_variables))
        )
        # List of values of State Variables at time 0
        self._initial_state: list[float] = initial_state

        # List of expected or recalculated values of Algebraic Variables at time 0
        # This list is further used in calculations
        self._initial_algebraic: list[float] = initial_algebraic

        # Groups nested lists by time_stamp
        self._initial_algebraic_original = copy.deepcopy(self._initial_algebraic)

        # Transforms nested lists in ca.MX or ca.DM array
        for index, column in enumerate(self._independent_variables):
            casadi_mx = ca.vcat(column)
            self._independent_variables[index] = casadi_mx

    def _set_integrator_settings(self, provided_settings: dict | None) -> None:
        """Sane default settings for integrators"""
        if provided_settings is not None:
            integrator_settings = provided_settings
        elif self.__integrator_name == "idas":
            integrator_settings = {
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
            integrator_settings = {
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
            integrator_settings = {
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

        self.__integrator_settings = integrator_settings

    def calculate_steady_state(self) -> dict[str, ca.DM]:
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

    def calculate_algebraic_initials(
        self, *, apply_intials: bool = False, analyze: bool = False
    ) -> None:
        if self.model.DAE:
            function = ca.Function(
                "eq_sys",
                [self.ode_system["x"], self.ode_system["z"], self.ode_system["p"]],
                [self.ode_system["ode"], self.ode_system["alg"]],
                ["x", "z", "p"],
                ["ode", "alg"],
            )

            res = self.rootfinder(
                self._initial_algebraic_original,
                ca.vertcat(
                    self._initial_state, self._guess_or_value_of_independent_variables
                ),
            )

            residual_original = function(
                x=self._initial_state,
                z=self._initial_algebraic_original,
                p=self._guess_or_value_of_independent_variables,
            )
            residual_calculated = function(
                x=self._initial_state,
                z=res,
                p=self._guess_or_value_of_independent_variables,
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

                # import pandas as pd

                # jac_func = function.jacobian()
                # jac = jac_func(
                #     x=self._initial_state,
                #     z=res,
                #     p=self._guess_or_value_of_independent_variables,
                # )

                # row_names = str(self.ode_system["ode"]).split()
                # row_names.extend(str(self.ode_system["alg"]).split())

                # col_names = str(self.ode_system["x"]).split()
                # col_names.extend(str(self.ode_system["z"]).split())
                # col_names.extend(str(self.ode_system["p"]).split())

                # jac = jac["jac"]
                # df = pd.DataFrame(jac.toarray(), index=row_names, columns=col_names)

            if apply_intials:
                self.logger.debug("Fixed algebraic intials")
                self._initial_algebraic = res

    def analyze_WIP(
        self, state_value: list[float] = None
    ) -> list[dict[str, ca.DM] | list[float]]:
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
            x=self._initial_state,
            z=self._initial_algebraic,
            p=self._independent_variables[0],
        )
        jacobian = function.factory(
            "jac_alg",
            function.name_in(),
            ["jac:alg:z", "jac:alg:x", "jac:ode:x", "jac:ode:z"],
        )
        check_jacobian = jacobian(
            x=self._initial_state,
            z=self._initial_algebraic,
            p=self._independent_variables[0],
        )

        check_alg = algebraic_eqsys(
            x=self._initial_state,
            z=self._initial_algebraic,
            p=self._independent_variables[0],
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
                self._initial_algebraic,
                ca.vertcat(state_value, self._independent_variables[0]),
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

    def _simulate_jac_dae(self) -> dict[str, ca.DM | ca.MX]:
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

    def _simulate_jac_ode(self) -> dict[str, ca.DM | ca.MX]:
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

    def _simulate_dae_calculate_algebraic(self) -> dict[str, ca.DM | ca.MX]:
        """Return dictionary with results "xf" - state,
        "zf" - algebraic
        """
        prev_time_step = 0
        res_states = []
        res_algebraic = []
        x_init = self._initial_state

        alg_init = self.rootfinder(
            self._initial_algebraic_original,
            ca.vertcat(
                x_init,
                self._independent_variables[0] * self.scaling,
            ),
        )

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

    def _simulate_t0(self) -> dict[str, ca.DM | ca.MX]:
        """Return dictionary with results "xf0" and "zf0" for state and
        algebraic variables at time 0"""
        prev_time_step = self.time_grid_relative[1]

        if self.model.DAE:
            res_integration = self.integrator_tau_with_t0(
                x0=self._initial_state,
                z0=self._initial_algebraic,
                p=ca.vertcat(
                    prev_time_step,
                    self._independent_variables[0] * self.scaling,
                ),
            )
        else:
            res_integration = self.integrator_tau_with_t0(
                x0=self._initial_state,
                p=ca.vertcat(
                    prev_time_step,
                    self._independent_variables[0] * self.scaling,
                ),
            )

        init_states = res_integration["xf"][:, 0]

        res = {"xf0": init_states}
        if self.model.DAE:
            init_algebraic = res_integration["zf"][:, 0]
            res["zf0"] = init_algebraic

        return res

    def _simulate_dae(self) -> dict[str, ca.DM | ca.MX]:
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

    def _simulate_ode(self) -> dict[str, ca.DM | ca.MX]:
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

    def generate_exp_data(
        self,
        algebraic: bool = False,
        recalculate_algebraic: bool = True,
        unfixed_variables: list[float] | np.ndarray | None = None,
    ) -> VariableList:
        """Runs simulation and returns results in VariableList class."""
        variables = VariableList()

        if recalculate_algebraic and self.model.DAE:
            self.calculate_algebraic_initials(apply_intials=True)

        result_simulation = self.simulate()
        result_initial = self._simulate_t0()
        if not algebraic or not self.model.DAE:
            result_varlist = [copy.deepcopy(self.model.varlist_state)]
            res_array = ca.horzcat(result_initial["xf0"], result_simulation["xf"])
        else:
            result_varlist = [
                copy.deepcopy(self.model.varlist_state),
                copy.deepcopy(self.model.varlist_algebraic),
            ]
            res_array = ca.vertcat(
                ca.horzcat(
                    result_initial["xf0"],
                    result_simulation["xf"],
                ),
                ca.horzcat(
                    result_initial["zf0"],
                    result_simulation["zf"],
                ),
            )

        if not isinstance(res_array, ca.DM):
            if unfixed_variables is None:
                raise ValueError("You need to supply values for unfixed variables")
            else:
                function = ca.Function("f", ca.symvar(res_array), [res_array])
                res_array = function(*unfixed_variables)

        shift_by = 0

        for variable_list in result_varlist:
            for count, var in enumerate(variable_list.values()):
                var.casadi_var = None
                new_var = copy.deepcopy(var)

                value = res_array[count + shift_by, :]
                # value is of ca.DM type and data is nested in first array
                value = value.toarray()[0]

                new_var.set_dataframe_from_value_and_time(
                    value, self.time_grid_relative, self.origin_ts
                )
                new_var.ignore_plotting = self.__input_variable_list[
                    var.name
                ].ignore_plotting

                variables.add_variable(new_var)
            shift_by = count + 1

        return variables

    def setup_time_grid(self, time_grid: np.ndarray) -> None:
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
            "Timegrid modified: \n self.timegrid \n {0} \n".format(
                self.time_grid_relative
            )
        )


class SimulatorNLE:

    supported_solvers: list[str] = ["ipopt", "rootfinder"]

    def __init__(
        self,
        model: Model,
        variable_list: VariableList,
        solver_settings: dict | None = None,
        solver_name: str = "rootfinder",
        *,
        use_bounds: bool = True,
    ):
        self.model: Model = model
        if solver_name not in self.supported_solvers:
            raise TypeError(
                f"Provided integrator name {solver_name} is not supported. Only theese are: {self.supported_solvers}."
            )
        if solver_name == "ipopt":
            print("IPOPT option for SimulatorNLE is not fully supported, use only if you know what you're doing")

        self.__solver_name: str = solver_name
        self.__input_variable_list: VariableList = copy.deepcopy(variable_list)

        self._set_solver_settings(solver_settings)

        self._setup_variables()
        self._reset_scaling()

        if self.__solver_name == "rootfinder":
            if use_bounds:
                self.solver_settings["constraints"] = self._rootfinder_bounds

            self.function: ca.Function = ca.Function(
                "f",
                [
                    self.model.varlist_algebraic.get_casadi_variables(),
                    self.model.varlist_independent.get_casadi_variables(),
                ],
                [self.model.equations_algebraic],
                ["x0", "p"],
                ["x"],
            )
            self.simulator: ca.Function = ca.rootfinder(
                "s", "nlpsol", self.function, self.solver_settings
            )
            self.call_arg: dict = {
                "x0": ca.DM(self._guess),
                "p": self._independent_variables * self.scaling,
            }
        elif self.__solver_name == "ipopt":
            self.simulator = ca.nlpsol(
                "solver",
                "ipopt",
                {
                    "x": self.model.varlist_algebraic.get_casadi_variables(),
                    "p": self.model.varlist_independent.get_casadi_variables(),
                    "g": self.model.equations_algebraic,
                    "f": (ca.sum1(self.model.equations_algebraic) ** 2),
                },
                self.solver_settings,
            )
            self.call_arg = {
                "x0": ca.DM(self._guess),
                "p": self._independent_variables * self.scaling,
                "lbg": 0,
                "ubg": 0,
            }
            if use_bounds:
                self.call_arg["lbx"] = self._lower_bound
                self.call_arg["ubx"] = self._upper_bound

        self.jacobian: ca.Function = self.simulator.jacobian()

    def _set_solver_settings(self, provided_settings: dict | None) -> None:
        """Set default settings, if None are provided"""
        if provided_settings is not None:
            solver_settings = provided_settings
        elif self.__solver_name == "rootfinder":
            solver_settings = {
                "nlpsol": "ipopt",
                "verbose": False,
                "print_in": False,
                "print_out": False,
                "expand": True,
                "nlpsol_options": {
                    "ipopt.hessian_approximation": "limited-memory",
                    "ipopt.max_iter": 300,
                    "ipopt.print_level": 0,
                    "print_time": False,
                },
            }
        elif self.__solver_name == "ipopt":
            solver_settings = {
                "verbose": False,
                "print_in": False,
                "print_out": False,
                "print_time": False,
                "expand": True,
                "ipopt": {
                    "hessian_approximation": "limited-memory",
                    "max_iter": 300,
                    "print_level": 0,
                },
            }

        self.solver_settings: dict = solver_settings

    def _setup_variables(self) -> None:
        guess = []
        lower_bound = []
        upper_bound = []
        rootfinder_bounds = []
        independent_variables = []
        for variable_name in self.model.varlist_all.keys():
            try:
                var = self.__input_variable_list[variable_name]
            except KeyError:
                continue

            if isinstance(var, VariableAlgebraic):
                guess.append(var.guess)
                if var.lower_bound is None:
                    lower_bound.append(-ca.inf)
                else:
                    lower_bound.append(var.lower_bound)
                if var.upper_bound is None:
                    upper_bound.append(ca.inf)
                else:
                    upper_bound.append(var.upper_bound)
            elif isinstance(var, VariableConstant):
                pass
            elif isinstance(var, (VariableControl, VariableParameter)):
                independent_variables.append(var.get_value_or_casadi())
            else:
                raise TypeError(f"{type(var)} is not supported")

        self._lower_bound: list[float] = lower_bound
        self._upper_bound: list[float] = upper_bound

        for lower_bound_i, upper_bound_i in zip(self._lower_bound, self._upper_bound):
            rootfinder_bound = 0
            if lower_bound_i == 0:
                rootfinder_bound = 1
            elif lower_bound_i > 0:
                rootfinder_bound = 2
            elif upper_bound_i == 0:
                rootfinder_bound = -1
            elif upper_bound_i < 0:
                rootfinder_bound = -2

            rootfinder_bounds.append(rootfinder_bound)

        self._guess: list[float] = guess
        self._rootfinder_bounds: list[int] = rootfinder_bounds
        self._independent_variables: ca.MX | ca.DM = ca.vcat(independent_variables)

        if isinstance(self._independent_variables, ca.MX):
            self.contains_unfixed = True
        elif isinstance(self._independent_variables, ca.DM):
            self.contains_unfixed = False
        else:
            raise NotImplementedError

    def generate_exp_data(self, unfixed_variables: list[float] = None) -> VariableList:
        res_array = self.simulate_sym_unfixed(unfixed_variables)

        variables = VariableList()

        for count, var in enumerate(self.model.varlist_algebraic.values()):
            new_var = copy.deepcopy(var)
            new_var.casadi_var = None
            new_var.lower_bound = self._lower_bound[count]
            new_var.upper_bound = self._upper_bound[count]
            new_var.set_dataframe_from_value_and_time([float(res_array[count])], [0])
            new_var.ignore_plotting = self.__input_variable_list[
                var.name
            ].ignore_plotting
            variables.add_variable(new_var)

        return variables

    def _reset_scaling(self) -> None:
        self.scaling: ca.DM = ca.DM.ones(self._independent_variables.size())

    def simulate_sym_unfixed(self, unfixed_variables: list[float] = None) -> ca.DM:
        """This is slower version of simulate_sym but it allows user to supply values
        for unfixed variables"""
        self.call_arg["p"] = self._independent_variables * self.scaling
        res_array = self.simulator.call(self.call_arg)["x"]

        if not isinstance(res_array, ca.DM):
            if unfixed_variables is None:
                raise ValueError("You need to supply values for unfixed variables")
            else:
                function = ca.Function("f", ca.symvar(res_array), [res_array])
                res_array = function(*unfixed_variables)

        return res_array

    def simulate_sym(self) -> dict[str, ca.MX | ca.DM]:
        self.call_arg["p"] = self._independent_variables * self.scaling

        res = self.simulator.call(self.call_arg)
        return res

    def calculate_jac(self) -> dict[str, ca.MX | ca.DM]:
        self.call_arg["p"] = self._independent_variables * self.scaling

        res = self.jacobian.call(self.call_arg)
        return res

    def get_variance_array(self) -> np.ndarray:
        variance_list = []
        for variable_name in self.model.varlist_all.keys():
            try:
                var = self.__input_variable_list[variable_name]
            except KeyError:
                continue

            if isinstance(var, VariableAlgebraic):
                variance_list.append(var.variance)

        variance = np.array(variance_list)
        return variance
