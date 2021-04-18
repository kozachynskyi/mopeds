import copy
import casadi as ca
import numpy as np
import logging

from par_est import (
    VariableList,
    Model,
    ExperimentData,
    Variable,
    VariableState,
    VariableAlgebraic,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableParameter,
    VariableConstant,
    BadVariableError,
)


class Simulator(object):
    def __init__(
        self,
        model: Model,
        input_time_grid,
        variable_list: VariableList,
        integrator_name="idas",
        integrator_settings=None,
    ):

        self.logger = logging.getLogger(__name__)
        self.logger.debug(
            "Creating Simulator object: \n timegrid \n {0} \n".format(input_time_grid)
        )
        self.__input_variable_list = copy.deepcopy(variable_list)
        self.model = model
        self.tau = ca.MX.sym("tau")
        self.scaling = None
        self.integrator_settings = None
        self.integrator_name = None
        self.setup_time_grid(input_time_grid)
        self.logger.debug(
            "Timegrid modified: \n self.timegrid \n {0} \n".format(self.time_grid)
        )

        if self.model.equations_algebraic is None:
            self.model.DAE = False
        else:
            self.model.DAE = True

        self.ode_system = {
            "x": self.model.varlist_state.get_casadi_var(),
            "p": ca.vertcat(self.model.varlist_independent.get_casadi_var()),
            "ode": self.model.equations_differential,
        }

        self.ode_system_tau = {
            "x": self.model.varlist_state.get_casadi_var(),
            "p": ca.vertcat(self.tau, self.model.varlist_independent.get_casadi_var()),
            "ode": self.model.equations_differential * self.tau,
        }

        if self.model.DAE:
            self.ode_system["alg"] = self.model.equations_algebraic
            self.ode_system_tau["alg"] = self.model.equations_algebraic
            self.ode_system["z"] = self.model.varlist_algebraic.get_casadi_var()
            self.ode_system_tau["z"] = self.model.varlist_algebraic.get_casadi_var()

        if integrator_name == "idas":
            self.integrator_name = "idas"
        elif integrator_name == "cvodes":
            self.integrator_name = "cvodes"
        else:
            self.integrator_name = "collocation"

        if integrator_settings is not None:
            self.integrator_settings = integrator_settings
        else:
            if self.integrator_name == "idas":
                self.integrator_settings = {
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
            elif self.integrator_name == "cvodes":
                self.integrator_settings = {
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
            else:
                self.integrator_settings = {
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

        # TODO This integrator is not used so far...
        self.integrator = ca.integrator(
            "integrator",
            "idas",
            self.ode_system,
            {"grid": self.time_grid, "output_t0": False, "print_stats": True},
        )

        self.integrator_tau = ca.integrator(
            "integrator_tau",
            self.integrator_name,
            self.ode_system_tau,
            self.integrator_settings,
        )

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

        # Arrays needed to initialize integrator.
        num_time_steps = len(self.time_grid) - 1
        self._variables = []
        self._initial_state = []
        self._initial_algebraic = []
        self._initial_algebraic_original = []
        self._variables_with_guess = []

        for var in self.__input_variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    try:
                        self._initial_state.append(var.value.value[0])
                    except Exception as e:
                        raise (BadVariableError(var)) from e
                elif isinstance(var, VariableAlgebraic):
                    self._initial_algebraic.append(var.guess)
                elif isinstance(var, VariableConstant):
                    pass
                elif isinstance(var, VariableParameter):
                    independent_variable = []
                    independent_variable.extend([var.get_value_based_on_fixed()] * num_time_steps)
                    if var.fixed:
                        self._variables_with_guess.append(var.value)
                    else:
                        self._variables_with_guess.append(var.guess)
                    self._variables.append(independent_variable)
                elif isinstance(var, VariableControl):
                    independent_variable = []
                    if isinstance(var, VariableControlPiecewiseConstant):
                        var_t0 = var.var_at_time(0)
                        if var_t0.fixed:
                            self._variables_with_guess.append(var_t0.value)
                        else:
                            self._variables_with_guess.append(var_t0.guess)

                        last_unfixed_variable = None
                        for time_stamp in self.time_grid:
                            var_at_timestamp = var.var_at_time(time_stamp)
                            # This if statement is required for OED in order to use casadi_var from previous step, if it was already used. Without it, control variable will be fixed to some value for given timestep
                            if var_at_timestamp.fixed:
                                if last_unfixed_variable is None:
                                    independent_variable.append(var_at_timestamp.get_value_based_on_fixed())
                                else:
                                    independent_variable.append(last_unfixed_variable.get_value_based_on_fixed())
                            else:
                                last_unfixed_variable = var_at_timestamp
                                independent_variable.append(last_unfixed_variable.get_value_based_on_fixed())

                    else:
                        if var.fixed:
                            self._variables_with_guess.append(var.value)
                        else:
                            self._variables_with_guess.append(var.guess)
                        independent_variable.extend([var.get_value_based_on_fixed()] * num_time_steps)
                    self._variables.append(independent_variable)


        self._variables = list(map(list, zip(*self._variables)))
        self._initial_algebraic_original = copy.deepcopy(self._initial_algebraic)

        # Transforms nested python list in ca.MX array
        for index, column in enumerate(self._variables):
            casadi_mx = ca.vcat(column)
            self._variables[index] = casadi_mx
        self._reset_scaling()

        if self.model.DAE is True:
            self.simulate = self._simulate_dae
            self.simulate_jac = self._simulate_jac_dae
        else:
            self.simulate = self._simulate_ode
            self.simulate_jac = self._simulate_jac_ode

    def _reset_scaling(self):
        self.scaling = ca.DM.ones(self._variables[0].size())

    def calculate_algebraic_initials(self, *, apply_intials=False):
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
        res = rf(
            self._initial_algebraic_original,
            ca.vertcat(self._initial_state, self._variables_with_guess),
        )

        residual_original = function(
            x=self._initial_state,
            z=self._initial_algebraic_original,
            p=self._variables_with_guess,
        )
        residual_calculated = function(
            x=self._initial_state, z=res, p=self._variables_with_guess
        )

        if apply_intials:
            residual_sum_original = ca.sum1(residual_original["alg"])
            residual_sum_calculated = ca.sum1(residual_calculated["alg"])
            self.logger.debug(
                f"Fixed algebraic intials. Residual before {residual_sum_original}, after {residual_sum_calculated}."
            )
            self._initial_algebraic = res

    def analyze_WIP(self, state_value=None):
        import par_est.tools as tools

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

        check_initials = function(
            x=self._initial_state, z=self._initial_algebraic, p=self._variables
        )
        jacobian = function.factory(
            "jac_alg",
            function.name_in(),
            ["jac:alg:z", "jac:alg:x", "jac:ode:x", "jac:ode:z"],
        )
        check_jacobian = jacobian(
            x=self._initial_state, z=self._initial_algebraic, p=self._variables
        )

        check_alg = algebraic_eqsys(
            x=self._initial_state, z=self._initial_algebraic, p=self._variables
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
            res = rf(self._initial_algebraic, ca.vertcat(state_value, self._variables))
        else:
            res = rf(
                self._initial_algebraic,
                ca.vertcat(self._initial_state, self._variables),
            )

        check_alg = function(x=self._initial_state, z=res, p=self._variables)
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

        for time_step, independent_variables in zip(self.time_grid[1:], self._variables):
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

        for time_step, independent_variables in zip(self.time_grid[1:], self._variables):
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

        for time_step, independent_variables in zip(self.time_grid[1:], self._variables):
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

        for time_step, independent_variables in zip(self.time_grid[1:], self._variables):
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

    def generate_exp_data(self, algebraic=False):
        """ Runs simulation and returns results in VariableList class."""
        variables = VariableList()
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

        convert_to_numpy = False
        if isinstance(res_array, ca.DM):
            convert_to_numpy = True

        shift_by = 0
        for variable_list in result_varlist:
            for count, var in enumerate(variable_list.values()):
                var.casadi_var = None
                new_var = copy.deepcopy(var)
                new_var.value = ExperimentData()
                if convert_to_numpy:
                    new_var.value.time = self.time_grid
                    new_var.value.value = res_array[count + shift_by, :].toarray()
                    if isinstance(var, VariableAlgebraic):
                        value_time_zero = var.guess
                    elif isinstance(var, VariableState):
                        value_time_zero = (
                            self.__input_variable_list[var.name].value.value[0],
                        )
                    else:
                        raise (NotImplementedError)

                    new_var.value.value = np.insert(
                        new_var.value.value,
                        0,
                        value_time_zero,
                    )
                else:
                    new_var.value.time = self.time_grid[1:]
                    new_var.value.value = res_array[count + shift_by, :]

                variables.add_variable(new_var)
            shift_by = count + 1

        return variables

    def setup_time_grid(self, time_grid):
        """ Time_grid provided by used may not take into account piecewise controls.
        Thus it might be needed to expand a time grid. """
        for var in self.__input_variable_list.values():
            if isinstance(var, VariableControlPiecewiseConstant):
                time_grid = np.append(time_grid, var.time)
        self.time_grid = np.unique(time_grid)
