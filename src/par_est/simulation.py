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
    VariableConstant
)


class Simulator(object):
    def __init__(
        self,
        model: Model,
        time_grid,
        variable_list: VariableList,
        integrator_name="idas",
        integrator_settings=None,
    ):

        self.logger = logging.getLogger(__name__)
        self.logger.debug(
            "Creating Simulator object: \n timegrid \n {0} \n".format(time_grid)
        )
        self.__input_variable_list = copy.deepcopy(variable_list)
        self.model = model
        self.tau = ca.MX.sym("tau")
        self.scaling = None
        self.integrator_settings = None
        self.integrator_name = None
        self.time_grid = time_grid

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

        steps_integrator = len(self.time_grid) - 1
        self.integrator_full = self.integrator_tau.mapaccum(
            "integrator_full", steps_integrator
        )

        if self.model.DAE:
            self.integrator_tau_jacobian = self.integrator_tau.factory(
                "integrator_tau_jacobian",
                self.integrator_full.name_in(),
                ["xf", "zf", "jac:xf:p"],
            )
        else:
            self.integrator_tau_jacobian = self.integrator_tau.factory(
                "integrator_tau_jacobian",
                self.integrator_full.name_in(),
                ["xf", "qf", "zf", "rxf", "rqf", "rzf", "jac:xf:p"],
            )

        self.integrator_full_jacobian = self.integrator_tau_jacobian.mapaccum(
            "jacobian", steps_integrator
        )

        # Arrays needed to initialize integrator.
        self._variables = []
        self._initial_state = []
        self._initial_algebraic = []

        for var in self.__input_variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    self._initial_state.append(var.value.value[0])
                elif isinstance(var, VariableAlgebraic):
                    self._initial_algebraic.append(var.guess)
                elif isinstance(var, VariableConstant):
                    pass
                else:
                    if var.fixed:
                        self._variables.append(var.value)
                    else:
                        self._variables.append(var.casadi_var)

        self._variables = ca.vcat(self._variables)
        self._reset_scaling()

    def _reset_scaling(self):
        self.scaling = ca.DM.ones(self._variables.size())

    def analyze_WIP(self, state_value=None):
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
        jacobian = function.factory("jac_alg", function.name_in(), ["jac:alg:z"])
        check_jacobian = jacobian(
            x=self._initial_state, z=self._initial_algebraic, p=self._variables
        )

        check_alg = algebraic_eqsys(
            x=self._initial_state, z=self._initial_algebraic, p=self._variables
        )

        # should fail by DAE index > 1
        ca.inv(check_jacobian["jac_alg_z"])

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

    def simulate(self, derivatives=False):
        """Return dictionary with results "xf" - state,
        "zf" - algebraic, "jac_xf_p" - derivatives.
        """
        map_num = len(self.time_grid) - 1
        initial_independent = ca.vertcat(
            ca.horzcat(*(self.time_grid[1:] - self.time_grid[:-1])),
            ca.repmat(self._variables * self.scaling, 1, map_num),
        )

        self.logger.debug(
            "Simulating: \n Initial States x0 \n {} \n Independent Variables p \n {} \n".format(
                self._initial_state, initial_independent
            )
        )
        if self.model.DAE:
            self.logger.debug(
                "Initial Algebraic z0 \n {} \n".format(self._initial_algebraic)
            )

        if not derivatives:
            if self.model.DAE:
                result_integration = self.integrator_full(
                    x0=self._initial_state,
                    z0=self._initial_algebraic,
                    p=initial_independent,
                )
            else:
                result_integration = self.integrator_full(
                    x0=self._initial_state, p=initial_independent
                )
        else:
            if self.model.DAE:
                result_integration = self.integrator_full_jacobian(
                    x0=self._initial_state,
                    z0=self._initial_algebraic,
                    p=initial_independent,
                )
            else:
                result_integration = self.integrator_full_jacobian(
                    x0=self._initial_state,
                    p=initial_independent,
                )

        return result_integration

    def generate_exp_data(self, algebraic=False):
        """ Runs simulation and returns results in VariableList class."""
        variables = VariableList()
        result_simulation = self.simulate()
        if not algebraic:
            result_varlist = [self.model.varlist_state]
            res_array = result_simulation["xf"]
        else:
            result_varlist = [self.model.varlist_state, self.model.varlist_algebraic]
            res_array = ca.vertcat(result_simulation["xf"], result_simulation["zf"])

        convert_to_numpy = False
        if isinstance(res_array, ca.DM):
            convert_to_numpy = True

        shift_by = 0
        for variable_list in result_varlist:
            for count, var in enumerate(variable_list.values()):
                new_var = copy.deepcopy(var)
                new_var.value = ExperimentData()
                if convert_to_numpy:
                    new_var.value.time = self.time_grid
                    new_var.value.value = res_array[count + shift_by, :].toarray()
                    if isinstance(var, VariableAlgebraic):
                        value_time_zero = var.guess
                    elif isinstance(var, VariableState):
                        value_time_zero = self.__input_variable_list[var.name].value.value[0],
                    else:
                        raise(NotImplementedError)

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


class SimulatorNLE:
    def __init__(self, model: Model, variable_list: VariableList):
        self.model = model
        self.__input_variable_list = copy.deepcopy(variable_list)

        self._variables = []
        self._guess = []

        # Fügt Variablen hinzu und sotiert nach festen/zu schätzenden Werten
        for var in self.__input_variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableAlgebraic):
                    self._guess.append(var.guess)
                else:
                    if var.fixed:
                        self._variables.append(var.value)
                    else:
                        self._variables.append(var.casadi_var)

        self._variables = ca.vcat(self._variables)

        self.function = ca.Function(
            "f",
            [
                self.model.varlist_algebraic.get_casadi_var(),
                self.model.varlist_independent.get_casadi_var(),
            ],
            [self.model.equations_algebraic],
            ["x0", "p"],
            ["r"],
        )
        self._reset_scaling()

    def generate_exp_data(self):
        opts = {}
        opts_ipopt = {}
        # opts["strategy"] = "linesearch"
        opts["nlpsol"] = "ipopt"
        opts_ipopt["print_time"] = False
        opts_ipopt["ipopt.print_level"] = 0
        opts["nlpsol_options"] = opts_ipopt
        # opts["abstol"] = 1e-14
        # opts["constraints"] = [2, -2]

        sim = ca.rootfinder("s", "nlpsol", self.function, opts)
        # breakpoint()

        res_array = sim(self._guess, self._variables)

        variables = VariableList()

        for variable_list in [self.model.varlist_algebraic]:
            for count, var in enumerate(variable_list.values()):
                new_var = copy.deepcopy(var)
                new_var.value = ExperimentData()
                new_var.value.time = 0
                new_var.value.value = res_array[count]
                variables.add_variable(new_var)

        return variables

    def _reset_scaling(self):
        self.scaling = ca.DM.ones(self._variables.size())

    def simulate_sym(self):
        opts = {}
        opts_ipopt = {}
        # opts["strategy"] = "linesearch"
        opts["nlpsol"] = "ipopt"
        opts_ipopt["print_time"] = False
        opts_ipopt["ipopt.print_level"] = 0
        opts["nlpsol_options"] = opts_ipopt
        # opts["strategy"] = "linesearch"
        # opts["abstol"] = 1e-14
        # opts["constraints"] = [2, -2]
        # opts["verbose"] = True
        # opts["print_in"] = True
        # opts["print_out"] = True
        # opts["print_stats"] = True

        sim = ca.rootfinder("s", "nlpsol", self.function, opts)
        arg = {"x0": ca.DM(self._guess), "p": self._variables}

        # res = sim(x0=self._guess, p=self._variables)
        res = sim.call(arg)
        return res
