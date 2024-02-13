from __future__ import annotations

import copy
import logging
from typing import Callable, cast
from numpy.typing import ArrayLike
from warnings import warn

import casadi as ca
import numpy as np

from mopeds import (
    BadVariableError,
    Model,
    VariableAlgebraic,
    VariableConstant,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableList,
    VariableParameter,
    VariableState,
    _ACADOS_SUPPORT,
    get_options,
    _consistent_scaling_decorator,
)

if _ACADOS_SUPPORT:
    from acados_template import AcadosModel
    from mopeds import casados_integrator


class SimulatorNLE:

    supported_solvers: list[str] = ["ipopt", "rootfinder"]

    def __init__(
        self,
        model: Model,
        variable_list: VariableList,
        solver_settings: dict | None = None,
        solver_name: str = "rootfinder",
        *,
        use_bounds: bool = None,
    ):
        self._created_with_options = get_options()
        if use_bounds is not None:
            warn("use_bounds argument is ignored", FutureWarning, 2)

        self.model: Model = model
        if solver_name not in self.supported_solvers:
            raise TypeError(
                f"Provided integrator name {solver_name} is not supported. Only theese are: {self.supported_solvers}."
            )
        if solver_name == "ipopt":
            self._call_simulator = self.__call_simulator_ipopt
        elif solver_name == "rootfinder":
            self._call_simulator = self.__call_simulator_rootfinder

        self._solver_name: str = solver_name
        self._input_variable_list: VariableList = copy.deepcopy(variable_list)
        self.model.subsitute_casadi_symbols(self._input_variable_list)

        if solver_settings is not None:
            self.solver_settings: dict = solver_settings
        else:
            self.solver_settings = self.get_default_simulator_settings()

        self._setup_variables()

        scaled_equations = ca.substitute(self.model.equations_algebraic, self._input_variable_list.get_casadi_variables(), self._input_variable_list.get_scaled_casadi_variables())
        self._model_equations = ca.cse(scaled_equations)

        if self._solver_name == "rootfinder":
            self.function: ca.Function = ca.Function(
                "f",
                [
                    self.model.varlist_algebraic(self._input_variable_list).get_casadi_variables(),
                    self.model.varlist_independent(self._input_variable_list).get_casadi_variables(),
                ],
                [self._model_equations],
                ["x0", "p"],
                ["x"],
            )
            self.simulator: ca.Function = ca.rootfinder(
                "s", "nlpsol", self.function, self.solver_settings
            )
            self.call_arg: dict = {
                "x0": ca.DM(self._guess),
                "p": self._independent_variables,
            }
        elif self._solver_name == "ipopt":
            self.simulator = ca.nlpsol(
                "solver",
                "ipopt",
                {
                    "x": self.model.varlist_algebraic(self._input_variable_list).get_casadi_variables(),
                    "p": self.model.varlist_independent(self._input_variable_list).get_casadi_variables(),
                    "g": self._model_equations,
                    "f": (ca.sum1(self._model_equations) ** 2),
                },
                self.solver_settings,
            )
            self.call_arg = {
                "x0": ca.DM(self._guess),
                "p": self._independent_variables,
                "lbg": 0,
                "ubg": 0,
            }
            self.call_arg["lbx"] = self._lower_bound
            self.call_arg["ubx"] = self._upper_bound

        self.jacobian: ca.Function = self.simulator.jacobian()

    def get_default_simulator_settings(self) -> None:
        """Set default settings, if None are provided"""
        if self._solver_name == "rootfinder":
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
        elif self._solver_name == "ipopt":
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

        return solver_settings

    @_consistent_scaling_decorator
    def _setup_variables(self) -> None:
        mapping_independent_variables = {}
        mapping_algebraic_variables = {}
        index_algebraic = 0
        index_independent = 0

        guess = []
        lower_bound = []
        upper_bound = []
        independent_variables = []
        for var in self.model.varlist(self._input_variable_list).values():
            if isinstance(var, VariableAlgebraic):
                mapping_algebraic_variables[var.name] = index_algebraic
                index_algebraic += 1
                guess.append(var.scale_from_original(var.guess))
                if var.lower_bound is None:
                    lower_bound.append(-ca.inf)
                else:
                    lower_bound.append(var.scale_from_original(var.lower_bound))
                if var.upper_bound is None:
                    upper_bound.append(ca.inf)
                else:
                    upper_bound.append(var.scale_from_original(var.upper_bound))
            elif isinstance(var, VariableConstant):
                pass
            elif isinstance(var, (VariableControl, VariableParameter)):
                mapping_independent_variables[var.name] = index_independent
                index_independent += 1
                independent_variables.append(var.get_value_or_casadi())
            else:
                raise TypeError(f"{type(var)} is not supported")

        self.mapping_independent_variables: dict[
            str, int
        ] = mapping_independent_variables
        self.mapping_algebraic_variables: dict[str, int] = mapping_algebraic_variables
        self._lower_bound: list[float] = lower_bound
        self._upper_bound: list[float] = upper_bound

        self._guess: list[float] = guess
        self._independent_variables: ca.MX | ca.DM = ca.vcat(independent_variables)

        if isinstance(self._independent_variables, ca.MX):
            self.contains_unfixed = True
        elif isinstance(self._independent_variables, ca.DM):
            self.contains_unfixed = False
        else:
            raise NotImplementedError

    @_consistent_scaling_decorator
    def change_independent_variables(self, ind_variables: dict[str, float]):
        """Use this method to change either Controls or Parameters of the Simulation. ind_variables is a dictionary
        with VariableNames as dict.keys(), and their respective values, as dict.values(). Example:
        {"e0_T": 373, "e0_p": 1e5}"""
        if self.contains_unfixed:
            raise NotImplementedError(
                "All variables should be fixed, to use this method"
            )
        for var_name, var_value in ind_variables.items():
            var = self._input_variable_list[var_name]
            index_var = self.mapping_independent_variables[var_name]
            self._independent_variables[index_var] = var.scale_from_original(var_value)

    @_consistent_scaling_decorator
    def generate_exp_data(self, unfixed_variables: dict[str, float] = None) -> VariableList:
        res_array = self.simulate_sym_unfixed(unfixed_variables)["x"].toarray()

        variables = VariableList()


        for var_name, res in zip(self.model.varlist_algebraic(self._input_variable_list).keys(), res_array):
            var = self._input_variable_list[var_name]
            new_var = copy.deepcopy(var)
            new_var.casadi_var = None

            value = var.scale_to_original(res)
            new_var.value = float(value[0])
            variables.add_variable(new_var)

        return variables

    def __call_simulator_rootfinder(self) -> ca.DM:
        """This method is needed to raise an error, if ipopt simulator fails to converge"""
        res = self.simulator.call(self.call_arg)
        return res

    def __call_simulator_ipopt(self) -> ca.DM:
        """This method is needed to raise an error, if ipopt simulator fails to converge"""
        res = self.simulator.call(self.call_arg)

        if isinstance(res["x"], ca.DM):
            if not self.simulator.stats()["success"]:
                raise ValueError(f"IPOPT failed as NLE solver:\n{self.simulator.stats()}")
        return res

    @_consistent_scaling_decorator
    def simulate_sym_unfixed(self, unfixed_variables: dict[str, float] = None) -> ca.DM:
        """This is slower version of simulate_sym but it allows user to supply values
        for unfixed variables"""
        self.call_arg["p"] = self._independent_variables
        res_array = self._call_simulator()["x"]

        if not isinstance(res_array, ca.DM):
            if unfixed_variables is None:
                raise ValueError("You need to supply values for unfixed variables")
            else:
                unfixed_symbols = ca.symvar(res_array)
                values = []
                for symbol in unfixed_symbols:
                    var_name = symbol.name()
                    values.append(self._input_variable_list[var_name].scale_from_original(unfixed_variables[var_name]))

                function = ca.Function("f", ca.symvar(*unfixed_symbols), [res_array])
                res_array = function(*values)

        res_dict = {"x": res_array}
        return res_dict

    def select_simulation_result(
        self, result: ca.DM | ca.MX, return_var_indexes: list[int]
    ) -> ca.DM | ca.MX:
        """Take result from self.simulate_sym() and return only a subset of results, defined by column
        index in return_var_names"""
        return result.get(False, return_var_indexes, 0)

    @_consistent_scaling_decorator
    def simulate(self, return_var_names: list[str] | None = None) -> ca.DM | ca.MX:
        """Wrapper for simulate_sym, that returns only results for variables, specified in return_var_names"""
        res = self.simulate_sym()["x"]
        if return_var_names is not None:
            return_var_index = []
            for var_name in return_var_names:
                return_var_index.append(self.mapping_algebraic_variables[var_name])
            res = self.select_simulation_result(res, return_var_index)
        return res

    def simulate_sym(self) -> dict[str, ca.MX | ca.DM]:
        self.call_arg["p"] = self._independent_variables

        res = self._call_simulator()
        return res

    def calculate_jac(self) -> dict[str, ca.MX | ca.DM]:
        self.call_arg["p"] = self._independent_variables

        res = self.jacobian.call(self.call_arg)
        return res
