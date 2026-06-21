from __future__ import annotations

import copy
from warnings import warn

import casadi as ca

from mopeds import (
    Model,
    VariableAlgebraic,
    VariableConstant,
    VariableControl,
    VariableList,
    VariableParameter,
    get_options,
    _consistent_scaling_decorator,
)


class SimulatorNLE:
    supported_solvers: list[str] = [
        "ipopt",
        "rootfinder",
        "newton",
        "fast_newton",
        "nlpsol",
    ]

    def __init__(
        self,
        model: Model,
        variable_list: VariableList,
        solver_settings: dict | None = None,
        solver_name: str = "nlpsol",
        *,
        use_bounds: bool = None,
    ):
        self._created_with_options = get_options()
        if use_bounds is not None:
            warn("use_bounds argument is ignored", FutureWarning, 2)

        if solver_name == "rootfinder":
            warn("solver_name=rootfinder is deprecated, use 'nlpsol'", FutureWarning, 2)
            solver_name = "nlpsol"

        self.model: Model = model
        if solver_name not in self.supported_solvers:
            raise TypeError(
                f"Provided integrator name {solver_name} is not supported. Only theese are: {self.supported_solvers}."
            )
        if solver_name == "ipopt":
            self._call_simulator = self._call_simulator_ipopt
        elif solver_name in ["newton", "nlpsol", "fast_newton"]:
            self._call_simulator = self._call_simulator_rootfinder

        self._solver_name: str = solver_name
        self._input_variable_list: VariableList = copy.deepcopy(variable_list)
        self.model.subsitute_casadi_symbols(self._input_variable_list)

        if solver_settings is not None:
            self.solver_settings: dict = solver_settings
        else:
            self.solver_settings = self.get_default_simulator_settings()

        self.__setup_variables()

        scaled_equations = ca.substitute(
            self.model.equations_algebraic,
            self._input_variable_list.get_casadi_variables(),
            self._input_variable_list.get_scaled_casadi_variables(),
        )
        self._model_equations = ca.cse(scaled_equations)

        self.function: ca.Function = ca.Function(
            "f",
            [
                self.model.varlist_algebraic(
                    self._input_variable_list
                ).get_casadi_variables(),
                self.model.varlist_independent(
                    self._input_variable_list
                ).get_casadi_variables(),
            ],
            [self._model_equations],
            ["x", "p"],
            ["rhs"],
        )
        if self._solver_name in ["newton", "nlpsol", "fast_newton"]:
            self.simulator: ca.Function = ca.rootfinder(
                "s", self._solver_name, self.function, self.solver_settings
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
                    "x": self.model.varlist_algebraic(
                        self._input_variable_list
                    ).get_casadi_variables(),
                    "p": self.model.varlist_independent(
                        self._input_variable_list
                    ).get_casadi_variables(),
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
        if self._solver_name == "nlpsol":
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
        elif self._solver_name in ["newton", "fast_newton"]:
            solver_settings = {
                "verbose": False,
                "print_in": False,
                "print_out": False,
                "expand": True,
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

    def __setup_variables(self) -> None:
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

        self.mapping_independent_variables: dict[str, int] = (
            mapping_independent_variables
        )
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
    def generate_exp_data(
        self, unfixed_variables: dict[str, float] = None
    ) -> VariableList:
        warn("generate_exp_data is deprecated. Use simulate()", FutureWarning)
        return self.simulate(unfixed_variables=unfixed_variables, return_varlist=True)[
            2
        ]

    def _call_simulator_rootfinder(self) -> ca.DM:
        """This method is needed to raise an error, if ipopt simulator fails to converge"""
        return self.simulator.call(self.call_arg)

    def _call_simulator_ipopt(self) -> ca.DM:
        """This method is needed to raise an error, if ipopt simulator fails to converge"""
        res = self.simulator.call(self.call_arg)

        if isinstance(res["x"], ca.DM):
            if not self.simulator.stats()["success"]:
                raise ValueError(
                    f"IPOPT failed as NLE solver:\n{self.simulator.stats()}"
                )
        return res

    def select_simulation_result(
        self, result: dict[str, ca.DM], return_var_names: list[str] | None
    ) -> ca.DM | ca.MX:
        """Take result from self.simulate_fast() and return only a subset of results, defined by
        list of variable names"""
        if return_var_names is not None:
            return_var_index = []
            for var_name in return_var_names:
                return_var_index.append(self.mapping_algebraic_variables[var_name])
            res_selected = result["x"].get(False, return_var_index, 0)
        else:
            res_selected = {"x": None}

        return res_selected

    @_consistent_scaling_decorator
    def simulate(
        self,
        *,
        return_var_names: list[str] | None = None,
        unfixed_variables: dict[str, float] | None = None,
        return_varlist: bool = True,
    ) -> (dict(str, ca.MX), ca.MX | None, VariableList | None):
        """Wrapper for simulate_fast, that returns scaled results.
        res is returned as dict to be consistent with Dynamic simulations."""
        res = self.simulate_fast()
        res_selected = self.select_simulation_result(res, return_var_names)

        if not isinstance(res["x"], ca.DM):
            if unfixed_variables is None:
                raise ValueError("You need to supply values for unfixed variables")
            else:
                unfixed_symbols = ca.symvar(res["x"])
                values = []
                for symbol in unfixed_symbols:
                    var_name = symbol.name()
                    values.append(
                        self._input_variable_list[var_name].scale_from_original(
                            unfixed_variables[var_name]
                        )
                    )

                function = ca.Function("f", unfixed_symbols, [res["x"]])
                res = {"x": function(*values)}
        else:
            res = {"x": res["x"]}

        if return_varlist:
            variables = VariableList()
            for var_name, res_var in zip(
                self.model.varlist_algebraic(self._input_variable_list).keys(),
                res["x"].toarray(),
            ):
                if return_var_names is not None:
                    if var_name not in return_var_names:
                        continue
                var = self._input_variable_list[var_name]
                new_var = copy.deepcopy(var)
                new_var.casadi_var = None

                value = var.scale_to_original(res_var)
                new_var.value = float(value[0])
                variables.add_variable(new_var)
        else:
            variables = None

        return res, res_selected, variables

    @_consistent_scaling_decorator
    def simulate_sym_unfixed(
        self, unfixed_variables: dict[str, float] = None
    ) -> dict[str, ca.DM]:
        warn(
            "simulate_sym_unfixed is deprecated. Use simulate(unfixed_variables=unfixed_variables, return_varlist=False)",
            FutureWarning,
        )
        return self.simulate(unfixed_variables=unfixed_variables, return_varlist=False)[
            0
        ]

    def simulate_sym(self) -> dict[str, ca.MX | ca.DM]:
        warn("simulate_sym is deprecated. Use simulate_fast", FutureWarning, 2)
        return self.simulate_fast()

    def simulate_fast(self) -> dict[str, ca.MX | ca.DM]:
        """Should be used internally for fast calculations. Output of the solution is not scaled"""
        self.call_arg["p"] = self._independent_variables

        res = self._call_simulator()
        return res

    def calculate_jac(self) -> dict[str, ca.MX | ca.DM]:
        self.call_arg["p"] = self._independent_variables

        res = self.jacobian.call(self.call_arg)
        return res
