import casadi as ca
import copy
from par_est import (
    VariableAlgebraic,
    VariableState,
    VariableParameter,
    VariableControl,
    VariableConstant,
    Variable,
    VariableList,
)


class Model(object):
    def __init__(self, variable_list):
        self.varlist_state = VariableList()
        self.varlist_algebraic = VariableList()
        self.varlist_independent = VariableList()
        self._varlist_constant = VariableList()
        self.varlist_all = VariableList()
        self.equations_differential = None
        self._equations_differential_constants = None
        self.equations_algebraic = None
        self._equations_differential_constants = None
        self.DAE = False
        self.values_constant = []

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    self.varlist_state.add_variable(VariableState(var.name))
                elif isinstance(var, VariableAlgebraic):
                    self.varlist_algebraic.add_variable(VariableAlgebraic(var.name))
                elif isinstance(var, VariableParameter) or isinstance(
                    var, VariableControl
                ):
                    self.varlist_independent.add_variable(type(var)(var.name))
                elif isinstance(var, VariableConstant):
                    self._varlist_constant.add_variable(VariableConstant(var.name))
                    self.values_constant.append(var.value)
                else:
                    raise VariableTypeError(var.name)
            else:
                raise VariableTypeError(var.name)

        self.varlist_all.update(self.varlist_state)
        self.varlist_all.update(self.varlist_algebraic)
        self.varlist_all.update(self.varlist_independent)
        self.varlist_all.update(self._varlist_constant)

    def add_equations_differential(self, equations):
        if self.equations_differential is None:
            self.equations_differential = ca.vcat(equations)
            self.substitute_constants()
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    def add_equations_algebraic(self, equations):
        if self.equations_algebraic is None:
            self.equations_algebraic = ca.vcat(equations)
            self.DAE = True
            self.substitute_constants()
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    def substitute_constants(self):
        self._equations_differential_constants = copy.deepcopy(
            self.equations_differential
        )
        self._equations_differential_constants = copy.deepcopy(
            self.equations_differential
        )
        self.varlist_all = VariableList()

        self.varlist_all.update(self.varlist_state)
        self.varlist_all.update(self.varlist_algebraic)
        self.varlist_all.update(self.varlist_independent)

        if self._varlist_constant:
            ode_system = {
                "x": self.varlist_state.get_casadi_var(),
                "p": ca.vertcat(self.varlist_independent.get_casadi_var()),
                "c": ca.vertcat(self._varlist_constant.get_casadi_var()),
                "ode": self.equations_differential,
            }

            if self.DAE:
                ode_system["alg"] = self.equations_algebraic
                ode_system["z"] = self.varlist_algebraic.get_casadi_var()

            if self.DAE:
                function_ode = ca.Function(
                    "eq_sys",
                    [
                        ode_system["x"],
                        ode_system["z"],
                        ode_system["p"],
                        ode_system["c"],
                    ],
                    [ode_system["ode"]],
                    ["x", "z", "p", "c"],
                    ["ode"],
                )

                function_alg = ca.Function(
                    "eq_sys",
                    [
                        ode_system["x"],
                        ode_system["z"],
                        ode_system["p"],
                        ode_system["c"],
                    ],
                    [ode_system["alg"]],
                    ["x", "z", "p", "c"],
                    ["alg"],
                )

                results_ode = function_ode(
                    ode_system["x"],
                    ode_system["z"],
                    ode_system["p"],
                    self.values_constant,
                )

                results_alg = function_alg(
                    ode_system["x"],
                    ode_system["z"],
                    ode_system["p"],
                    self.values_constant,
                )
                self.equations_differential = results_ode
                self.equations_algebraic = results_alg

            else:
                function = ca.Function(
                    "eq_sys",
                    [ode_system["x"], ode_system["p"], ode_system["c"]],
                    [ode_system["ode"]],
                    ["x", "p", "c"],
                    ["ode"],
                )
                results = function(
                    ode_system["x"], ode_system["p"], self.values_constant
                )

                self.equations_differential = results


class VariableTypeError(Exception):
    def __init__(self, name):
        message = f"Not a supported par_est_casadi variable class! Wrong variable with name: {name}"
        super().__init__(message)
