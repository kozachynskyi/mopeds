from __future__ import annotations

import casadi as ca
from dataclasses import dataclass
from functools import cached_property

from mopeds import VariableList, VariableConstant

@dataclass
class DummyClass():
    casadi_var: ca.MX


class Model(object):
    """Model class is used to get lists of variables used in this model,
    create equations and determine if model is DAE or ODE.
    """

    def __init__(self, variable_list: VariableList, name: str = "default") -> None:
        self.equations_differential: ca.MX = None
        self.equations_algebraic: ca.MX = None
        self.DAE: bool = False

        casadi_vars = []
        for var in variable_list.values():
            casadi_vars.append(var.casadi_var)

        # This dictionary is used to to consistently iterate over all variables in Simulation and Optimization that use the same model
        self.variables_all = dict(zip(variable_list.keys(), casadi_vars))

        self.name = name

    def add_equations_differential(self, equations: list[ca.MX]) -> None:
        if self.equations_differential is None:
            self.equations_differential = ca.vcat(equations)
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    def add_equations_algebraic(self, equations: list[ca.MX]) -> None:
        if self.equations_algebraic is None:
            self.equations_algebraic = ca.vcat(equations)
            self.DAE = True
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    @cached_property
    def varlist_all(self):
        return_dict = {}
        for name, casadi_var in self.variables_all.items():
            return_dict[name] = DummyClass(casadi_var)
        return return_dict

    def subsitute_casadi_symbols(self, variable_list: VariableList) -> VariableList:
        """Replace casadi symbols of the valist with the ones from model"""
        for var in variable_list.values():
            if not isinstance(var, VariableConstant):
                var.casadi_var = self.variables_all[var.name]
        return variable_list

    def varlist(self, variable_list: VariableList) -> VariableList:
        """Returns model-ordered varible list with variables included in model"""
        return variable_list._get_sorted_varlist(list(self.variables_all.keys()), raise_error=False)

    def varlist_state(self, variable_list: VariableList) -> VariableList:
        return self.varlist(variable_list).get_state()

    def varlist_algebraic(self, variable_list: VariableList) -> VariableList:
        return self.varlist(variable_list).get_algebraic()

    def varlist_independent(self, variable_list: VariableList) -> VariableList:
        return self.varlist(variable_list).get_independent()
