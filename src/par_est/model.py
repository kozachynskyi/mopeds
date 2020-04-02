# import copy
# from collections import OrderedDict
# from datetime import datetime, timedelta

import casadi as ca
# import matplotlib.cm as cm
# import numpy as np
# from matplotlib import pyplot as plt
# from opcua import ua
# from opcua.ua import NumericNodeId
# from optipal.client import OptiPALClient

from par_est import VariableParameter, VariableAlgebraic, VariableState, Variable, VariableList


class Model(object):

    """Docstring for model. """

    def __init__(self, variable_list):
        """TODO: to be defined. """
        self.varlist_state = VariableList()
        self.varlist_algebraic = VariableList()
        self.variables = VariableList()
        self._all_variables = VariableList()
        self.differential_equations = None
        self.algebraic_equations = None

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    self.varlist_state.add_variable(VariableState(var.name))
                elif isinstance(var, VariableAlgebraic):
                    self.varlist_algebraic.add_variable(VariableAlgebraic(var.name))
                else:
                    self.variables.add_variable(VariableParameter(var.name))
            else:
                raise (ValueError)

        self._all_variables.update(self.varlist_state)
        self._all_variables.update(self.varlist_algebraic)
        self._all_variables.update(self.variables)

    def add_differential_equations(self, equations):
        if self.differential_equations is None:
            self.differential_equations = ca.vcat(equations)
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    def add_algebraic_equations(self, equations):
        if self.algebraic_equations is None:
            self.algebraic_equations = ca.vcat(equations)
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)
