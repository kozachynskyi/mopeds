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

from par_est import Parameter_variable, Algebraic_variable, State_variable, Variable, VariableList


class Model(object):

    """Docstring for model. """

    def __init__(self, variable_list):
        """TODO: to be defined. """
        self.states = VariableList()
        self.algebraic_variables = VariableList()
        self.variables = VariableList()
        self._all_variables = VariableList()
        self.differential_equations = None
        self.algebraic_equations = None

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, State_variable):
                    self.states.add_variable(State_variable(var.name))
                elif isinstance(var, Algebraic_variable):
                    self.algebraic_variables.add_variable(Algebraic_variable(var.name))
                else:
                    self.variables.add_variable(Parameter_variable(var.name))
            else:
                raise (ValueError)

        self._all_variables.update(self.states)
        self._all_variables.update(self.algebraic_variables)
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
