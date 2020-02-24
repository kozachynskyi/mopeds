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

from par_est import Parameter_variable, State_variable, Variable, VariableList


class Model(object):

    """Docstring for model. """

    def __init__(self, variable_list):
        """TODO: to be defined. """
        self.states = VariableList()
        self.variables = VariableList()
        self._all_variables = VariableList()
        self.equations = None

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, State_variable):
                    self.states.add_variable(State_variable(var.name))
                else:
                    self.variables.add_variable(Parameter_variable(var.name))
            else:
                raise (ValueError)

        self._all_variables.update(self.states)
        self._all_variables.update(self.variables)

    def add_equations(self, equations):
        if self.equations is None:
            self.equations = ca.vcat(equations)
        else:
            raise (NotImplementedError)
