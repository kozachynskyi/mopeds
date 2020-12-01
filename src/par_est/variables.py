from collections import OrderedDict
from datetime import datetime, timedelta
import copy

import casadi as ca
import numpy as np
from matplotlib import pyplot as plt

try:
    from opcua import ua
    from opcua.ua import NumericNodeId
    from optipal.client import OptiPALClient
except Exception:
    pass


class Variable(object):
    def __init__(self, name):
        self.name = name
        self.casadi_var = ca.MX.sym(self.name)
        self.fixed = False
        self.opc_ua_id = None
        self.starting_value = None
        self.value = None
        self.guess = None
        self.lower_bound = None
        self.upper_bound = None

    @classmethod
    def get_subclasses(cls):
        for subclass in cls.__subclasses__():
            yield from subclass.get_subclasses()
            yield subclass


class VariableState(Variable):
    def __init__(self, name, starting_value=None, opc_ua_id=None):
        super().__init__(name)
        self.starting_value = starting_value
        self.value = ExperimentData()
        self.opc_ua_id = opc_ua_id


class VariableAlgebraic(Variable):
    def __init__(self, name, guess=None, opc_ua_id=None):
        super().__init__(name)
        self.guess = guess
        self.opc_ua_id = opc_ua_id
        self.value = ExperimentData()


class VariableParameter(Variable):
    def __init__(self, name, value=None, lb=None, ub=None):
        super().__init__(name)
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub


class VariableControl(Variable):
    def __init__(self, name, value=None, lb=None, ub=None, opc_ua_id=None):
        super().__init__(name)
        self.value = value
        self.lower_bound = lb
        self.upper_bound = ub
        self.opc_ua_id = opc_ua_id


class VariableList(OrderedDict):
    def __init__(self):
        super().__init__()

    def add_variable(self, variable: Variable):
        if variable.name in self:
            raise SameVariableNameError(variable.name)
        else:
            self.update({variable.name: variable})

    def get_variable_name(self):
        names = []
        for var in self.values():
            names.append(var.name)
        return names

    def get_casadi_var(self):
        casadi_vars = []
        for var in self.values():
            casadi_vars.append(var.casadi_var)
        return ca.vcat(casadi_vars)

    def get_data_opcua(self, time_start: datetime, time_stop: datetime):
        client = OptiPALClient("opc.tcp://admin@localhost:4840")  # type: OptiPALClient
        client.connect()
        try:
            ns_working = client.get_working_ns_idx()
            for var in self.values():
                values_opcua = []
                time_opcua = []
                if isinstance(var, VariableState):
                    sensor = client.get_node(NumericNodeId(var.opc_ua_id, ns_working))
                    process_value = client.get_child_simple(sensor, ["d:ProcessValue"])
                    results = process_value.read_raw_history(
                        time_start, time_stop, 1000
                    )
                    var.value = ExperimentData()

                    for result in results:
                        if not time_opcua:
                            time_opcua.append(0.0)
                            time_zero = result.SourceTimestamp
                            var.starting_value = result.Value.Value
                        else:
                            time_from_ref = (
                                result.SourceTimestamp - time_zero
                            ).total_seconds()
                            time_opcua.append(time_from_ref)

                        values_opcua.append(result.Value.Value)

                    var.value.value = np.array(values_opcua)
                    var.value.time = np.array(time_opcua)
        finally:
            client.disconnect()

    def set_variable_list_fixed(self, fix_list=None):
        self._list_fixation(fix_list, True)

    def set_variable_list_unfixed(self, unfix_list=None):
        self._list_fixation(unfix_list, False)

    def _list_fixation(self, fixation_list, val):
        if fixation_list == None:
            for var in self.values():
                var.fixed = val
        else:
            for var in self.values():
                if var.name in fixation_list:
                    var.fixed = val

    def set_starting_values(self, values_simulation, NLE_Flag=False):
        for key, var in values_simulation.items():
            if NLE_Flag == True:
                var.guess = var.value.value
                self[key] = var
            elif NLE_Flag == False:
                var.starting_value = var.value.value
                self[key] = var

    def set_bounds(self, val=0.25, emerg_val=None):
        for var in self.values():
            if isinstance(var, VariableParameter) and var.fixed == False:
                if var.value > 0:
                    var.lower_bound = var.value - var.value * val
                    var.upper_bound = var.value + var.value * val
                elif var.value < 0:
                    var.lower_bound = var.value + var.value * val
                    var.upper_bound = var.value - var.value * val
                elif var.value == 0:
                    if emerg_val is None:
                        # Setting bounds for val == 0 without emerg_val is not implemented
                        raise (NotImplementedError)
                    else:
                        var.lower_bound = - emerg_val
                        var.upper_bound = emerg_val
                else:
                    # Setting bounds for arrays is not implemented
                    raise (NotImplementedError)
                var.guess = var.lower_bound

    def write_data_opcua(self, time_start: datetime):
        client = OptiPALClient("opc.tcp://admin@localhost:4840")  # type: OptiPALClient
        client.connect()
        try:
            time_zero = time_start
            ns_working = client.get_working_ns_idx()
            for var in self.values():
                if isinstance(var, VariableState):
                    sensor = client.get_node(NumericNodeId(var.opc_ua_id, ns_working))
                    process_value = client.get_child_simple(sensor, ["d:ProcessValue"])
                    for value, time in zip(var.value.value, var.value.time):
                        datavalue = ua.DataValue(value)
                        datavalue.SourceTimestamp = time_zero + timedelta(seconds=time)
                        process_value.set_attribute(ua.AttributeIds.Value, datavalue)
        finally:
            client.disconnect()

    def plot_states(self, as_one_plot=False, algebraic=False):
        # Choose only state variables
        plot_varlist = VariableList()
        for var in self.values():
            if var.value.value is None:
                raise PlottingError("variables")
            elif var.value.time is None:
                raise PlottingError("time grid")

            if isinstance(var, VariableState):
                plot_varlist.add_variable(var)
            elif isinstance(var, VariableAlgebraic):
                if algebraic is True:
                    plot_varlist.add_variable(var)

        if as_one_plot is True:
            for var in plot_varlist.values():
                plt.plot(var.value.time, var.value.value, label=var.name)
            plt.legend()
        else:
            figure, axes_array = plt.subplots(len(plot_varlist))
            for var, ax in zip(plot_varlist.values(), axes_array):
                ax.plot(var.value.time, var.value.value, label=var.name)
                ax.legend()
        plt.show()


class ExperimentData(object):
    def __init__(self):
        self.time = None
        self.value = None


class SameVariableNameError(Exception):
    def __init__(self, name):
        message = f"There is already an existing variable with the same name! Wrong variable with name: {name}"
        super().__init__(message)


class PlottingError(Exception):
    def __init__(self, error_switch):
        if error_switch == "variables":
            message = "There are no variables to plot!"
        elif error_switch == "time grid":
            message = "There is no time grid to plot against!"
        else:
            message = "Plotting not possible"
        super().__init__(message)
