from collections import OrderedDict
from datetime import datetime, timedelta
from typing import List, Union

import casadi as ca
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

try:
    from opcua import ua
    from opcua.ua import NumericNodeId
    from optipal.client import OptiPALClient
except Exception:
    pass

ORIGIN_TS: pd.Timestamp = pd.Timestamp(year=1970, month=1, day=1)
""" Indicats a default zero timestamp for data, if date is irrelevant.
Chosen DateTime is the same, that is used by pd.to_datetime() by default.
"""


class Variable(object):
    def __init__(self, name, lb=None, ub=None):
        self.name: str = name
        self.casadi_var: ca.MX.sym = ca.MX.sym(self.name)
        # fixed is property in order to deal with VariableControlPiecewiseConstant properly
        self._fixed: Union[bool, List[bool]] = True
        self.opc_ua_id: Union[None, int] = None
        if not isinstance(self, VariableControlPiecewiseConstant):
            self.dataframe: pd.DataFrame = None
        self.guess = None
        self.lower_bound = lb
        self.upper_bound = ub
        self.variance = 1.0
        # attibute used to decide if variable should be plotted
        self.ignore_plotting = True

    @classmethod
    def get_subclasses(cls):
        for subclass in cls.__subclasses__():
            yield from subclass.get_subclasses()
            yield subclass

    def plot(self, ax=None):
        axis = self.dataframe.plot(ax=ax)
        plt.show()
        return axis

    def __repr__(self):
        return f"{self.name}\n{type(self)}\n{self.value}"

    def get_value_or_casadi(self) -> Union[float, ca.MX]:
        """Return either value at time=0 or casadi_variable.
        Used in Simulator for readability and less if statements.
        """
        if self.fixed:
            return self.value[0]
        else:
            return self.casadi_var

    def get_value_or_guess(self) -> float:
        """Return guess or value at time zero. Used further for
        readability"""
        if self.fixed:
            return self.value[0]
        else:
            return self.guess

    @property
    def value(self) -> List:
        """Returns a list with values of variables"""
        return self.dataframe[self.name].tolist()

    @property
    def time_absolute(self) -> pd.Series:
        """Returns a list which contains time_stamps with date and time"""
        return self.dataframe.index

    @property
    def time_relative(self) -> List:
        """Returns a list which contains timestamps in seconds.
        First time is considered to be zero second"""
        return (self.dataframe.index - self.dataframe.index[0]).total_seconds().tolist()

    def is_value_consistent(self) -> None:
        """Returns True if self.value is consistent or raise Error.

        Checks if index of self.value is increasing and unique,
        ensuring that first element in index is always time=0.
        Checks if any element in Data is Nan.

        Raises:
            BadVariableError: with descriptive text.
        """
        if isinstance(self.dataframe, pd.DataFrame):
            if not self.dataframe.index._is_strictly_monotonic_increasing:
                raise BadVariableError(self, "Value index is not unique or not sorted")
            if self.name not in self.dataframe.columns:
                raise BadVariableError(
                    self, "Column name in Variable.value dosn't equal Variable.name"
                )
            if not isinstance(
                self, (VariableAlgebraic, VariableControlPiecewiseConstant)
            ):
                if self.dataframe[self.name].hasnans:
                    raise BadVariableError(self, "Variable value has Nan")
        else:
            raise BadVariableError(self, "Value of Variable is of wrong type")

    @property
    def origin_ts(self) -> Union[None, pd.Timestamp]:
        """Propoerty that return the first Timestamp in self.value.

        Can be used to compare if Variables have same origin in .value.
        Does check self.value for consistensy.
        Returns:
            Union[None, pd.Timestamp]:
            None is self.value is None or Timestamp that corresponds to time=0
        """
        self.is_value_consistent()

        if self.value is None:
            return None
        else:
            if isinstance(self, VariableControlPiecewiseConstant):
                return self.time_absolute[0]
            else:
                return self.dataframe.index[0]

    @property
    def fixed(self):
        # VariableControlPiecewiseConstant is fixed only if all variables inside are fixed
        if isinstance(self, VariableControlPiecewiseConstant):
            fixed_list = []
            for variable in self.variable_list.values():
                fixed_list.append(variable.fixed)
            self._fixed = all(fixed_list)
        return self._fixed

    @fixed.setter
    def fixed(self, state):
        if isinstance(self, VariableControlPiecewiseConstant):
            for var in self.variable_list.values():
                var.fixed = state
        self._fixed = state

    @property
    def get_constraint_idas(self):
        """Constrain the solution y=[x,z].  0 (default): no constraint on yi,
        1: yi >= 0.0, -1: yi <= 0.0, 2: yi > 0.0, -2: yi < 0.0."}},"""

        if self.lower_bound == 0:
            constraint = 1
        elif self.lower_bound > 0:
            constraint = 2
        elif self.upper_bound == 0:
            constraint = -1
        elif self.upper_bound < 0:
            constraint = -2
        else:
            constraint = 0

        return constraint

    @property
    def lower_bound(self):
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, lower_bound):
        if lower_bound is None or lower_bound == -1e9:
            self._lower_bound = -ca.inf
        else:
            self._lower_bound = lower_bound

    @property
    def upper_bound(self):
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, upper_bound):
        if upper_bound is None or upper_bound == 1e9:
            self._upper_bound = ca.inf
        else:
            self._upper_bound = upper_bound

    def _dataframe_from_value(self, value: Union[None, float], origin=ORIGIN_TS):
        df = pd.DataFrame(
            [value],
            index=[origin],
            columns=[self.name],
            dtype="float64",
        )
        return df

    def set_dataframe_from_value_and_time(
        self, value: List[float], time_relative: List[float], origin="unix"
    ):
        if not len(value) == len(time_relative):
            raise ValueError(
                f"Value and time must have same length. Supplied Value:\n{value}\nTime:\n{time_relative}"
            )
        if not time_relative[0] == 0:
            raise ValueError("Time vector should start with 0, you supplied:\n{time}")

        time_series = pd.to_datetime(time_relative, unit="s", origin=origin)
        if isinstance(self, VariableControlPiecewiseConstant):
            raise NotImplementedError
        else:
            dataframe = pd.DataFrame(
                value, index=time_series, columns=[self.name], dtype="float64"
            )
            self.dataframe = dataframe


class VariableState(Variable):
    def __init__(
        self,
        name,
        starting_value: Union[float, None] = None,
        lb=None,
        ub=None,
        opc_ua_id=None,
    ):
        super().__init__(name, lb, ub)
        # Assuming that State Variables are always to be plotted
        self.ignore_plotting = False
        self.dataframe = self._dataframe_from_value(starting_value)
        self.opc_ua_id = opc_ua_id


class VariableAlgebraic(Variable):
    def __init__(self, name, guess=None, lb=None, ub=None, opc_ua_id=None):
        super().__init__(name, lb, ub)
        self.guess = guess
        self.opc_ua_id = opc_ua_id
        self.dataframe = self._dataframe_from_value(None)


class VariableParameter(Variable):
    def __init__(self, name, value=None, lb=None, ub=None):
        super().__init__(name, lb, ub)
        self.guess = value
        self.dataframe = self._dataframe_from_value(value)


class VariableControl(Variable):
    def __init__(self, name, value=None, lb=None, ub=None, opc_ua_id=None):
        super().__init__(name, lb, ub)
        if not isinstance(self, VariableControlPiecewiseConstant):
            self.dataframe = self._dataframe_from_value(value)
        self.guess = value
        self.opc_ua_id = opc_ua_id


class VariableControlPiecewiseConstant(VariableControl):
    """self.time - [time_stamps] list with time points of all variables in self.variables_list."""

    def __init__(self, name, value=None, lb=None, ub=None, opc_ua_id=None):
        super().__init__(name)
        self.variable_list = VariableList()
        var_t0 = VariableControl(name + "_t0", value, lb, ub, opc_ua_id)
        var_t0.fixed = True
        self.variable_list.add_variable(var_t0)

    @property
    def value(self):
        values = []
        for var in self.variable_list.values():
            values.extend(var.value)
        return values

    @property
    def time_absolute(self) -> pd.Series:
        time_list = []
        for variable in self.variable_list.values():
            time_list.append(variable.time_absolute[0])
        time_series = pd.Series(time_list)
        return time_series

    @property
    def time_relative(self) -> List[float]:
        time_series = self.time_absolute
        return (time_series - time_series.iloc[0]).dt.total_seconds().tolist()

    def to_dictionary(self):
        time_var_dict = dict(zip(self.time_relative, list(self.variable_list.values())))
        return time_var_dict

    def get_variable_at_time_absolute(self, time_stamp_absolute) -> VariableControl:
        index = pd.Index(self.time_absolute).get_loc(
            time_stamp_absolute, method="ffill"
        )
        return list(self.variable_list.values())[index]

    def get_variable_at_time_relative(self, time_stamp_relative) -> VariableControl:
        index = pd.Index(self.time_relative).get_loc(
            time_stamp_relative, method="ffill"
        )
        return list(self.variable_list.values())[index]

    def get_value_or_casadi(self, time_grid_relative) -> List:
        """This method is used to avoid following problem: if current Control is fixed at given time_stamp, simulator
        should use either - a fixed value, provided with Variable, or a value of a Control Variable from previous timestamp.
        Input:
                        t0      t1      t2      t3
        Value / Var     20      var_t1  var_t2  20
        Fixed / Unfixed Fixed   Unfixed Unfixed Fixed
        Result:
        Simulate with   20      var_t1  var_t2  var_t2
        """
        independent_variable = []
        last_unfixed_variable = None

        for time_stamp in time_grid_relative:
            var_at_timestamp = self.get_variable_at_time_relative(time_stamp)
            # This if statement is required for OED in order to use casadi_var from previous step, if it was already used. Without it, control variable will be fixed to some value for given timestep
            if var_at_timestamp.fixed:
                if last_unfixed_variable is None:
                    independent_variable.append(var_at_timestamp.get_value_or_casadi())
                else:
                    independent_variable.append(
                        last_unfixed_variable.get_value_or_casadi()
                    )
            else:
                last_unfixed_variable = var_at_timestamp
                independent_variable.append(last_unfixed_variable.get_value_or_casadi())

        return independent_variable

    def expand_horizon(self, times, values):
        if not len(times) == len(values):
            raise ValueError(
                "Length of times and values vector should be same. You supplied:\ntimes\n{times}\nvalues\n{values}"
            )
        if not len(self.time_relative) == 1:
            raise NotImplementedError(
                "Cannot be used to expand already expanded variable"
            )
        for index, (time, value) in enumerate(zip(times, values), 1):
            var = VariableControl(
                f"{self.name}_t{index}",
                value,
                self.variable_list.index(0).lower_bound,
                self.variable_list.index(0).upper_bound,
                self.opc_ua_id,
            )
            var.fixed = True
            var.dataframe = var._dataframe_from_value(
                value, self.time_absolute[0] + timedelta(seconds=time)
            )
            self.variable_list.add_variable(var)

    def set_horizon(self, times, values):
        """Used when control at time 0 should also be rewritten"""
        raise NotImplementedError

    @property
    # WIP, not tested
    def dataframe(self):
        values = []
        times = []
        for var in self.variable_list.values():
            values.append(var.value[0])
            times.append(var.time_absolute[0])

        dataframe = pd.DataFrame(
            values, index=times, columns=[self.name], dtype="float64"
        )

        return dataframe


class VariableConstant(Variable):
    def __init__(self, name, value=None, opc_ua_id=None):
        super().__init__(name)
        self.casadi_var = value
        self.dataframe = self._dataframe_from_value(value)
        self.opc_ua_id = opc_ua_id
        self.fixed = True


class VariableList(OrderedDict):
    def __init__(self):
        super().__init__()

    def __repr__(self):
        if bool(self):
            types = [type(item) for item in list(self.values())]
            counter_types = {x: types.count(x) for x in types}
            list_names = {var_type: [] for var_type in counter_types.keys()}
            message = f"Var list has {sum(counter_types.values())} variables:\n"
            for var in self.values():
                list_names[type(var)].extend([var.name])
            for var_type in counter_types.keys():
                if "VariableConstant" in str(var_type) or "VariableAlgebraic" in str(
                    var_type
                ):
                    print_list_names = str()
                else:
                    print_list_names = f":\n{list_names[var_type]}"
                message = (
                    message
                    + f"{var_type} of length {counter_types[var_type]}{print_list_names}\n"
                )
        else:
            message = f"Empty {type(self)}"
        return message

    def get_common_origin(
        self, strict=False, variable_type=Variable
    ) -> Union[pd.Timestamp, bool]:
        """Returns a common Timestamp of State, Algebraic, and Control variables. If no common origin exists - return ORIGIN_TS, strict is False

        Args:
            strict: if no common origin is found, return False instead of ORIGIN_TS
        """
        list_of_origins = []
        for variable in self.values():
            if isinstance(variable, variable_type):
                list_of_origins.append(variable.time_absolute[0])

        if len(set(list_of_origins)) < 2:
            return list_of_origins[0]
        else:
            if strict:
                return False
            else:
                return ORIGIN_TS

    @property
    def dataframe(self) -> pd.DataFrame:
        data_frame = pd.DataFrame()
        for var in self.values():
            data_frame = data_frame.join(var.dataframe, how="outer")
        return data_frame

    def index(self, var_index: int) -> Variable:
        """Return variable at given index (if VariableList was a List).

        Primary way to index Variables in VariableList is name of the variable.
        This method is used for debugging, and should not be used by inexperienced users.
        Args:
            var_index (int): var_index

        Returns:
            Variable: Variable that correspons to given index.
        """
        var: Variable = list(self.values())[var_index]
        return var

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

    def get_casadi_variables(self) -> ca.MX:
        """Returns a concatanated vector of all variables in a variable_list."""
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
                time_opcua: List[float] = []
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
        if fixation_list is None:
            for var in self.values():
                var.fixed = val
        else:
            for var in self.values():
                if var.name in fixation_list:
                    var.fixed = val

    def set_bounds(self, val=0.25, emerg_val=None):
        for var in self.values():
            if isinstance(var, VariableParameter) and var.fixed is False:
                value = var.value[0]
                if value > 0:
                    var.lower_bound = value * (1 - val)
                    var.upper_bound = value * (1 + val)
                elif value < 0:
                    var.lower_bound = value * (1 + val)
                    var.upper_bound = value * (1 - val)
                elif value == 0:
                    if emerg_val is None:
                        # Setting bounds for val == 0 without emerg_val is not implemented
                        raise (NotImplementedError)
                    else:
                        var.lower_bound = -emerg_val
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

    def plot(self, as_one_plot=False, algebraic=False, prefix=None, **kwargs):
        """Plots variables that are not ignored via var.ignore_plotting
        If as_one_plot is True, plot every variable on separate plot
        If algebraic is True, plot als algebraic variables

        prefix is used to append name to a variable name
        **kwargs are matplotlib options, for example marker='o'"""
        plot_varlist = self._get_varlist_to_plot(algebraic)

        if "subplots" not in kwargs:
            kwargs["subplots"] = True

        if as_one_plot is True:
            axes = []
            for var in plot_varlist.values():
                axes.append(var.plot())
            axes = np.array(axes)

        else:
            dataframe = plot_varlist.dataframe
            if prefix is not None:
                dataframe = dataframe.add_prefix(prefix)
            axes = dataframe.plot(**kwargs)

        plt.show()
        return axes

    def _get_varlist_to_plot(self, algebraic=False):
        """Return varlist that has only "plottable" variables"""
        plot_varlist = VariableList()
        for var in self.values():
            if not var.ignore_plotting:
                if isinstance(var, VariableState):
                    plot_varlist.add_variable(var)
                elif isinstance(var, VariableAlgebraic):
                    if algebraic is True:
                        plot_varlist.add_variable(var)
        return plot_varlist


class SameVariableNameError(Exception):
    def __init__(self, name):
        message = f"There is already an existing variable with the same name! Wrong variable with name: {name}"
        super().__init__(message)


class BadVariableError(Exception):
    def __init__(self, variable, message=None):
        if message is None:
            message = "Failed while using this variable:"
        message = message + f"\n{variable}"
        super().__init__(message)


class PlottingError(BadVariableError):
    pass
