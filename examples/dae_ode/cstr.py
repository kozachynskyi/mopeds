import copy
from datetime import datetime, timedelta  # noqa: F401

import numpy as np

import mopeds
import mopeds.examples

if __name__ == "__main__":

    variable_list, m = mopeds.examples.cstr_ode()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 40)
    time_grid = np.insert(time_grid, 0, 0)

    # Generate experimental data
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    sim_fixed = mopeds.Simulator(m, time_grid, var_list_fixed)
    res = sim_fixed.simulate()
    var_list_exp = sim_fixed.simulate()[2]

    # start_time = datetime(2018, 1, 1, 1, 0, 0, 0) + timedelta(days=1)
    # end_time = start_time + timedelta(seconds=1000)
    # var_list3 = copy.deepcopy(var_list_fixed)
    # var_list_exp.write_data_opcua(start_time)
    # var_list3.get_data_opcua(start_time, end_time)
