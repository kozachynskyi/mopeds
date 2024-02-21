import copy
import pickle

import numpy as np

import mopeds
import mopeds.examples

if __name__ == "__main__":
    """This example shows how one can store Data in pickle files and
    check them later. Proposed method makes all objects that depend on
    casadi Types unusable, because casadi Types cannot be pickled.

    It just means that you can still plot VariableList, but cannot use it
    for further simulations, because corresponding casadi.MX variables
    are gone.
    """

    variable_list, m = mopeds.examples.cstr_dae()
    time_grid = np.linspace(10, 10000, 40)
    time_grid = np.insert(time_grid, 0, 0)

    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    sim_fixed = mopeds.Simulator(m, time_grid, var_list_fixed)
    res_simple = sim_fixed.simulate_fast()
    res = sim_fixed.simulate(algebraic=True)[2]
    variable = variable_list["e0_T"]

    objects_to_pickle_names = [
        "variable",
        "variable_list",
        "sim_fixed",
        "res_simple",
        "res",
    ]

    objects_to_store_dict = dict(
        zip(objects_to_pickle_names, [eval(x) for x in objects_to_pickle_names])
    )

    # Create pickle file were data is stored and "dump" dict
    file_write = open("tmp.pkl", "wb")
    pickler = mopeds.MXPickler(file_write)
    pickler.dump(objects_to_store_dict)
    file_write.close()

    # Open file to read data back
    file_read = open("tmp.pkl", "rb")
    loaded_dict = pickle.load(file_read)
    loaded_dict["res"].plot()
    file_read.close()

    # Notice that all Casadi Objects were transformed to Strings and are not usable anymore
    print(f"Variable.casadi_var type before pickling {type(variable.casadi_var)}")
    print(
        f"Variable.casadi_var type after pickling {type(loaded_dict['variable'].casadi_var)}"
    )
