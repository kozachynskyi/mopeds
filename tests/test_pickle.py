import par_est
import par_est.examples
import numpy as np
import pickle
import copy
import pathlib


def test_pickling_objects(tmp_path):
    variable_list, model = par_est.examples.cstr_dae()
    variable = variable_list["e0_T"]
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    simulation = par_est.Simulator(model, time_grid, var_list_fixed)
    var_list_exp = simulation.generate_exp_data()

    for key, var in var_list_exp.items():
        variable_list[key] = var
    for var in variable_list.values():
        var.fixed = True
    variable_list["e0_E_r1"].fixed = False

    pe = par_est.ParameterEstimation(model, [variable_list])
    oed = par_est.OptimalExperimentalDesign(model, [variable_list], time_grid)
    for object_current in [variable, variable_list, model, simulation, pe, oed]:
        file_write = open(tmp_path / "tmp.pkl", "wb")
        pickler = par_est.MXPickler(file_write)
        pickler.dump(object_current)
        file_write.close()

        # Open file to read data back
        file_read = open(tmp_path / "tmp.pkl", "rb")
        loaded_object = pickle.load(file_read)
        assert type(loaded_object) == type(object_current)


if __name__ == "__main__":
    test_pickling_objects(pathlib.Path.cwd())
