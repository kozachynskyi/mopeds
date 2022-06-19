import par_est
import par_est.examples
import numpy as np
import pickle
import copy
import pathlib
import pytest
import casadi as ca


@pytest.mark.parametrize("piecewise", [True, False])
def test_pickling_objects(tmp_path, piecewise):
    variable_list, model = par_est.examples.cstr_dae(piecewise)
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


@pytest.mark.parametrize("piecewise", [True, False])
def test_varlist_simulation_reusability(tmp_path, piecewise):
    """ Test if simulator created from varlist and pickled/unpickled varlist provides same results."""
    variable_list, model = par_est.examples.pendulum_dae_1(piecewise)
    time_grid = np.linspace(0, 1, 3)
    variable_list["g"].value = 12.0
    simulation = par_est.Simulator(model, time_grid, variable_list)
    res_before_pickle = simulation.simulate()

    file_write = open(tmp_path / "tmp.pkl", "wb")
    pickler = par_est.MXPickler(file_write)
    pickler.dump(variable_list)
    file_write.close()
    file_read = open(tmp_path / "tmp.pkl", "rb")
    variable_list_after = pickle.load(file_read)
    simulation = par_est.Simulator(model, time_grid, variable_list_after)
    res_after_pickle = simulation.simulate()

    assert np.isclose(
        ca.vertcat(res_before_pickle["xf"], res_before_pickle["zf"]), ca.vertcat(res_after_pickle["xf"], res_after_pickle["zf"])
    ).all()

    variable_list, model = par_est.examples.pendulum_dae_1(piecewise, variable_list)
    simulation = par_est.Simulator(model, time_grid, variable_list)
    res_after_pickle = simulation.simulate()

    assert np.isclose(
        ca.vertcat(res_before_pickle["xf"], res_before_pickle["zf"]), ca.vertcat(res_after_pickle["xf"], res_after_pickle["zf"])
    ).all()


if __name__ == "__main__":
    test_pickling_objects(pathlib.Path.cwd())
    # test_varlist_reusability(pathlib.Path.cwd(), True)
