import mopeds
import mopeds.examples
import numpy as np
import pickle
import copy
import pathlib
import pytest
import casadi as ca


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
def test_pickling_objects(tmp_path, piecewise, dae, use_constant):
    variable_list, model = mopeds.examples.cstr(piecewise, dae, use_constant)
    variable = variable_list["e0_T"]
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    simulation = mopeds.Simulator(model, time_grid, var_list_fixed)
    var_list_exp = simulation.generate_exp_data()

    for key, var in var_list_exp.items():
        variable_list[key] = var
    for var in variable_list.values():
        var.fixed = True
    variable_list["e0_E_r1"].fixed = False

    pe = mopeds.ParameterEstimation(model, [variable_list])
    oed = mopeds.OptimalExperimentalDesign(model, [variable_list], time_grid)
    for object_current in [variable, variable_list, model, simulation, pe, oed]:
        file_write = open(tmp_path / "tmp.pkl", "wb")
        pickler = mopeds.MXPickler(file_write)
        pickler.dump(object_current)
        file_write.close()

        # Open file to read data back
        file_read = open(tmp_path / "tmp.pkl", "rb")
        loaded_object = pickle.load(file_read)
        assert type(loaded_object) == type(object_current)


@pytest.mark.parametrize("piecewise", [True, False])
def test_varlist_simulation_reusability(tmp_path, piecewise):
    """Test if simulator created from varlist and pickled/unpickled varlist provides same results."""
    variable_list, model = mopeds.examples.pendulum_dae_1(piecewise)
    time_grid = np.linspace(0, 1, 3)
    variable_list["g"].value = 12.0
    simulation = mopeds.Simulator(model, time_grid, variable_list)
    res_before_pickle = simulation.simulate_sym()

    file_write = open(tmp_path / "tmp.pkl", "wb")
    pickler = mopeds.MXPickler(file_write)
    pickler.dump(variable_list)
    file_write.close()
    file_read = open(tmp_path / "tmp.pkl", "rb")
    variable_list_after = pickle.load(file_read)
    simulation = mopeds.Simulator(model, time_grid, variable_list_after)
    res_after_pickle = simulation.simulate_sym()

    assert np.isclose(
        ca.vertcat(res_before_pickle["xf"], res_before_pickle["zf"]),
        ca.vertcat(res_after_pickle["xf"], res_after_pickle["zf"]),
    ).all()

    variable_list, model = mopeds.examples.pendulum_dae_1(piecewise, variable_list)
    simulation = mopeds.Simulator(model, time_grid, variable_list)
    res_after_pickle = simulation.simulate_sym()

    assert np.isclose(
        ca.vertcat(res_before_pickle["xf"], res_before_pickle["zf"]),
        ca.vertcat(res_after_pickle["xf"], res_after_pickle["zf"]),
    ).all()


if __name__ == "__main__":
    test_pickling_objects(pathlib.Path.cwd())
    # test_varlist_reusability(pathlib.Path.cwd(), True)
