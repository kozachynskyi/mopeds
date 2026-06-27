import mopeds
import numpy as np
import mopeds.examples
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
def test_dae_initials_calculation(piecewise):
    varlist, model = mopeds.examples.empy_dae(piecewise)
    varlist["X1"].value = 1
    varlist["C"].fixed = True
    varlist["P"].fixed = True
    time_grid = np.array([1, 2])

    sim = mopeds.Simulator(model, time_grid, varlist)
    assert sim._initial_algebraic[0] == 0
    sim.calculate_algebraic_initials(apply_intials=True)
    assert sim._initial_algebraic[0] == -1


@pytest.mark.parametrize("piecewise", [True, False])
def test_pendulum_dae(piecewise):
    varlist, model = mopeds.examples.pendulum_dae_1(piecewise)

    time_grid = np.linspace(0, 1, 3)
    for var in varlist.values():
        var.fixed = True
    with mopeds.options(variable_scaling=False):
        sim = mopeds.Simulator(model, time_grid, varlist)
        res_tau = sim.simulate(algebraic=True)[2].dataframe
        res_tau.drop(columns="L", inplace=True, errors="ignore")
        res = np.array(
            [
                [3.42289, 4.68624],
                [1.96674, 2.34688],
                [3.6447, 1.74332],
                [-1.84705, -6.30866],
                [1.16669, -1.11495],
            ]
        )
        assert np.isclose(res_tau.iloc[1:].T, res).all()

        res_tau = sim.simulate(return_var_names=["u", "v"])[2].dataframe
        res_tau.drop(columns="L", inplace=True, errors="ignore")
        assert np.isclose(res_tau.iloc[1:].T, res[[1, 3]]).all()
        res_tau = sim.simulate(return_var_names=["u", "v"], return_varlist=False)[2]
        assert res_tau is None


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_cstr(piecewise, dae, use_constant, scaling):
    with mopeds.options(variable_scaling=scaling):
        variable_list, m = mopeds.examples.cstr(piecewise)
        # Create time-grid. Zero should be first
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)

        for var in variable_list.values():
            if isinstance(
                var,
                (
                    mopeds.VariableControl,
                    mopeds.VariableParameter,
                    mopeds.VariableControlPiecewiseConstant,
                ),
            ):
                var.fixed = False

        for i in range(5):
            for j in range(2):
                if i == 1:
                    variable_list["e0_U"].fixed = True
                elif i == 2:
                    variable_list["e0_T_in"].fixed = True
                elif i == 3:
                    variable_list["e0_c_i2"].fixed = True
                elif i == 4:
                    for var in variable_list.values():
                        var.fixed = True

                sim = mopeds.Simulator(m, time_grid, variable_list, simulate_jac=True)
                if j == 0:
                    res_simple = sim.simulate_fast()
                else:
                    res_simple = sim.simulate_jac()

                if i == 4:
                    res = sim.simulate()[2]
                    if piecewise:
                        assert len(res) == 7
                    else:
                        assert len(res) == 5
                    if sim.model.DAE:
                        res = sim.simulate(algebraic=True)[2]
                        if piecewise:
                            assert len(res) == 8
                        else:
                            assert len(res) == 6
                else:
                    with pytest.raises(ValueError):
                        res = sim.simulate()[2]

                if j == 1:
                    assert res_simple["jac_xf_p"].size() == (5, 76)
                assert res_simple["xf"].size() == (5, 4)

                if i == 0:
                    assert sim._independent_variables[0][6].is_symbolic()
                elif i == 1:
                    assert not sim._independent_variables[0][6].is_symbolic()
                    assert sim._independent_variables[0][15].is_symbolic()
                elif i == 2:
                    assert not sim._independent_variables[0][15].is_symbolic()
                elif i == 3:
                    assert sim._independent_variables[0][10].is_symbolic()
                elif i == 4:
                    assert not sim._independent_variables[0][15].is_symbolic()


@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_piecewise(dae, use_constant, scaling):
    with mopeds.options(variable_scaling=scaling):
        time_grid = np.array([0, 10, 2000, 4000, 10000])
        time_grid_piecewise = np.array([0, 10, 10000])

        variable_list, m = mopeds.examples.cstr(False, dae, use_constant)
        for var in variable_list.values():
            var.fixed = True
        sim = mopeds.Simulator(m, time_grid, variable_list, simulate_jac=True)

        res = sim.simulate_jac()

        variable_list, m = mopeds.examples.cstr(True, dae, use_constant)
        T_in = variable_list["e0_T_in"]
        T_in.expand_horizon([2000, 4000], [373, 373])
        for var in variable_list.values():
            var.fixed = True
        sim = mopeds.Simulator(m, time_grid_piecewise, variable_list, simulate_jac=True)

        res_piecewise = sim.simulate_jac()

        print(res["xf"])
        assert np.isclose(res["xf"], res_piecewise["xf"]).all()
        assert np.isclose(res["jac_xf_p"], res_piecewise["jac_xf_p"]).all()
        if m.DAE:
            assert np.isclose(res["zf"], res_piecewise["zf"]).all()


@pytest.mark.parametrize("piecewise", [True, False])
@pytest.mark.parametrize("dae", [True, False])
@pytest.mark.parametrize("use_constant", [True, False])
@pytest.mark.parametrize("scaling", [True, False])
def test_steadystate(piecewise, dae, use_constant, scaling):
    with mopeds.options(variable_scaling=scaling):
        variable_list, m = mopeds.examples.cstr(piecewise, dae, use_constant)
        # Create time-grid. Zero should be first
        time_grid = np.linspace(10, 100000, 4)
        time_grid = np.insert(time_grid, 0, 0)

        sim = mopeds.Simulator(m, time_grid, variable_list)

        sim_res = sim.simulate_fast()

        steady_state = sim.calculate_steady_state()

        assert np.isclose(sim_res["xf"][:, -1], steady_state[0:5]).all()

        if m.DAE:
            assert np.isclose(sim_res["zf"][:, -1], steady_state[5]).all()


if __name__ == "__main__":
    pass
    # test_pendulum_dae(True)
    # test_cstr(True, True, True, False)
    # test_steadystate(True, True, True, True)
    # test_dae_initials_calculation(True)
    test_piecewise(False, True, True)
    # test_constraints_idas()
