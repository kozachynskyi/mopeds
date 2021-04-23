import par_est
import numpy as np
import logging
import par_est.examples
import casadi as ca
import copy
import pytest


@pytest.mark.parametrize("piecewise", [True, False])
def test_dae_initials_calculation(piecewise):
    varlist, model = par_est.examples.empy_dae(piecewise)
    varlist["X1"].value.value = [1]
    varlist["C"].fixed = True
    varlist["P"].fixed = True
    time_grid = [1, 2]

    sim = par_est.Simulator(model, time_grid, varlist)
    assert sim._initial_algebraic[0] == 0
    sim.calculate_algebraic_initials(apply_intials=True)
    assert sim._initial_algebraic[0] == -1


@pytest.mark.parametrize("piecewise", [True, False])
def test_pendulum_dae(piecewise):
    varlist, model = par_est.examples.pendulum_dae_1(piecewise)

    time_grid = np.linspace(0, 1, 3)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res_tau = sim.simulate()
    res = sim.integrator(
        x0=sim._initial_state,
        z0=sim._initial_algebraic,
        p=ca.vertcat(sim._variables[0] * sim.scaling),
    )

    assert np.isclose(
        ca.vertcat(res_tau["xf"], res_tau["zf"]), ca.vertcat(res["xf"], res["zf"])
    ).all()


@pytest.mark.parametrize("piecewise", [True, False])
def test_cstr(piecewise):
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
        par_est.examples.cstr_ode_constant,
    ]:
        variable_list, m = cstr_model(piecewise)
        # Create time-grid. Zero should be first
        time_grid = np.linspace(10, 10000, 4)
        time_grid = np.insert(time_grid, 0, 0)

        for var in variable_list.values():
            if isinstance(var, (par_est.VariableControl, par_est.VariableParameter, par_est.VariableControlPiecewiseConstant)):
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

                sim = par_est.Simulator(m, time_grid, variable_list)
                if j == 0:
                    res_simple = sim.simulate()
                else:
                    res_simple = sim.simulate_jac()

                res = sim.generate_exp_data()
                if j == 1:
                    assert res_simple["jac_xf_p"].size() == (5, 76)
                assert res_simple["xf"].size() == (5, 4)
                assert len(res) == 5

                if i == 0:
                    assert sim._variables[0][6].is_symbolic()
                elif i == 1:
                    assert not sim._variables[0][6].is_symbolic()
                    assert sim._variables[0][15].is_symbolic()
                elif i == 2:
                    assert not sim._variables[0][15].is_symbolic()
                elif i == 3:
                    assert sim._variables[0][10].is_symbolic()
                elif i == 4:
                    assert not sim._variables[0][15].is_symbolic()


def test_piecewise():
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
        par_est.examples.cstr_ode_constant,
    ]:
        time_grid = [0, 10, 2000, 4000, 10000]
        time_grid_piecewise = [0, 10, 10000]

        variable_list, m = cstr_model(False)
        for var in variable_list.values():
            var.fixed = True
        sim = par_est.Simulator(m, time_grid, variable_list)

        res = sim.simulate_jac()

        variable_list, m = cstr_model(True)
        T_in = variable_list["e0_T_in"]
        T_in.expand_horizon([2000, 4000], [373, 373])
        for var in variable_list.values():
            var.fixed = True
        sim = par_est.Simulator(m, time_grid_piecewise, variable_list)

        res_piecewise = sim.simulate_jac()

        assert np.isclose(res["xf"], res_piecewise["xf"]).all()
        assert np.isclose(res["jac_xf_p"], res_piecewise["jac_xf_p"]).all()
        if m.DAE:
            assert np.isclose(res["zf"], res_piecewise["zf"]).all()

@pytest.mark.parametrize("piecewise", [True, False])
def test_steadystate(piecewise):
    for cstr_model in [
        par_est.examples.cstr_ode,
        par_est.examples.cstr_dae,
        par_est.examples.cstr_dae_constant,
        par_est.examples.cstr_ode_constant,
    ]:
        variable_list, m = cstr_model(piecewise)
        # Create time-grid. Zero should be first
        time_grid = np.linspace(10, 100000, 4)
        time_grid = np.insert(time_grid, 0, 0)

        sim = par_est.Simulator(m, time_grid, variable_list)

        sim_res = sim.simulate()

        steady_state = sim.calculate_steady_state()

        assert np.isclose(sim_res["xf"][:,-1], steady_state[0:5]).all()

        if m.DAE:
            assert np.isclose(sim_res["zf"][:,-1], steady_state[5]).all()


if __name__ == "__main__":
    pass
    # test_pendulum_dae()
    # test_cstr()
    # test_dae_initials_calculation()
    # test_piecewise()
