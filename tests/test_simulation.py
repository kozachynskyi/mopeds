import par_est
import numpy as np
from conftest import cstr_model_ode, pendulum_dae_1
import pytest
import casadi as ca


def test_tau():
    varlist, model = pendulum_dae_1()

    time_grid = np.linspace(0, 1, 40)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res_tau = sim.simulate()
    res = sim.integrator(
        x0=sim._initial_states,
        z0=sim._initial_alg,
        p=ca.vertcat(sim._variables * sim.scaling),
    )
    print(res_tau - ca.vertcat(res["xf"], res["zf"]))


def try_map():
    varlist, model = cstr_model_ode()
    time_grid = np.linspace(0, 5, 6)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)


    res0 = sim.simulate()
    accum = sim.integrator_tau.mapaccum('simulator', 2)

    res = accum(
            x0=sim._initial_states,
            z0=sim._initial_alg,
            p=ca.horzcat(ca.vertcat(ca.DM(1), sim._variables * sim.scaling), ca.vertcat(ca.DM(3), sim._variables * sim.scaling)),
            )
    breakpoint()


@pytest.mark.skip(reason="WIP")
def test_cstr_ode():
    variable_list, m = cstr_model_ode()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)

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
                res_simple, jac_simple = sim.simulate(True)

            res = sim.generate_exp_data()
            if j == 1:
                assert jac_simple.size() == (20, 19)
            assert res_simple.size() == (5, 4)
            assert len(res) == 5

            if i == 0:
                assert sim._variables[6].is_symbolic()
            elif i == 1:
                assert not sim._variables[6].is_symbolic()
                assert sim._variables[15].is_symbolic()
            elif i == 2:
                assert not sim._variables[15].is_symbolic()
            elif i == 3:
                assert sim._variables[10].is_symbolic()
            elif i == 4:
                assert not sim._variables[15].is_symbolic()


if __name__ == "__main__":
    try_map()
    # test_tau()
    # test_cstr_ode()
