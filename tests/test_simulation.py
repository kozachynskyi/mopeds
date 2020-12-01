import par_est
import numpy as np
import logging
import par_est.examples
import casadi as ca


def test_pendulum_dae():
    varlist, model = par_est.examples.pendulum_dae_1()

    time_grid = np.linspace(0, 1, 3)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res_tau = sim.simulate()
    res = sim.integrator(
        x0=sim._initial_state,
        z0=sim._initial_algebraic,
        p=ca.vertcat(sim._variables * sim.scaling),
    )

    assert np.isclose(
        ca.vertcat(res_tau["xf"], res_tau["zf"]), ca.vertcat(res["xf"], res["zf"])
    ).any()

def test_vle_nle():
    variable_list, model = par_est.examples.vle_nle_problem()

    variable_list.set_variable_list_fixed()
    variable_list['x'].value = 0.5 
    sim = par_est.simulation.SimulatorNLE(model, variable_list)
    res = sim.simulate_sym()
    true_answer_T = 359.451

    logging.warning(
            f"Model.NLE: {model}, Result: {res['r']}, Expecting: {true_answer_T}"
        )
    assert np.isclose(res["r"], ca.DM(true_answer_T), rtol=0, atol=1.0e-3)


def test_cstr():
    for cstr_model in [par_est.examples.cstr_ode, par_est.examples.cstr_dae, par_est.examples.cstr_dae_constant, par_est.examples.cstr_ode_constant]:
        variable_list, m = cstr_model()
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
                    res_simple = sim.simulate(True)

                res = sim.generate_exp_data()
                if j == 1:
                    assert res_simple["jac_xf_p"].size() == (5, 76)
                assert res_simple["xf"].size() == (5, 4)
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
    test_pendulum_dae()
    test_cstr()
    test_vle_nle()
