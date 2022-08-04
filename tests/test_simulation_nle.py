import numpy as np
import logging
import par_est.examples
import casadi as ca


def test_vle_nle():
    variable_list, model = par_est.examples.vle_nle_problem()

    variable_list.set_variable_list_fixed()
    variable_list["x"].value = 0.5
    for i in range(2):
        if i == 0:
            sim = par_est.simulation.SimulatorNLE(model, variable_list)
        else:
            sim = par_est.simulation.SimulatorNLE(
                model, variable_list
            )
        res = sim.simulate_sym()
        true_answer_T = 359.451

        logging.warning(
            f"Model.NLE: {model}, Result: {res['x']}, Expecting: {true_answer_T}"
        )
        assert np.isclose(res["x"], ca.DM(true_answer_T), rtol=0, atol=1.0e-3)


def test_utilities_methods():
    variable_list, model = par_est.examples.simple_mixer()

    sim = par_est.SimulatorNLE(model, variable_list)
    res_full = sim.simulate()

    assert res_full[0] == sim.simulate(["e0_F_s2"])
    assert res_full[1] == sim.simulate(["e0_F_s4"])

    sim.change_independent_variables({"e0_F_s1": 10, "e0_F_s3": 15})
    res_full = sim.simulate()
    assert np.array_equal(res_full, ca.DM([10., -5., 10.]))


if __name__ == "__main__":
    pass
    # test_vle_nle()
    # test_utilities_methods()
