import numpy as np
import logging
import mopeds.examples
import casadi as ca


def test_vle_nle():
    for i in [True, False]:
        mopeds.set_options(variable_scaling=i)
        variable_list, model = mopeds.examples.vle_nle_problem()

        variable_list.set_variable_list_fixed()
        variable_list["x"].value = 0.5
        sim = mopeds.SimulatorNLE(model, variable_list)
        res = sim.simulate_sym()
        true_answer_T = 359.451

        variable_list["a1"].fixed = False
        sim_unfixed = mopeds.SimulatorNLE(model, variable_list)
        res_unfixed = sim_unfixed.simulate_sym_unfixed({"a1": 5.24125})

        logging.warning(
            f"Model.NLE: {model}, Result: {res['x']}, Expecting: {true_answer_T}"
        )
        assert np.isclose(res["x"], ca.DM(true_answer_T), rtol=0, atol=1.0e-3)
        assert np.isclose(res_unfixed["x"], ca.DM(true_answer_T), rtol=0, atol=1.0e-3)


def test_utilities_methods():
    variable_list, model = mopeds.examples.simple_mixer()

    sim = mopeds.SimulatorNLE(model, variable_list)
    res_full = sim.simulate()

    assert res_full[0] == sim.simulate(["e0_F_s2"])
    assert res_full[1] == sim.simulate(["e0_F_s4"])

    sim.change_independent_variables({"e0_F_s1": 10, "e0_F_s3": 15})
    res_full = sim.simulate()
    assert np.array_equal(res_full, ca.DM([10.0, -5.0, 10.0]))


if __name__ == "__main__":
    pass
    test_vle_nle()
    # test_utilities_methods()
