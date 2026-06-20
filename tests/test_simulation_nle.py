import numpy as np
import logging
import mopeds.examples
import casadi as ca


def test_vle_nle():
    for solver in ["rootfinder", "nlpsol", "newton", "fast_newton", "ipopt"]:
        for i in [True, False]:
            mopeds.set_options(variable_scaling=i)
            variable_list, model = mopeds.examples.vle_nle_problem()

            variable_list.set_variable_list_fixed()
            variable_list["x"].value = 0.5
            sim = mopeds.SimulatorNLE(model, variable_list, solver_name=solver)
            res = sim.simulate_fast()
            true_answer_T = 359.451

            variable_list["a1"].fixed = False
            sim_unfixed = mopeds.SimulatorNLE(model, variable_list, solver_name=solver)
            res_unfixed = sim_unfixed.simulate(unfixed_variables={"a1": 5.24125})[0]

            logging.warning(
                f"Model.NLE: {model}, Result: {res['x']}, Expecting: {true_answer_T}"
            )
            assert np.isclose(res["x"], ca.DM(true_answer_T), rtol=0, atol=1.0e-3)
            assert np.isclose(res_unfixed["x"], ca.DM(true_answer_T), rtol=0, atol=1.0e-3)


def test_utilities_methods():
    variable_list, model = mopeds.examples.simple_mixer()

    sim = mopeds.SimulatorNLE(model, variable_list)
    res_full = sim.simulate()[0]

    assert res_full["x"][0] == sim.simulate(return_var_names=["e0_F_s2"])[1]
    assert res_full["x"][1] == sim.simulate(return_var_names=["e0_F_s4"])[1]

    sim.change_independent_variables({"e0_F_s1": 10, "e0_F_s3": 15})
    res_full = sim.simulate()[0]["x"]
    assert np.array_equal(res_full, ca.DM([10.0, -5.0, 10.0]))


if __name__ == "__main__":
    pass
    test_vle_nle()
    # test_utilities_methods()
