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


if __name__ == "__main__":
    pass
    # test_vle_nle()
