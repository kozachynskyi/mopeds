import logging
import casadi as ca
import numpy as np

import copy
import par_est.examples
import par_est
import par_est.tools


def test_pe():
    """Test that ParameterEstimationNLE on NLE always yields same result.
    Helpfull to see if any drastic changes in calculation were made
    """
    variable_list, model = par_est.examples.vle_nle_problem()

    var_list_fixed = copy.deepcopy(variable_list)
    var_list_fixed.set_variable_list_fixed()

    var_list_fixed["x"].value = 0.5
    variable_list_optimizer = par_est.tools.generate_varlist_with_data_NLE(
        model, var_list_fixed
    )
    variable_list_optimizer.set_variable_list_unfixed(["a2"])
    variable_list_optimizer.set_bounds(emerg_val=50)

    for i in range(3):
        if i == 0:
            simulator_name = "rootfinder"
        else:
            simulator_name = "ipopt"

        if i == 2:
            for var in variable_list_optimizer.values():
                if isinstance(var, par_est.VariableAlgebraic):
                    var.lower_bound = 350
                    var.upper_bound = 380

        pe = par_est.optimization.ParameterEstimationNLE(
            model,
            [variable_list_optimizer, variable_list_optimizer],
            simulator_name=simulator_name,
        )

        for switch in [True, False]:
            res = pe.optimize(switch)
            answer_f = 0
            answer_param = [5.19625]

            logging.warning(
                f"Model.NLE: {model}, Result: {res['f']}, Expecting: {answer_f}"
            )
            assert np.isclose(res["f"], ca.DM(answer_f), rtol=0, atol=1.0e-9)

            logging.warning(
                f"Model.NLE: {model}, Result: {res['x']}, Expecting: {answer_param}"
            )
            assert np.all(np.isclose(res["x"], ca.DM(answer_param), rtol=0, atol=1.0e-9))

        for sim in pe.list_simulators:
            if i == 1:
                assert sim._lower_bound[0] == -ca.inf
                assert sim._upper_bound[0] == ca.inf
                assert sim._rootfinder_bounds[0] == 0
            if i == 2:
                assert sim._lower_bound[0] == 350
                assert sim._upper_bound[0] == 380
                assert sim._rootfinder_bounds[0] == 2


if __name__ == "__main__":
    pass
    test_pe()
