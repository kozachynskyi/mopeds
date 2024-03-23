import mopeds
import pandas as pd
import numpy as np


if __name__ == "__main__":
    mopeds.set_options(variable_scaling=True)
    varlist_original , model = mopeds.examples.vle_wilson()

    DIRECT_FAILS = True
    if DIRECT_FAILS:
        controls = {"e0_P": [1,3,2], "e0_x_c1": [0.02, 0.99, 3]}
    else:
        controls = {"e0_P": [1,3,2], "e0_x_c1": [0.02, 0.99, 3]}

    sim = mopeds.SimulatorNLE(model, varlist_original)
    res = sim.simulate()[2]

    exp_data, true_parameters = mopeds.tools.generate_artificial_data_from_grid_nle(model, varlist_original, control_bounds=controls, measurement_names=["e0_T", "e0_y_c1"], perturbate=False)

    if True:
        parameters = ["e0_greek_lambdaA_c1_j2", "e0_greek_lambdaA_c2_j1"]
        parameters_guess = {"e0_greek_lambdaA_c1_j2": 0.02, "e0_greek_lambdaA_c2_j1": -0.5}
        for vl in exp_data:
            for par_name in parameters:
                vl[par_name].fixed = False
                vl[par_name].guess = parameters_guess[par_name]

    if False:
        controls_scaling = {"e0_P": [1,3,2], "e0_x_c1": [0.02, 0.99, 5]}
        mopeds.tools.analyze_scaling_nle(model, varlist_original, controls_scaling)

    pe = mopeds.ParameterEstimationNLE(model, exp_data)
    res = pe.optimize(direct_optimization=True)
    v = pe.check_result_bounds(res)
    print(res["x_dict"])

