import copy

import mopeds
from mopeds.examples import vle_nle_problem
from mopeds.tools import generate_artificial_data_from_grid_nle

if __name__ == "__main__":
    variable_list, model = vle_nle_problem()

    param_list = [
        # "a1",
        # "b1",
        # "c1",
        "a2",
        # "b2",
        # "c2",
    ]  # List of parameters for parameter estimation

    # Set parameters and controls to fixed state so their values are used for
    # simulation
    var_list_fixed = copy.deepcopy(variable_list)
    var_list_fixed["T"].variance = (0.1) ** 2

    # Fix all variables
    var_list_fixed.set_variable_list_fixed()

    # Define changing mole fraction x
    x_bounds = {"x": [0, 1, 9]}

    variable_list_optimizer, true_parameters = generate_artificial_data_from_grid_nle(
        model, var_list_fixed, x_bounds
    )

    for varlist in variable_list_optimizer:
        varlist.set_variable_list_unfixed(param_list)

    # Define ParameterEstimationNLE
    pe = mopeds.ParameterEstimationNLE(model, variable_list_optimizer)

    # Solve parameter estimation problem
    res = pe.optimize(False)
