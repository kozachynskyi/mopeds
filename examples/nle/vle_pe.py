import copy

import par_est
from par_est.examples import vle_nle_problem
from par_est.tools import generate_varlist_with_data_NLE

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

    # Fix all variables
    var_list_fixed.set_variable_list_fixed()

    # Preallocate list for optimizer_lists
    optimizer_list = []

    # Define changing mole fraction x
    x = [0.00, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00]

    for val in x:
        var_list_fixed["x"].dataframe.iloc[0] = val  # adjust the value for x

        # Generate data using the fixed values for the parameter. The output
        # is used as starting values for parameter estimation
        variable_list_optimizer = generate_varlist_with_data_NLE(model, var_list_fixed)

        # Unfix parameters for parameter estimation
        variable_list_optimizer.set_variable_list_unfixed(param_list)

        # Set upper and lower bounds
        variable_list_optimizer.set_bounds(val=0.01, emerg_val=50)

        # Append variable_list_optimizer to optimizer_list
        optimizer_list.append(variable_list_optimizer)

    # Define ParameterEstimationNLE
    pe = par_est.ParameterEstimationNLE(model, optimizer_list)

    # Solve parameter estimation problem
    res = pe.optimize(False)
