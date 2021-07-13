import numpy as np

import par_est
import par_est.examples

if __name__ == "__main__":

    # Model setup
    piecewiseswitch = True
    variable_list, m = par_est.examples.cstr_ode(piecewiseswitch)
    for var in variable_list.values():
        var.fixed = True

    variable_list["e0_U"].fixed = False
    variable_list["e0_T_in"].fixed = False
    variable_list["e0_T"].variance = 0.1

    time_grid1 = np.linspace(0, 1000, 4)

    # Here I setup horizon for artificial data
    if piecewiseswitch:
        variable_list["e0_T_in"].expand_horizon([241., 482., 723.], [363, 373, 383])

    # Here data is generated
    data1 = par_est.tools.generate_varlist_with_data(variable_list, m, time_grid1)

    # Now I need to reset a "e0_T_in" variable! Else - error raised
    # Value will be used as a guess for every horizon
    # Lower and upper bound are set also for all horizons
    data1["e0_T_in"] = par_est.VariableControlPiecewiseConstant("e0_T_in", 373.0, 353.0, 393.0)
    data1["e0_T_in"].fixed = False

    # Time_grid is created based on available data
    # Than divided in 4 equal horizonts and for each of this
    # Horizons a new decision variable is created
    # Tested only for one VariableControlPiecewiseConstant

    mpc = par_est.ModelPredictiveControl(m, [data1], number_of_time_horizonts=4)
    print(f"Optimizer calculates control variables at time:\n {mpc.time_grid_controls}")

    # Because optimizer is formulated as Single Shooting, ipopt is not efficient. So it may get stuck and may require many iterations.
    mpc.solver_settings["ipopt"]["max_iter"] = 50

    res = mpc.optimize()
    print(res)
