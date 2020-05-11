import par_est.examples
import par_est
import casadi as ca
import numpy as np
import par_est.tools


def jacobian_manipulation():
    """ Vizualizes how OED manipulates jacobian to get covariance. """

    num_time = 3
    varlist, model = par_est.examples.cstr_ode()

    time_grid = np.linspace(0, 1000, num_time + 1)
    for var in varlist.values():
        var.fixed = True
    sim = par_est.Simulator(model, time_grid, varlist)
    res = sim.simulate(True)

    num_param = 19

    res_jacobian = res["jac_xf_p"]

    index_all_states = list(range(len(model.varlist_state)))

    split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param

    list_jacobian_at_timepoint = ca.horzsplit(res_jacobian, split_vector)

    for jacobian in list_jacobian_at_timepoint:
        selected = jacobian.get(False, index_all_states, [0, 2])
        par_est.tools.plot_arrays([jacobian, selected])


if __name__ == "__main__":
    jacobian_manipulation()
