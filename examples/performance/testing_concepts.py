import casadi as ca
import numpy as np

import mopeds
import mopeds.examples
import mopeds.tools


def jacobian_manipulation():
    """Vizualizes how OED manipulates jacobian to get covariance."""

    num_time = 3
    varlist, model = mopeds.examples.cstr_ode()

    time_grid = np.linspace(0, 1000, num_time + 1)
    for var in varlist.values():
        var.fixed = True
    sim = mopeds.Simulator(model, time_grid, varlist, simulate_jac=True)
    res = sim.simulate_jac()

    num_param = 19

    res_jacobian = res["jac_xf_p"]

    index_all_states = list(range(len(model.varlist_state)))

    split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param

    list_jacobian_at_timepoint = ca.horzsplit(res_jacobian, split_vector)

    for jacobian in list_jacobian_at_timepoint:
        selected = jacobian.get(False, index_all_states, [0, 2])
        mopeds.utilities.plot_arrays([jacobian, selected])


def covariance_matrix_summation():
    for i in range(10):
        A_1 = np.random.rand(25, 2)
        # A_2 = np.random.rand(5,2)
        A_2 = np.zeros((5, 2))
        A = np.vstack([A_1, A_2])

        COV_1 = A_1.T @ A_1
        COV_2 = A_2.T @ A_2
        COV = A.T @ A
        COV_ALL = COV_1 + COV_2

        print("COV_1:\n", COV_1 / COV)
        # print(COV_2 / COV)


if __name__ == "__main__":
    # jacobian_manipulation()
    covariance_matrix_summation()
