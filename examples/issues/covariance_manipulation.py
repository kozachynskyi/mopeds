import casadi as ca
import numpy as np

import par_est
import par_est.examples

""" Jacobian returned by simulation has a size (NumState X NumOfParam * NumOfTimepoints). Once can get a covariance matrix of the parameters in different ways, and this example shows that something that works on paper (simple linear algebra) doesn't work with casadi.

It can have something to do with precision of linear algebra, but from theoretical perspective, covariance_difference_fail in the end should be zero.
Now it's close to zero, but not zero. This tells something about precision of mathematics in python / casadi / numpy.
"""

num_time = 10
varlist, model = par_est.examples.cstr_ode()

time_grid = np.linspace(0, 1000, num_time + 1)
for var in varlist.values():
    var.fixed = True
sim = par_est.Simulator(model, time_grid, varlist, simulate_jac=True)
res = sim.simulate_jac()

num_param = 19
num_state = 5
covariance_measurement = np.eye(num_state)

res_jacobian = res["jac_xf_p"]
covariance_full = res_jacobian.T @ covariance_measurement @ res_jacobian

covariance_measurement_reshape = np.kron(
    np.eye(num_time, dtype=int), covariance_measurement
)

split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param

list_jacobian_at_timepoint = ca.horzsplit(res_jacobian, split_vector)

jacobian_reshape = None
covariance_full_sum = None

for count_time_step in range(num_time):
    index_from = count_time_step * num_param
    index_till = num_param * count_time_step + (num_param)

    jac_at_timepoint = res_jacobian[:, index_from:index_till][:, 1:]
    cov_at_timepoint = jac_at_timepoint.T @ covariance_measurement @ jac_at_timepoint

    assert (
        list_jacobian_at_timepoint[count_time_step][:, 1:] - jac_at_timepoint
    ).is_zero()

    if jacobian_reshape is None:
        jacobian_reshape = jac_at_timepoint
    else:
        jacobian_reshape = ca.vertcat(jacobian_reshape, jac_at_timepoint)

    covariance_difference = (
        covariance_full[index_from + 1 : index_till, index_from + 1 : index_till]
        - cov_at_timepoint
    )
    assert covariance_difference.is_zero()

    if covariance_full_sum is None:
        covariance_full_sum = cov_at_timepoint
    else:
        covariance_full_sum = covariance_full_sum + cov_at_timepoint

covariance_reshape = (
    jacobian_reshape.T  # type: ignore
    @ covariance_measurement_reshape
    @ jacobian_reshape
)

covariance_difference_fail = covariance_full_sum - covariance_reshape

print(covariance_difference_fail)
