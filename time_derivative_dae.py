import copy

import numpy as np
from matplotlib import pyplot as plt

import par_est
import par_est.examples

plt.ion()
import pandas as pd
import casadi as ca

class SimulatorCustom(par_est.Simulator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        if self._Simulator__integrator_name == "collocation":
            self.integrator_tau_jac_zf = self.integrator_tau.factory(
                "integrator_tau_jacobian",
                self.integrator_tau.name_in(),
                ["xf", "zf", "qf", "rxf", "rzf", "rqf", "jac:xf:p", "jac:zf:p"],
            )
            self.simulate_jac_zf = self._simulate_jac_dae_zf

    def _simulate_jac_dae_zf(self):
        """Return dictionary with results "xf" - state,
        "zf" - algebraic, "jac_xf_p" - derivatives.
        """
        prev_time_step = 0
        res_states = []
        res_algebraic = []
        res_jacobian = []
        res_jacobian_zf = []
        x_init = self._initial_state
        alg_init = self._initial_algebraic

        for time_step, independent_variables in zip(
            self.time_grid_relative[1:], self._independent_variables
        ):
            res_integration = self.integrator_tau_jac_zf(
                x0=x_init,
                z0=alg_init,
                p=ca.vertcat(
                    time_step - prev_time_step, independent_variables * self.scaling
                ),
            )

            prev_time_step = time_step
            x_init = res_integration["xf"]
            alg_init = res_integration["zf"]

            res_states.append(res_integration["xf"])
            res_algebraic.append(res_integration["zf"])
            res_jacobian.append(res_integration["jac_xf_p"])
            res_jacobian_zf.append(res_integration["jac_zf_p"])

        res_states = ca.hcat(res_states)
        res_algebraic = ca.hcat(res_algebraic)
        res_jacobian = ca.hcat(res_jacobian)
        res_jacobian_zf = ca.hcat(res_jacobian_zf)

        res = {"xf": res_states, "zf": res_algebraic, "jac_xf_p": res_jacobian, "jac_zf_p": res_jacobian_zf}
        return res


if __name__ == "__main__":

    piecewiseswitch = False
    variable_list, m = par_est.examples.cstr_dae(piecewiseswitch)

    # Create time-grid. Zero should be first
    time_grid = np.linspace(0, 1000, 5)
    time_grid2 = np.linspace(0, 1002, 5)

    # Create simulation Object
    sim_fixed = SimulatorCustom(
    # sim_fixed = par_est.Simulator(
        m, time_grid, variable_list, use_idas_constraints=False, simulate_jac=True
        , integrator_name="collocation"
    )
    res = sim_fixed.generate_exp_data(algebraic=True)

    sim_fixed2 = SimulatorCustom(
    # sim_fixed2 = par_est.Simulator(
        m, time_grid2, variable_list, use_idas_constraints=False, simulate_jac=True
    )
    res2 = sim_fixed2.generate_exp_data(algebraic=True)

    # jac = sim_fixed.simulate_jac_zf()["jac_xf_p"]
    jac_zf = sim_fixed.simulate_jac_zf()["jac_zf_p"]
    jac = sim_fixed.simulate_jac_zf()["jac_xf_p"]
    r = sim_fixed.simulate_jac_zf()
    for i, (t1, t2) in enumerate(zip(time_grid, time_grid2)):
        if i == 0:
            pass
        else:
            print((res["e0_T"].value[i] - res2["e0_T"].value[i])/(t1-t2))
            print(jac[0,19*(i-1)])

            print((res["e0_c_tot"].value[i] - res2["e0_c_tot"].value[i])/(t1-t2))
            print(jac_zf[0,19*(i-1)])
