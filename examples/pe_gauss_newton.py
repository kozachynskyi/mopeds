import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples
import casadi as ca

plt.ion()

class PE_GN(mopeds.ParameterEstimation):
    def _optimize(self, scale: bool = None, direct_optimization: bool = False, *, reuse_solver: bool = False) -> dict[str, ca.DM | ca.MX]:

        """Runs optimizer, uses scaling if needed. Returned values is scaled back.
        Scaling should be done before setting a solver and solver settings."""

        x = self.varlist_decision.get_casadi_variables()
        F = self._objective()[1]
        J = ca.jacobian(F, x)
        p = ca.MX.sym("x",0,1)
        lam_f = ca.MX.sym("x")
        lam_g = ca.MX.sym("x",0,1)
        GN = ca.Function('GN',[x,p,lam_f,lam_g],[lam_f*ca.triu(2*ca.mtimes(J.T,J))])

        self.solver_settings["hess_lag"] = GN
        # self.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"


        self.solver: ca.Function = ca.nlpsol(
            "solver",
            self.solver_name,
            {
                "x": self.varlist_decision.get_casadi_variables(),
                "f": self._objective()[0],
                "p": self._nlpsol_p_mx,

            },
            self.solver_settings,
        )

        lb_scaled = self.lower_bound
        ub_scaled = self.upper_bound

        # Scaling of negative numbers requires a switch bounds
        for index, (lb, ub) in enumerate(zip(lb_scaled, ub_scaled)):
            if lb > ub:
                lb_scaled[index] = ub
                ub_scaled[index] = lb

        res_solver = self.solver(
            x0=self.guess,
            lbx=lb_scaled,
            ubx=ub_scaled,
            p=self._nlpsol_p_values
        )

        res_solver["x"] = res_solver["x"]

        res_dict = {}
        for solution, var_name in zip(
            res_solver["x"].toarray(), list(self.varlist_decision.keys())
        ):
            res_dict[var_name] = float(solution[0])

        res_solver["x_dict"] = res_dict

        return res_solver


if __name__ == "__main__":
    mopeds.set_options(variable_scaling=False)

    piecewiseswitch = False
    variable_list, m = mopeds.examples.cstr(piecewiseswitch)
    for var in variable_list.values():
        var.fixed = True

    variable_list["e0_U"].fixed = False
    variable_list["e0_c_p"].fixed = True
    # variable_list["e0_E_r1"].fixed = False
    variable_list["e0_T_in"].fixed = False
    variable_list["e0_T"].variance = 1

    # Create time-grid. Zero should be first
    time_grid1 = np.linspace(0, 100, 400)
    time_grid2 = np.linspace(0, 1000, 8)

    e0_T_in = variable_list["e0_T_in"]
    variable_list["e0_F"].fixed = False
    if isinstance(e0_T_in, mopeds.VariableControlPiecewiseConstant):
        e0_T_in.expand_horizon([10, 723], [363, 453])
        e0_T_in.variable_list.index(0).fixed = False
        e0_T_in.variable_list.index(1).fixed = True
        e0_T_in.variable_list.index(2).fixed = False

    data1 = mopeds.tools.generate_varlist_with_data(variable_list, m, time_grid1, True)
    data2 = mopeds.tools.generate_varlist_with_data(variable_list, m, time_grid2, True)
    # data1.show()

    # If data is not available for all simulated points, PE works
    e0_T = data2["e0_T"]
    e0_T.dataframe = e0_T._dataframe_from_value(e0_T.value[0])

    e0_c_i1_df = data2["e0_c_i1"].dataframe
    e0_c_i1_df.drop(e0_c_i1_df.index[2:], inplace=True)

    # Perturbate alg variable from "ideal solution" to see its affect on PE
    a = data2["e0_c_tot"].dataframe
    data2["e0_c_tot"].dataframe = data2["e0_c_tot"].dataframe * 1.05

    # pe_state = mopeds.ParameterEstimation(m, [data1, data2])
    pe = mopeds.ParameterEstimation(m, [data1])
    pegn = PE_GN(m, [data1])
    # a = pe_state.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_c_p": 3.5, "e0_E_r1": 9.6e4})
    pe.solver_settings["ipopt"]["hessian_approximation"] = "exact"
    print(pe.optimize())
    print(pegn.optimize())
    pe.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"
    print(pe.optimize())
