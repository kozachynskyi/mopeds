import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples
import casadi as ca
from warnings import warn
from functools import wraps, partial, cached_property

plt.ion()


class PE_GN(mopeds.ParameterEstimation):
    def optimize(
        self,
        scale=None,
        objective_function="wls",
        direct_optimization=False,
        *,
        reuse_solver=False,
        gauss_newton=False,
    ):
        if objective_function == "wls":
            self._objective = partial(
                self._objective_wls, direct_optimization=direct_optimization
            )
        elif objective_function == "ols":
            self._objective = partial(
                self._objective_ols, direct_optimization=direct_optimization
            )
        elif objective_function == "fair":
            self._objective = partial(
                self._objective_fair, direct_optimization=direct_optimization
            )
        elif objective_function == "tikh":
            self._objective = partial(
                self._objective_tikhonov, direct_optimization=direct_optimization
            )
        else:
            raise NotImplementedError(
                f"Objective function '{objective_function}' is not supported"
            )

        return self._optimize(
            scale,
            direct_optimization=direct_optimization,
            reuse_solver=reuse_solver,
            gauss_newton=gauss_newton,
        )

    def _optimize(
        self,
        scale: bool = None,
        direct_optimization: bool = False,
        *,
        reuse_solver: bool = False,
        gauss_newton=False,
    ) -> dict[str, ca.DM | ca.MX]:
        """Runs optimizer, uses scaling if needed. Returned values is scaled back."""
        if scale is not None:
            warn("Scale argument is deprecated", FutureWarning, 5)

        if direct_optimization:
            varlist_decision = self.varlist_decision_direct
            if gauss_newton:
                raise NotImplementedError
        else:
            varlist_decision = self.varlist_decision

        self.nlpsol_dict = {
            "x": varlist_decision.get_casadi_variables(),
            "f": self._objective()[0],
            "p": self._nlpsol_p_mx,
        }

        if direct_optimization:
            self.nlpsol_dict["g"] = self.nlpsol_g_direct

        # source: casadi test/python/nlp.py  test_gauss_newton_sqpmethod Commit: 3d820e62cb588e
        if gauss_newton:
            x = self.varlist_decision.get_casadi_variables()
            F = self._unscale_residuals(self._objective()[1])
            J = ca.jacobian(F, x)
            p = self._nlpsol_p_mx
            lam_f = ca.MX.sym("x")
            lam_g = ca.MX.sym("x", 0, 1)
            GN = ca.Function(
                "GN", [x, p, lam_f, lam_g], [lam_f * ca.triu(2 * ca.mtimes(J.T, J))]
            )

            self.solver_settings["hess_lag"] = GN
        else:
            self.solver_settings.pop("hess_lag", None)

        print(self.solver_settings)
        if not (hasattr(self, "solver") and reuse_solver):
            self.solver: ca.Function = ca.nlpsol(
                "solver",
                self.solver_name,
                self.nlpsol_dict,
                self.solver_settings,
            )

        if direct_optimization:
            self.nlpsol_args = {
                "x0": self.guess_direct,
                "lbx": self.lower_bound_direct,
                "ubx": self.upper_bound_direct,
                "lbg": [0] * self.nlpsol_g_direct.shape[0],
                "ubg": [0] * self.nlpsol_g_direct.shape[0],
            }

        else:
            self.nlpsol_args = {
                "x0": self.guess,
                "lbx": self.lower_bound,
                "ubx": self.upper_bound,
            }
        self.nlpsol_args["p"] = self._nlpsol_p_values

        res_solver = self.solver.call(self.nlpsol_args)

        res_solver["x_unscaled"] = res_solver["x"].toarray()
        res_solver["x_all"] = np.asarray(
            varlist_decision.scale_to_original(res_solver["x"])
        )
        if direct_optimization:
            res_solver["x"] = res_solver["x_all"][: len(self.varlist_decision)]
        else:
            res_solver["x"] = res_solver["x_all"]

        res_dict = {}
        for solution, var_name in zip(
            res_solver["x"], list(self.varlist_decision.keys())
        ):
            res_dict[var_name] = float(solution[0])

        if direct_optimization:
            res_dict_all = {}
            for solution, var_name in zip(
                res_solver["x_all"], list(varlist_decision.keys())
            ):
                res_dict_all[var_name] = float(solution[0])
        else:
            res_dict_all = res_dict

        res_solver["x_dict"] = res_dict
        res_solver["x_dict_all"] = res_dict
        self.reset_acados()

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
    # pe.solver_settings["ipopt"]["hessian_approximation"] = "exact"
    # print(pe.optimize())
    print(pegn.optimize(gauss_newton=True)["x"])
    print(pegn.optimize(gauss_newton=False)["x"])
    # pe.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"
    # print(pe.optimize())
