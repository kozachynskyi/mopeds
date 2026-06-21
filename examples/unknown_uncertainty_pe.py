import numpy as np
from matplotlib import pyplot as plt

import mopeds
import mopeds.examples
import casadi as ca

plt.ion()


class PE(mopeds.ParameterEstimation):
    def _objective_ols(self):
        """Objective function is Bayesian with Unknown general covariance Bard Page65"""
        residuals = (self.simulate_all_mx - self.array_data) * self.array_data_mask
        num_exp = len(self.list_simulators)
        z_matrix = residuals.T @ residuals
        objective = (num_exp / 2) * ca.log(
            z_matrix[0, 0] * z_matrix[1, 1] - z_matrix[1, 0] * z_matrix[0, 1]
        )
        return objective, residuals


if __name__ == "__main__":
    piecewiseswitch = False
    variable_list, m = mopeds.examples.cstr_dae(piecewiseswitch)
    for var in variable_list.values():
        var.fixed = True

    variable_list["e0_U"].fixed = False
    variable_list["e0_c_p"].fixed = True
    # variable_list["e0_E_r1"].fixed = False
    variable_list["e0_T_in"].fixed = False
    variable_list["e0_T"].variance = 1

    # Create time-grid. Zero should be first
    time_grid1 = np.linspace(0, 1000, 4)
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
    pe_state = PE(m, [data1])
    # a = pe_state.calculate_sensitivity_and_fim({"e0_U": 1.4, "e0_c_p": 3.5, "e0_E_r1": 9.6e4})
    print(pe_state.optimize(True, "ols"))
    print(pe_state.optimize(True, "wls"))
