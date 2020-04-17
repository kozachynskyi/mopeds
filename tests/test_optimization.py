import par_est
import logging
import casadi as ca
import numpy as np
import pytest

import copy
from conftest import cstr_model_ode, pendulum_dae_1
import logging


@pytest.mark.skip(reason="WIP")
def test_ode():
    var_list, model = cstr_model_ode()
    time_grid = np.linspace(10, 10000, 3)
    time_grid = np.insert(time_grid, 0, 0)

    x = np.array([0,1,2,3])
    xx = np.array([[0,1,2,3],[4,5,6,7]])
    xxx = np.vstack((xx,[8,8,8,8]))
    y = np.array([0,1])
    z = np.array([2,3])
    mask_y = np.isin(x,y)
    mask_z = np.isin(x,z)
    maskmask = np.vstack((mask_y,mask_z))
    masknonzero = np.nonzero(maskmask)
    print(xx.shape)
    print(maskmask.shape)
    print(xx[maskmask])

    M = ca.SX([[0,1,2,3],[4,5,6,7]])
    print(M.get(True,ca.Slice(0,2),ca.Slice(0,3)))
    exit()

    var_list_fixed = copy.deepcopy(var_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(
        model, time_grid, var_list_fixed
    ).generate_exp_data()

    for key, var in var_list_exp.items():
        var_list[key] = var

    for var in var_list.values():
        var.fixed = True

    var_list["e0_E_r1"].fixed = False

    var_list["e0_T"].value.value = var_list["e0_T"].value.value[0:2]
    var_list["e0_T"].value.time = var_list["e0_T"].value.time[0:2]
    # var_list["e0_c_i1"].value.value = var_list["e0_c_i1"].value.value[0:2]
    # var_list["e0_c_i1"].value.time = var_list["e0_c_i1"].value.time[0:2]

    pe = par_est.ParameterEstimation(model, [var_list])

    res = pe.optimize()
    logging.warning(f"{res['f']}")
    assert np.isclose(res["f"], ca.DM(1.93785816e-15), rtol=0, atol=1.0e-21)

    res = pe.optimize(False)
    logging.warning(f"{res['f']}")
    assert np.isclose(res["f"], ca.DM(1.78126185e-11), rtol=0, atol=1.0e-17)

    var_list["e0_c_in_i1"].fixed = False

    oed = par_est.OptimalExperimentalDesign(model, [var_list], time_grid)
    res = oed.optimize()
    fim = oed.get_fim_matrix()

    assert fim[0].size() == (20, 1)
    assert fim[1].size() == (1, 1)

    logging.warning(f"{res['f']}")
    assert np.isclose(res["f"], ca.DM(45.16749745), rtol=0, atol=1.0e-8)

    res = oed.optimize(False)

    logging.warning(f"{res['f']}")
    assert np.isclose(res["f"], ca.DM(45.16749745), rtol=0, atol=1.0e-8)


def not_test_optimizer():
    variable_list, m = cstr_model_ode()

    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 1000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    # Generate experimental data for parameter estimation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True
    var_list_exp = par_est.Simulator(m, time_grid, var_list_fixed).get_symbolic_result()

    # Replace empty state variables with results from simulation
    for key, var in var_list_exp.items():
        variable_list[key] = var

    for i in range(5):
        if i == 1:
            variable_list["e0_U"].fixed = True
        elif i == 2:
            variable_list["e0_T_in"].fixed = True
        elif i == 3:
            variable_list["e0_c_i2"].fixed = True
        elif i == 4:
            for var in variable_list.values():
                var.fixed = True
            variable_list["e0_U"].fixed = False
            variable_list["e0_T_in"].fixed = False

        for j in range(2):
            if j == 0:
                pe = par_est.ParameterEstimation(m, [variable_list])
            else:
                pe = par_est.ParameterEstimation(m, [variable_list, variable_list])

            pe.solver_settings = {
                "verbose": False,
                "ipopt": {"max_iter": 1},
            }

            oed = par_est.OptimalExperimentalDesign(m, [variable_list], time_grid)
            oed.solver_settings = {
                "verbose": False,
                "ipopt": {"hessian_approximation": "limited-memory", "max_iter": 1},
            }

            for ij in range(2):
                if ij == 0:
                    res_pe = pe.optimize()
                    if j == 0:
                        res_oed = oed.optimize()
                else:
                    res_pe = pe.optimize(False)
                    if j == 0:
                        res_oed = oed.optimize(False)

            oed.get_fim_matrix()

            if i == 0:
                assert res_pe["x"].size() == (11, 1)
                assert res_oed["x"].size() == (7, 1)
            elif i == 1:
                assert res_pe["x"].size() == (10, 1)
                assert res_oed["x"].size() == (7, 1)
            elif i == 2:
                assert res_pe["x"].size() == (10, 1)
                assert res_oed["x"].size() == (6, 1)
            elif i == 3:
                assert res_pe["x"].size() == (10, 1)
                assert res_oed["x"].size() == (6, 1)
            elif i == 4:
                assert res_pe["x"].size() == (1, 1)
                assert res_oed["x"].size() == (1, 1)


if __name__ == "__main__":
    logging.basicConfig(format='%(name)s:%(levelname)s:%(message)s', level=logging.DEBUG)
    test_ode()
    # not_test_optimizer()
