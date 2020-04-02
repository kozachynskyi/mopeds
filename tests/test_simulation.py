import par_est
import numpy as np
from conftest import cstr_model_ode
import pytest


@pytest.mark.skip(reason="WIP")
def test_cstr_ode():
    variable_list, m = cstr_model_ode()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 4)
    time_grid = np.insert(time_grid, 0, 0)

    for i in range(5):
        for j in range(2):
            if i == 1:
                variable_list["e0_U"].fixed = True
            elif i == 2:
                variable_list["e0_T_in"].fixed = True
            elif i == 3:
                variable_list["e0_c_i2"].fixed = True
            elif i == 4:
                for var in variable_list.values():
                    var.fixed = True

            sim = par_est.Simulator(m, time_grid, variable_list)
            if j == 0:
                res_simple = sim.simulate()
            else:
                res_simple, jac_simple = sim.simulate(True)

            res = sim.generate_exp_data()
            if j == 1:
                assert jac_simple.size() == (20, 19)
            assert res_simple.size() == (5, 4)
            assert len(res) == 5

            if i == 0:
                assert sim._variables[6].is_symbolic()
            elif i == 1:
                assert not sim._variables[6].is_symbolic()
                assert sim._variables[15].is_symbolic()
            elif i == 2:
                assert not sim._variables[15].is_symbolic()
            elif i == 3:
                assert sim._variables[10].is_symbolic()
            elif i == 4:
                assert not sim._variables[15].is_symbolic()


if __name__ == "__main__":
    test_cstr_ode()
