import numpy as np
import pytest

import par_est


def test_artificial_data_generator_nle():
    variable_list, model = par_est.examples.simple_mixer()
    control_bounds = {"e0_F_s1": [17, 20, 3]}
    res, par = par_est.tools.generate_varlist_with_data_NLE(model, variable_list, control_bounds, perturbate=False)

    for varlist, expected in zip(res, [17, 18.5, 20]):
        assert varlist["e0_F_s2"].value[0] == expected

    for varlist, expected in zip(res, [4, 5.5, 7]):
        assert varlist["e0_F_s4"].value[0] == expected

    res, par = par_est.tools.generate_varlist_with_data_NLE(model, variable_list, control_bounds, perturbate=True, rng=np.random.default_rng(0), measurement_names=["e0_F_s2"])

    for varlist, expected in zip(res, [17.125730221093395, 18.367895136708697, 20.64042265044328]):
        assert np.isclose(varlist["e0_F_s2"].value[0], expected)
        assert np.isnan(varlist["e0_F_s4"].value[0])

def test_startpoint_generation():
    for method in ["lhs", "hammersley"]:
        bound = np.array([[1e-3, 1e3]])
        par_est.utilities.make_startpoints(bound, 10, method)
        with pytest.raises(ValueError):
            bound = np.array([[0, 1e3]])
            par_est.utilities.make_startpoints(bound, 10)


if __name__ == "__main__":
    test_artificial_data_generator_nle()
    test_startpoint_generation()
