import numpy as np
import pytest

import mopeds


def test_artificial_data_generator_nle():
    variable_list, model = mopeds.examples.simple_mixer()
    control_bounds = {"e0_F_s1": [17, 20, 3]}
    res, par, _ = mopeds.tools.generate_artificial_data_from_grid_nle(model, variable_list, control_bounds, perturbate=False)

    for varlist, expected in zip(res, [17, 18.5, 20]):
        assert varlist["e0_F_s2"].value[0] == expected

    for varlist, expected in zip(res, [4, 5.5, 7]):
        assert varlist["e0_F_s4"].value[0] == expected

    res, par, _ = mopeds.tools.generate_artificial_data_from_grid_nle(model, variable_list, control_bounds, perturbate=True, rng=np.random.default_rng(0), measurement_names=["e0_F_s2"])

    for varlist, expected in zip(res, [17.125730221093395, 18.367895136708697, 20.64042265044328]):
        assert np.isclose(varlist["e0_F_s2"].value[0], expected)
        assert np.isnan(varlist["e0_F_s4"].value[0])

def test_startpoint_generation():
    for method in ["lhs", "hammersley"]:
        bound = np.array([[1e-3, 1e3]])
        mopeds.utilities.make_startpoints(bound, 10, method)
        with pytest.raises(ValueError):
            bound = np.array([[0, 1e3]])
            mopeds.utilities.make_startpoints(bound, 10)

@pytest.mark.parametrize("scaling", [True, False])
def test_error_analyzer(scaling):
    """Test if df_params and list_predictions is correctly calculated, independently of the scaling"""
    with mopeds.options(variable_scaling=scaling):
        vl_original, model = mopeds.examples.linear_example()
        vl_original["y"].variance = 0.5**2
        # prediction_grid = {"u": [0, 1, 20], "v": [3, 4, 20]}
        prediction_grid = {"u": [0, 1, 3], "v": [3, 4, 1]}
        meas_grid = {"u": [0, 1, 2], "v": [3, 4, 2]}

        selected_parameters = ["a", "b", "c", "d"]
        rng = np.random.default_rng(1)

        analyzer = mopeds.tools.ErrorAnalyzer(vl_original, model, prediction_grid, meas_grid, selected_parameters, None, rng=rng)

        # fmt: off
        expected_df = np.array([[0.51852532, 2.46688906, 3.12389498, 3.32869702], [0.56753121, 2.66689296, 3.01510804, 3.7935085 ], [0.29663089, 2.430173  , 3.15858234, 4.12330323], [1.08274401, 1.40411709, 2.93449814, 4.7935197 ], [1.17740749, 1.7414517 , 3.05451991, 3.68074771], [1.75871993, 2.71211575, 2.80390916, 3.00223604], [1.11758469, 1.01196753, 2.90063113, 5.49061778], [2.09875925, 2.75065818, 2.8117853 , 2.7443598 ], [1.11590142, 1.23649892, 2.94801049, 5.00849873], [1.51800603, 1.38736098, 2.97966344, 4.31926724]])
        expected_predictions = np.array([[[ 9.89021025,  0.        ,  9.37168493], [11.95582903,  1.23344453,  9.37168493], [15.68579633,  2.46688905,  9.37168493]], [[ 9.61285534,  0.        ,  9.04532413], [11.89467895,  1.33344648,  9.04532413], [16.07325681,  2.66689296,  9.04532413]], [[ 9.7723779 ,  0.        ,  9.47574701], [12.01829021,  1.2150865 ,  9.47574701], [16.32585413,  2.430173  ,  9.47574701]], [[ 9.88623842,  0.        ,  8.80349441], [11.78667689,  0.70205855,  8.80349441], [16.08387521,  1.40411709,  8.80349441]], [[10.34096722,  0.        ,  9.16355973], [12.13187999,  0.87072585,  9.16355973], [15.76316662,  1.7414517 ,  9.16355973]], [[10.17044739,  0.        ,  8.41172747], [12.27706428,  1.35605788,  8.41172747], [15.88479919,  2.71211575,  8.41172747]], [[ 9.8194781 ,  0.        ,  8.7018934 ], [11.69811631,  0.50598377,  8.7018934 ], [16.32206341,  1.01196753,  8.7018934 ]], [[10.53411516,  0.        ,  8.43535591], [12.5955342 ,  1.37532909,  8.43535591], [16.02913313,  2.75065818,  8.43535591]], [[ 9.9599329 ,  0.        ,  8.84403148], [11.83030705,  0.61824946,  8.84403148], [16.20493055,  1.23649892,  8.84403148]], [[10.45699635,  0.        ,  8.93899033], [12.23049365,  0.69368049,  8.93899033], [16.16362458,  1.38736097,  8.93899033]]])
        # fmt: on

        analyzer.parameter_covariance_mc(plot=False, num_samples=10)
        assert np.isclose(expected_df, analyzer.df_params).all()

        analyzer.model_prediction_error_mc(plot=False)
        assert np.isclose(expected_predictions, np.array(analyzer.list_predictions)).all()


if __name__ == "__main__":
    pass
    # test_artificial_data_generator_nle()
    # test_startpoint_generation()
    test_error_analyzer(True)
