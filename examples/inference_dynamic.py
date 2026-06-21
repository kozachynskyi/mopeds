import mopeds
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # VAR_LIST, MODEL, EXP_DATA = mopeds.examples.puromycin_model()
    # VAR_LIST["f"].variance = 20.93  ** 2
    # controls_list =  [
    #     [0.02, 76],
    #     [0.02, 47],
    #     [0.06, 97],
    #     [0.06, 107],
    #     [0.11, 123],
    #     [0.11, 139],
    #     [0.22, 159],
    #     [0.22, 152],
    #     [0.56, 191],
    #     [0.56, 201],
    #     [1.10, 207],
    #     [1.10, 200],
    # ]

    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.bod_model()
    VAR_LIST["f"].variance = 0.2  #  ** 2
    controls_list = [
        [1, 8.3],
        [2, 10.3],
        [3, 19.0],
        [4, 16.0],
        [5, 15.6],
        [7, 19.8],
    ]

    cl = []
    for ci in controls_list:
        cl.append({"x": ci[0]})

    list_params = []

    exp_data, true_params, _ = mopeds.tools.generate_artificial_data_nle(
        MODEL, VAR_LIST, cl, perturbate=False
    )
    exp_data[0]["theta1"].fixed = False
    exp_data[0]["theta2"].fixed = False
    pe = mopeds.ParameterEstimationNLE(MODEL, exp_data)
    all_df = pd.concat([vl.dataframe for vl in exp_data])
    rng = np.random.default_rng(1)

    for i in range(10):
        perturb_data = rng.normal(all_df["f"], VAR_LIST["f"].variance ** 0.5)
        pe.array_data = perturb_data

        res = pe.optimize(direct_optimization=True, reuse_solver=False)
        list_params.append(res["x_dict"])

    list_prediction = []
    list_x = []

    grid_num = 20
    ub = 8

    x = np.linspace(0, ub, grid_num)
    for par in list_params:
        for par_name, par_value in par.items():
            VAR_LIST[par_name].value = par_value
            prediction, _, _ = mopeds.tools.generate_artificial_data_from_grid_nle(
                MODEL, VAR_LIST, {"x": [0.0, ub, grid_num]}, perturbate=False
            )
            y_data = []
            for vl in prediction:
                y_data.append(float(vl["f"].value[0]))
            # v = mopeds.ParameterEstimationNLE(MODEL, prediction)
            list_prediction.append(y_data)
            list_x.append(x)

    df = pd.DataFrame(np.squeeze(np.array(list_prediction)), columns=x)
    ax = df.mean().plot()
    # print(df.std())
    ax = (df.mean() + 2 * df.std()).plot(ax=ax)
    ax = (df.mean() - 2 * df.std()).plot(ax=ax)

    # df.T.plot()
    plt.show()
    # ax = df.min().plot(ax=ax)
    # ax = df.max().plot(ax=ax)
    # plt.show()
    breakpoint()

    breakpoint()
    list_of_params = list(["theta1", "theta2"])

    dict_of_params = {}
    for param in list_of_params:
        dict_of_params[param] = float(data[0][param].value[0])

    x_bounds = [0.0, 1.2]
    y_bounds = [0.0, 250]
