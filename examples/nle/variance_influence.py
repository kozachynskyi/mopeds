import mopeds
import numpy as np
import matplotlib.pyplot as plt


def get_model():
    variable_list = mopeds.VariableList()
    variable_list.add_variable(mopeds.VariableAlgebraic("y", 10.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("x", 1, 20))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("a", 10.2))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("b", 2.0))  # noqa: E501

    m = mopeds.Model(variable_list)

    x = m.varlist_all["x"].casadi_var  # noqa: E501
    y = m.varlist_all["y"].casadi_var  # noqa: E501
    a = m.varlist_all["a"].casadi_var  # noqa: E501
    b = m.varlist_all["b"].casadi_var  # noqa: E501

    EQ_alg1 = y - (x**b + x * a)

    list_algebraic_equations = [EQ_alg1]  # noqa: E501

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


if __name__ == "__main__":
    fig, axes = plt.subplots(2, 2)

    for i in [True, False]:
        for j in [True, False]:
            REAL_ERROR_ABSOLUTE = i
            ASSUMPTION_ABSOLUTE = j
            RELATIVE_ERROR = 0.15
            STD = 200 * RELATIVE_ERROR

            vl, m = get_model()
            bounds = {"x": [1, 30, 10]}
            exp_data, true_params, _ = (
                mopeds.tools.generate_artificial_data_from_grid_nle(
                    m, vl, bounds, perturbate=False
                )
            )

            x_values = []
            y_std_real = []
            y_std_assumed = []
            rng = np.random.default_rng(3)
            for vl_i in exp_data:
                std_rel = (vl_i["y"].value[0] * RELATIVE_ERROR) / 2

                if REAL_ERROR_ABSOLUTE:
                    vl_i["y"].value = rng.normal(vl_i["y"].value[0], STD)
                    y_std_real.append(STD)
                else:
                    vl_i["y"].value = rng.normal(vl_i["y"].value[0], std_rel)
                    y_std_real.append(std_rel)

                if ASSUMPTION_ABSOLUTE:
                    vl_i["y"].variance = (STD) ** 2
                    y_std_assumed.append(STD)
                else:
                    vl_i["y"].variance = (std_rel) ** 2
                    y_std_assumed.append(std_rel)

                vl_i["a"].fixed = False
                vl_i["b"].fixed = False
                x_values.append(vl_i["x"].value[0])

            pe = mopeds.ParameterEstimationNLE(m, exp_data)
            res = pe.optimize()
            v = pe.calculate_objective_and_residual(res["x_dict"])

            ax = axes[int(i), int(j)]
            ax.errorbar(
                v["y"].flatten(),
                pe.array_data.flatten(),
                yerr=y_std_real,
                ls="",
                marker="o",
            )
            ax.errorbar(
                v["y"].flatten(), pe.array_data.flatten(), xerr=y_std_assumed, ls=""
            )

            if ASSUMPTION_ABSOLUTE:
                ax.plot(pe.array_data, pe.array_data - STD)
                ax.plot(pe.array_data, pe.array_data + STD)
            else:
                ax.plot(pe.array_data, pe.array_data * (1 - RELATIVE_ERROR))
                ax.plot(pe.array_data, pe.array_data * (1 + RELATIVE_ERROR))

            if REAL_ERROR_ABSOLUTE:
                title = "Data has ABSOLUTE error\n"
            else:
                title = "Data has RELATIVE error\n"

            if ASSUMPTION_ABSOLUTE:
                title += "You assumed ABSOLUTE error\n"
            else:
                title += "You assumed RELATIVE error\n"

            ax.set_title(
                title + f"Real params {true_params}\n Estimated {res['x_dict']}"
            )

            ax.plot(pe.array_data, pe.array_data, ls="dashed")
    plt.show()
    breakpoint()
