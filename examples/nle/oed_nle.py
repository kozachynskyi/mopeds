import par_est


def isomerization():
    VAR_LIST, MODEL, EXP_DATA = par_est.examples.isomerization_model()
    VAR_LIST["x1"].fixed = False
    VAR_LIST["x1"].lower_bound = 0.1
    VAR_LIST["x1"].upper_bound = 3
    VAR_LIST["theta1"].fixed = False
    oed = par_est.OptimalExperimentalDesign_NLE(MODEL, [VAR_LIST])
    oed.optimize()
    print(oed.calculate_objective_and_jacobian({"x1": 1}))
    breakpoint()

    # Example isomerization 1 Bates Page 56 Table 2.2
    res_x = pe.optimize()["x_dict"]
    print("expected 35.92, 0.0708, 0.0377, 0.167")
    print("Estimated par:\n", res_x)

    # Plot of convidence 95%-region page 58, Figure 2.18
    # They are not the same
    print("expected std: 8.21, 0.1783, 0.0988, 0.415")
    res = pe.parameter_analysis(res_x)


if __name__ == "__main__":
    isomerization()
