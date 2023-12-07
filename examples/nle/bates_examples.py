import numpy as np

import mopeds


def isomerization():
    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.isomerization_model()
    pe = mopeds.ParameterEstimationNLE(MODEL, EXP_DATA)

    # Example isomerization 1 Bates Page 56 Table 2.2
    res_x = pe.optimize()["x_dict"]
    print("expected 35.92, 0.0708, 0.0377, 0.167")
    print("Estimated par:\n", res_x)

    # Plot of convidence 95%-region page 58, Figure 2.18
    # They are not the same
    print("expected std: 8.21, 0.1783, 0.0988, 0.415")
    res = pe.parameter_analysis(res_x)


def bod():
    VAR_LIST, MODEL, EXP_DATA = mopeds.examples.bod_model()
    pe = mopeds.ParameterEstimationNLE(MODEL, EXP_DATA)

    # Example BOD 5 Bates Page 54
    res = pe.optimize()["x_dict"]
    print("expected par 19.143, 0.5311")
    print("Estimated par:\n", res)

    print("S2 expected 6.498")
    print("S^2: ", pe.calculate_objective_and_residual(res, "ols")["f"] / 4)

    # Plot of convidence 95%-region page 55
    res = pe.parameter_analysis(res)


def puromycin():
    varlist, model, data = mopeds.examples.puromycin_model()
    pe = mopeds.ParameterEstimationNLE(model, data["Treated"])

    # Puromycin 6 Bates Page 51
    res_x = pe.optimize(objective_function="ols")["x_dict"]
    print("expected = {'theta1': 221.7, 'theta2': 0.0641}")
    print("res_x: ", res_x)
    print("expected 119.5")
    print("s2 = ", pe.calculate_objective_and_residual(res_x)["f"] / 10)

    # Puromycin 7 Bates Page 53
    print("expected standard error 6.95 and 8.28e-3")
    res_sens = pe.calculate_sensitivity_and_fim(res_x)
    print("Rinv =\n", np.linalg.inv(np.linalg.qr(res_sens["jac_full"])[1]))

    res = pe.parameter_analysis(res_x)


def spmma():
    # Doesn't work as in book
    varlist, model, data = mopeds.examples.spmma()
    pe = mopeds.ParameterEstimationNLE(model, data)

    # Puromycin 6 Bates Page 51
    res_x = pe.optimize(objective_function="ols")["x_dict"]
    a = pe.calculate_objective_and_residual(res_x)


# bod()
# spmma()
# puromycin()
isomerization()
