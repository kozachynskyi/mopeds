import numpy as np

import par_est


def main():
    varlist, model, data = par_est.examples.puromycin_model()

    print(par_est.SimulatorNLE(model, varlist).simulate())

    pe = par_est.ParameterEstimationNLE(model, data["Treated"])

    # Puromycin 6 Bates Page 51
    res_x = pe.optimize(objective_function="ols")["x_dict"]
    print("res_x: ", res_x)
    print("s2 = ", pe.calculate_objective_and_residual(res_x)["f"] / 10)

    # Puromycin 6 Bates Page 53
    res_sens = pe.calculate_sensitivity_and_fim(res_x)
    print("Rinv =\n", np.linalg.inv(np.linalg.qr(res_sens["jac_full"])[1]))
    breakpoint()


if __name__ == "__main__":
    main()
