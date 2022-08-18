import matplotlib.pyplot as plt

import par_est

VAR_LIST, MODEL, EXP_DATA = par_est.examples.bod_model()

dict_of_params = {
    "theta1": 19.143,
    "theta2": 0.5311,
}

dict_of_controls = {
    "x": [0.0, 8.0, 90 + 1],
}

dict_of_responses = {
    "f": 1e1,
}

pe = par_est.ParameterEstimationNLE(MODEL, EXP_DATA)
exp_inference_results, exp_data, sim_data = pe.calculate_inference_bounds(
    dict_of_params, dict_of_responses, dict_of_controls
)

plt.plot(exp_data["x"], exp_data["f"], "ko", label="Exp. data")
plt.plot(
    exp_inference_results["x"],
    exp_inference_results["f"]["simulation"],
    "k-",
    label="Sim. data",
)
plt.plot(
    exp_inference_results["x"],
    exp_inference_results["f"]["lower bound"],
    "b:",
    label="Lower bound",
)
plt.plot(
    exp_inference_results["x"],
    exp_inference_results["f"]["upper bound"],
    "r:",
    label="Upper bound",
)
plt.legend()
plt.xlim(0, 8)
plt.xlabel("Time")
plt.ylim(-5, 30)
plt.ylabel("Oxygen demand")
plt.grid()
plt.show()

dict_of_artificial_controls = {
    "x": [0.0, 8.0, 10 + 1],
}
(
    artificial_inference_results,
    artificial_data,
    artificial_data_sim,
) = pe.calculate_inference_bounds(
    dict_of_params, dict_of_responses, dict_of_controls, dict_of_artificial_controls
)

plt.plot(artificial_data["x"], artificial_data["f"], "ko", label="Artifial data")
plt.plot(
    artificial_inference_results["x"],
    artificial_inference_results["f"]["simulation"],
    "k-",
    label="Sim. data",
)
plt.plot(
    artificial_inference_results["x"],
    artificial_inference_results["f"]["lower bound"],
    "b:",
    label="Lower bound",
)
plt.plot(
    artificial_inference_results["x"],
    artificial_inference_results["f"]["upper bound"],
    "r:",
    label="Upper bound",
)
plt.legend()
plt.xlim(0, 8)
plt.xlabel("Time")
plt.ylim(-5, 30)
plt.ylabel("Oxygen demand")
plt.grid()
plt.show()

print(exp_inference_results["f"]["s"])
print(artificial_inference_results["f"]["s"])
