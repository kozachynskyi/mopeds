import matplotlib.pyplot as plt

import par_est

VAR_LIST, MODEL, EXP_DATA = par_est.examples.puromycin_model()

data = EXP_DATA["Treated"]
list_of_params = list(["theta1", "theta2"])

dict_of_params = {}
for param in list_of_params:
    dict_of_params[param] = float(data[0][param].value[0])

dict_of_controls = {
    "x": [0.0, 1.2, 120 + 1],
}

dict_of_responses = {
    "f": 5e1,
}

pe = par_est.ParameterEstimationNLE(MODEL, data)
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
plt.xlim(0, 1.2)
plt.xlabel("Concentration")
plt.ylim(0, 250)
plt.ylabel("Velocity")
plt.grid()
plt.show()

dict_of_artificial_controls = {
    "x": [0.0, 1.2, 10 + 1],
}
(
    artificial_inference_results,
    artificial_data,
    artificial_sim_data,
) = pe.calculate_inference_bounds(
    dict_of_params,
    dict_of_responses,
    dict_of_controls,
    dict_of_artificial_controls,
    seed=42,
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
plt.xlim(0, 1.2)
plt.xlabel("Concentration")
plt.ylim(0, 250)
plt.ylabel("Velocity")
plt.grid()
plt.show()

print(exp_inference_results["f"]["s"])
print(artificial_inference_results["f"]["s"])
