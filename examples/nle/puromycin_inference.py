import par_est
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from scipy.stats import f

VAR_LIST, MODEL, EXP_DATA = par_est.examples.puromycin_model()

data = EXP_DATA["Treated"]
time_exp = []
fun = []
fun_sim = []
for experiment in data:
    time_exp.append(experiment.dataframe["x"].to_numpy())
    fun.append(experiment.dataframe["f"].to_numpy())
   
concentration_grid = np.linspace(0, 1.2, 49, endpoint=True)  
simulation_var_lists = [] 
for concentration in concentration_grid:
    simulation_var_list = deepcopy(VAR_LIST)
    simulation_var_list['x'].value = concentration
    simulation_var_list['f'].value = 1
    for var in simulation_var_list.values():
        var.fixed = True
        
    res = par_est.SimulatorNLE(MODEL, simulation_var_list).simulate_sym()   
    simulation_var_list["theta1"].fixed = False
    simulation_var_list["theta2"].fixed = False
    simulation_var_lists.append(simulation_var_list)
    fun_sim.append(res["x"].toarray().flatten()) 
    
pe = par_est.ParameterEstimationNLE(MODEL, data)
pe_grid = par_est.ParameterEstimationNLE(MODEL, simulation_var_lists)
param_dict = {}
for param_name in list(pe.varlist_parameter.keys()
):  
    param_dict[param_name] = float(data[0].dataframe[param_name].to_numpy()[0])

OLS = pe.calculate_objective_and_residual(param_dict)["f"]
sens_info = pe.calculate_sensitivity_and_fim(param_dict,list(pe.varlist_parameter.keys()))
sens_info_grid = pe_grid.calculate_sensitivity_and_fim(param_dict,list(pe.varlist_parameter.keys()))
len_exp = len(data)
len_param = len(param_dict)
jac_grid = sens_info_grid["jac_full"]
jac_grid_sorted = sens_info_grid["jac_sorted"]["f"]
print(jac_grid_sorted)
DOF = len_exp-len_param
s = np.sqrt(OLS/DOF)
R = np.linalg.qr(sens_info["jac_full"], mode="reduced")[1]
fisher95 = f(len_param, DOF).ppf(0.95)
val = jac_grid@np.linalg.inv(R)
bound_array = s*np.linalg.norm(jac_grid@np.linalg.inv(R), axis=1)*np.sqrt(len_param*fisher95)
concentration_grid = np.hstack(concentration_grid)
fun_sim = np.hstack(fun_sim)
    
plt.plot(time_exp, fun, "ko", label="Exp. data")
plt.plot(concentration_grid, fun_sim, "k-", label="Sim. data")
plt.plot(concentration_grid, fun_sim-bound_array, "b:", label="Lower bound")
plt.plot(concentration_grid, fun_sim+bound_array, "r:", label="Upper bound")
plt.legend()
plt.xlim(0, 1.2)
plt.xlabel("Concentration")
plt.ylim(0, 250)
plt.ylabel("Velocity")
plt.grid()
plt.show()