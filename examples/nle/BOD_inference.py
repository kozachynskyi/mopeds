import par_est
import matplotlib.pyplot as plt
from copy import deepcopy
import numpy as np
from scipy.stats import f

VAR_LIST, MODEL, EXP_DATA = par_est.examples.bod_model()

VAR_LIST["theta1"].value = 19.143
VAR_LIST["theta2"].value = 0.5311

time_exp = []
fun = []
fun_sim = []
for experiment in EXP_DATA:
    experiment["theta1"].value = 19.143
    experiment["theta2"].value = 0.5311
    time_exp.append(experiment.dataframe["x"].to_numpy())
    fun.append(experiment.dataframe["f"].to_numpy())

simulation_var_lists = []   
time_grid = np.linspace(0, 8, 33, endpoint=True)   
for time in time_grid:
    simulation_var_list = deepcopy(VAR_LIST)
    simulation_var_list['x'].value = time
    simulation_var_list['f'].value = 1
    for var in simulation_var_list.values():
        var.fixed = True
        
    res = par_est.SimulatorNLE(MODEL, simulation_var_list).simulate_sym()   
    fun_sim.append(res["x"].toarray().flatten()) 
    
    simulation_var_list["theta1"].fixed = False
    simulation_var_list["theta2"].fixed = False
    simulation_var_lists.append(simulation_var_list)
    
pe = par_est.ParameterEstimationNLE(MODEL, EXP_DATA)
pe_grid = par_est.ParameterEstimationNLE(MODEL, simulation_var_lists)
param_dict = {}
for param_name in list(pe.varlist_parameter.keys()
):  
    param_dict[param_name] = float(EXP_DATA[0].dataframe[param_name].to_numpy()[0])

OLS = pe.calculate_objective_and_residual(param_dict)["f"]
jac = pe.calculate_sensitivity_and_fim(param_dict,list(pe.varlist_parameter.keys()))["jac_sorted"]["f"]
jac_grid = pe_grid.calculate_sensitivity_and_fim(param_dict,list(pe.varlist_parameter.keys()))["jac_sorted"]["f"]
len_exp = len(EXP_DATA)
len_param = len(param_dict)
DOF = len_exp-len_param
s = np.sqrt(OLS/DOF)
R = np.linalg.qr(jac, mode="reduced")[1]
fisher95 = f(len_param, DOF).ppf(0.95)
bound_array = s*np.linalg.norm(jac_grid@np.linalg.inv(R), axis=1)*np.sqrt(len_param*fisher95)
time_grid = np.hstack(time_grid)
fun_sim = np.hstack(fun_sim)
    
plt.plot(time_exp, fun, "ko", label="Exp. data")
plt.plot(time_grid, fun_sim, "k-", label="Sim. data")
plt.plot(time_grid, fun_sim-bound_array, "b:", label="Lower bound")
plt.plot(time_grid, fun_sim+bound_array, "r:", label="Upper bound")
plt.legend()
plt.xlim(0, 8)
plt.xlabel("Time")
plt.ylim(-5, 30)
plt.ylabel("Oxygen demand")
plt.grid()
plt.show()