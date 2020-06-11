import os
import csv
import par_est
import copy

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

from par_est.simulation import SimulatorNLE
from par_est.optimization import ParameterEstimationNLE

import outfileGB

# This allows for a more robust script execution
DATA_PATH = os.path.join(os.path.dirname(__file__), "experimental_data.csv")

results = []
with open(DATA_PATH) as csvfile:
    reader = csv.reader(csvfile, delimiter=",")  # change contents to floats
    for column in reader:  # each row is a list
        results.append(column)
    print(results)

# Why not to already use a varlist from your model?
var_list, model = outfileGB.initialize_problem()

# This part is hard to explain breefly. The way list intializtion is tricky in python
# you can get something that you do not expect, if you are not carefull
experiments_list_together = []
for _ in range(len(results[1:])):
    experiments_list_together.append(copy.deepcopy(var_list))

# Do not check for names, varlist already knows which variable type is specific variable name
for n, experiments in enumerate(experiments_list_together):
    for name in results[0]:
        if name in experiments.get_variable_name():
            current_variable = experiments[name]
            # When you read from csv, you get string, you need to use float()
            experiment_value = float(results[n + 1][results[0].index(name)])
            if isinstance(current_variable, par_est.VariableState):
                current_variable.value.value = experiment_value
            else:
                current_variable.value = experiment_value


print(experiments_list_together[0]["e0_T_j2"].value)


initialGuess = SimulatorNLE(model, experiments_list_together[0])
resInit = initialGuess.generate_exp_data()

for key, var in resInit.items():
    var.starting_value = var.value.value
    experiments_list_together[0][key] = var

for var in experiments_list_together[0].values():
    var.fixed = True
    if isinstance(var, par_est.VariableParameter):
        var.lower_bound = var.value - var.value*0.05
        var.upper_bound = var.value + var.value*0.05
        var.guess = var.lower_bound
experiments_list_together[0]["e0_A_r3_i3"].fixed = False 

# for experiments in experiments_list_together:
#     for key, var in resInit.items():
#         var.starting_value = var.value.value
#         experiments[key] = var

# for experiments in experiments_list_together:
#     for var in experiments.values():
#         var.fixed = True
#         if isinstance(var, par_est.VariableParameter):
#             var.lower_bound = var.value - var.value*0.05
#             var.upper_bound = var.value + var.value*0.05
#             var.guess = var.lower_bound

        
# for n, experiments in enumerate(experiments_list_together):
#     experiments_list_together[n]["e0_A_r3_i3"].fixed = False
   
    
peGB = ParameterEstimationNLE(model, [experiments_list_together[0]])

#res = peGB.optimize(False)