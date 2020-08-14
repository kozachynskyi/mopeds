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

import import_data_nle

variable_list, model, experiments_list_together = import_data_nle.import_model_and_exp_data()

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
