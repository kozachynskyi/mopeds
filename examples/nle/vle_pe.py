import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est
from par_est.simulation import SimulatorNLE
from par_est.optimization import ParameterEstimationNLE

def generate_exp_data_list_NLE(model, var_list_fixed):
        # Create simulation Object
        sim_fixed = SimulatorNLE(model, var_list_fixed)
        # Run simulation and connect results with actual state variables
        val_fix = sim_fixed.generate_exp_data()
        # Copy variable_list
        variable_list_optimizer = copy.deepcopy(var_list_fixed)
        # Set startings values
        variable_list_optimizer.set_starting_values(val_fix, True)

        return variable_list_optimizer


def initialize_problem():
    # Id. VLE of EtOH and Water

    # Variables
    variable_list = par_est.variables.VariableList()  # Preallocate variable_list

    # Define variables
    #     T in K
    #     x in 1
    #     P in Pa
    #     # EtOH = 1,      H2O = 2
    #     a = [5.24125,    5.19625] # a in 1
    #     b = [1592.864,   1730.630]# b in K
    #     c = [-46.9659,   -39.7239] # c in K

    variable_list.add_variable(par_est.VariableState("T", 373))
    variable_list.add_variable(par_est.VariableControl("x", 0.5))
    variable_list.add_variable(par_est.VariableControl("P", 1e5))
    variable_list.add_variable(par_est.VariableParameter("a1", 5.24125))
    variable_list.add_variable(par_est.VariableParameter("a2", 5.19625))
    variable_list.add_variable(par_est.VariableParameter("b1", 1592.864))
    variable_list.add_variable(par_est.VariableParameter("b2", 1730.630))
    variable_list.add_variable(par_est.VariableParameter("c1", -46.9659))
    variable_list.add_variable(par_est.VariableParameter("c2", -39.7239))

    model = par_est.Model(variable_list)  # adding all variables to the model

    # Equations
    RES = model.varlist_all['P'].casadi_var- (model.varlist_all["x"].casadi_var*10**(model.varlist_all['a1'].casadi_var-model.varlist_all['b1'].casadi_var/(model.varlist_all['c1'].casadi_var+model.varlist_all["T"].casadi_var))*1E5 + (1-model.varlist_all["x"].casadi_var)*10**(model.varlist_all['a2'].casadi_var-model.varlist_all['b2'].casadi_var/(model.varlist_all['c2'].casadi_var+model.varlist_all["T"].casadi_var))*1E5)
    model.add_equations_differential([RES])  # adding the equations to model

    return variable_list, model


if __name__ == "__main__":

    variable_list, model = initialize_problem()

    param_list = [
        "a1",
        "b1",
        "c1",
        "a2",
        # "b2",
        # "c2",
    ]  # List of parameters for parameter estimation

    # Set parameters and controls to fixed state so their values are used for
    # simulation
    var_list_fixed = copy.deepcopy(variable_list)

    # Fix all variables
    var_list_fixed.set_variable_list_fixed()

    # Preallocate list for optimizer_lists
    optimizer_list = []

    # Define changing mole fraction x
    x = [0.00, 0.10, 0.25, 0.50, 0.75, 0.90, 1.00]


    for val in x:
        var_list_fixed['x'].value = val # adjust the value for x

        # Generate data using the fixed values for the parameter. The output
        # is used as starting values for parameter estimation
        variable_list_optimizer = generate_exp_data_list_NLE(model, var_list_fixed)

        # Unfix parameters for parameter estimation
        variable_list_optimizer.set_variable_list_unfixed(param_list)

        # Set upper and lower bounds
        variable_list_optimizer.set_bounds(emerg_val=50)

        # Append variable_list_optimizer to optimizer_list
        optimizer_list.append(variable_list_optimizer)


    # Define ParameterEstimationNLE
    pe = ParameterEstimationNLE(
        model, optimizer_list
    )

    # Solve parameter estimation problem
    res = pe.optimize(False)
