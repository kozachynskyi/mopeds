import sys
import re
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format="%(message)s")

if len(sys.argv) == 2:
    path_to_file = sys.argv[1]
else:
    print("please give a  path")
    exit()

# Using readlines()
file1 = open(path_to_file, "r")
Lines = file1.readlines()

output = []

output.append(
    """import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est


def initialize_problem():

    variable_list = par_est.VariableList()

    # fmt: off"""
)

counter = 0

for line in Lines:
    fields = line.strip().split()
    if "MASS(" in line:
        state_var = int(line.split("eye(", 1)[1][0])

for line in Lines:
    fields = line.strip().split()
    if "Y_INIT(" in line:
        counter += 1
        if counter > state_var:
            output.append(
                f'    variable_list.add_variable(par_est.VariableAlgebraic("{fields[4]}", {fields[2][:-1]}))'
            )
        else:
            output.append(
                f'    variable_list.add_variable(par_est.VariableState("{fields[4]}", {fields[2][:-1]}))'
            )
    elif "PARAMS(" in line and "%" in line:
        output.append(
            f'    variable_list.add_variable(par_est.VariableParameter("{fields[4]}", {fields[2][:-1]}))'
        )

output.append("\n    m = par_est.Model(variable_list)\n")

counter = 0
diff_equations = ""
alg_equations = ""

for line in Lines:
    if "DYDX(" in line:
        counter += 1
        new_line = f"    dydx{counter} = "
        if counter > state_var:
            alg_equations += f"dydx{counter} ,"
        else:
            diff_equations += f"dydx{counter} ,"

        line.strip()
        line = re.sub("DYDX\(\d*\) =", "", line)
        while "power" in line:
            line = re.sub(r"(power)(\()(\(.*?\)),", "\g<3> ** (", line)
        line = re.sub(r"e0_\w*", r"m.varlist_all['\g<0>'].casadi_var", line)
        line = re.sub("exp", "ca.exp", line)
        line = re.sub("log", "ca.log", line)
        line = re.sub(r"\^", r"**", line)
        line = re.sub(";", "", line)

        new_line += line
        output.append(new_line)

output.append(
    f"""
    # fmt: on

    m.add_equations_differential([{diff_equations}])
    m.add_equations_algebraic([{alg_equations}])

    return variable_list, m

if __name__ == "__main__":

    variable_list, m = initialize_problem()
    # Create time-grid. Zero should be first
    time_grid = np.linspace(10, 10000, 40)
    time_grid = np.insert(time_grid, 0, 0)

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = par_est.Simulator(m, time_grid, var_list_fixed)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data()
    # res.plot_states()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="\t")"""
)

# print("\n".join(output))

with open("outfile.py", "w") as outfile:
    outfile.write("\n".join(output))
