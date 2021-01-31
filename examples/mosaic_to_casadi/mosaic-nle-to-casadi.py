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
from par_est.simulation import SimulatorNLE


def initialize_problem():

    variable_list = par_est.VariableList()

    # fmt: off """
)

counter = 0

for line in Lines:
    fields = line.strip().split()
    logging.debug(fields)
    if "X_ITER(" in line and "%" in line:
        output.append(
            f'    variable_list.add_variable(par_est.VariableState("{fields[4]}", {fields[2][:-1]}))'
        )
    elif "PARAMS(" in line and "%" in line:
        output.append(
            f'    variable_list.add_variable(par_est.VariableParameter("{fields[4]}", {fields[2][:-1]}))'
        )

output.append("\n    m = par_est.Model(variable_list)\n")

diff_equations = ""
counter = 0

# for line in Lines:
#     fields = line.strip().split()
#     if "Y(" in line:
#         logging.debug(fields)
#         counter += 1
#         new_line = f"    dydx{counter} = "
#         diff_equations += f"dydx{counter} ,"


#         for field in fields[2:]:
#             while "power" in field:
#                 field = re.sub(r"(power)(\()(\(.*?\)),","\g<3> ** (",field)

#             field = re.sub(r"e0_\w*", r"m.varlist_all['\g<0>'].casadi_var", field)
#             field = re.sub("exp", "ca.exp", field)
#             field = re.sub("log", "ca.log", field)
#             field = re.sub("\^", "\*\*", field)
#             field = re.sub(";", "", field)
#             # field = field.replace("exp", "ca.exp")
#             # field = field.replace("^", "**")
#             # field = field.replace(";", "")

#             new_line += f" {field} "
#         output.append(new_line)

for line in Lines:
    if "Y(" in line:
        logging.debug(line)
        counter += 1
        new_line = f"    dydx{counter} = "
        diff_equations += f"dydx{counter} ,"

        line = re.sub("Y\(\d*\) =", "", line)
        while "power" in line:
            line = re.sub(r"(power)(\()(\(.*?\)),", "\g<3> ** (", line)
        line = re.sub(r"e0_\w*", r"m.varlist_all['\g<0>'].casadi_var", line)
        line = re.sub("exp", "ca.exp", line)
        line = re.sub("log", "ca.log", line)
        line = re.sub("\^", "\*\*", line)
        line = re.sub(";", "", line)
        # field = field.replace("exp", "ca.exp")
        # field = field.replace("^", "**")
        # field = field.replace(";", "")

        new_line += line

        output.append(new_line)

output.append(
    f"""
    # fmt: on"

    m.add_equations_differential([{diff_equations}])

    return variable_list, m

if __name__ == "__main__":

    variable_list, m = initialize_problem()

    # Set parameters and controls to fixed state so their values are used for simulation
    var_list_fixed = copy.deepcopy(variable_list)
    for var in var_list_fixed.values():
        var.fixed = True

    # Create simulation Object
    sim_fixed = SimulatorNLE(m, var_list_fixed)
    # Run simulation and get simple results as array of numbers, but information about state variables and timestamp is lost
    res_simple = sim_fixed.simulate_sym()
    # Run simulation and connect results with actual state variables, which can be plotted based on available data
    res = sim_fixed.generate_exp_data()
    # res.plot_states()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="\t")"""
)

# print("\n".join(output))

with open("outfile.py", "w") as outfile:
    outfile.write("\n".join(output))
