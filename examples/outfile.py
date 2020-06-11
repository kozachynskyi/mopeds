import copy
from datetime import datetime, timedelta

import casadi as ca
import matplotlib.cm as cm
import numpy as np
from matplotlib import pyplot as plt

import par_est


def initialize_problem():

    variable_list = par_est.VariableList()

    # fmt: off 
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_par4", 3.8E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_par5", -65610.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_par6", -26.64))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_x_i1_j1", 0.2))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i2", -1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i2", -3.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_x_i2_j1", 0.6))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i2", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_x_i3_j1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i1", 3.7E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i1", -0.0169))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i7", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i7", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i1", 5.45E-12))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i1", -2.44E-8))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i4", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i4", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i4", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i5", -1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i5", -1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i5", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i6", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i6", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i3", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r2_i3", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i3", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i3", -2.50993E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i3", 2.2E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_par2", 3.707))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_par3", -0.002783))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_par1", 4019.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i3", -241830.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i4", -1.54836E-11))
    variable_list.add_variable(par_est.Parameter_variable("e0_R", 8.314))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i3", -0.074))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i3", 42.04061276))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i4", -0.00322))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i4", 29.59614774))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i4", 2.10847E-8))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i4", 6.07E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_x_i4_j1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_x_i5_j1", 0.2))
    variable_list.add_variable(par_est.Parameter_variable("e0_x_i6_j1", 1.0E-16))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i1", 31.5))
    variable_list.add_variable(par_est.Parameter_variable("e0_T_j2", 500.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i3", 1.07368E-10))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i2", -7.52E-12))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i2", 2.7E-8))
    variable_list.add_variable(par_est.Parameter_variable("e0_T_f", 298.15))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i2", 26.2))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i2", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i2", -3.17E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i2", 0.0162))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_eq_r1_par1", -2073.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_eq_r1_par2", 2.029))
    variable_list.add_variable(par_est.Parameter_variable("e0_V", 0.1))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i6", -200940.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i7", 3.87561E-10))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_eq_r2_par1", 3066.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i6", -0.121))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_eq_r2_par2", -10.592))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i6", 45.78156724))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i7", -0.0526))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i7", 45.14491836))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i7", -8.56412E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i7", 6.29E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i7", -184100.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_p_j2", 50.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_r1", 1.07))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i5", 1.93148E-11))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i5", -1.32293E-8))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i4", -110530.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_i5", 19.6353642))
    variable_list.add_variable(par_est.Parameter_variable("e0_h_f_i5", -393520.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i5", -4.18E-5))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_i5", 0.0718))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_i6", 5.57E-4))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par1_i6", 2.65055E-10))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_i6", -6.5285E-7))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r3_i3", 92000.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_i6", 223.2))
    variable_list.add_variable(par_est.Parameter_variable("e0_F_j1", 1.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r3", -55060.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_r3_i6", 105100.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3", 85190.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_T_j1", 293.15))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r3_i6", -2.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_nu_r1_i7", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_rho_r1", 1775.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_par1_r1", 40000.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_greek_rho_r3", 100.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par3_r1", 0.499))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_par3_r1", 17197.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par2_r1", 3453.38))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_par2_r1", 0.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par5_r1", 1.22E10))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_par5_r1", -98084.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_par4_r1", 6.62E-11))
    variable_list.add_variable(par_est.Parameter_variable("e0_E_par4_r1", 124119.0))
    variable_list.add_variable(par_est.Parameter_variable("e0_A_r3_i3", 0.5498))

    m = par_est.Model(variable_list)


    # fmt: on"

    m.add_equations_differential([])
    m.add_equations_algebraic([])

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
    res.plot_states()
    # np.savetxt("exp.txt", res.toarray().T, delimiter="	")