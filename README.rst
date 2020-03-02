=======
par_est
=======

A project to play around with Parameter Estimation (PE) and Optimal Experimental Design (OED) for ODE and DAE system of equations.

Math is handled via casadi.

TODO
===========

State variable starting value is redundant. It should be taken/stored from Experimental Data.
Control variable can hold only one control value for one experiment. No control change in time and additional experiments can be added.
Determine how to store time_grid - list or numpy array
in Optimizer if guess is 0 it's set to 1, which shouldn't hold all the time. Consider to determine scalingn based on lb and ub

Parameter estimation sequantially solves simulations. That can be speed up by multiprocessing.

OED has builtin variances, it should be taken from state variables

Structure
========

Experimental Data is object that stores value and time arrays.

Variables can be State, Parameter and Control. They have name and casadi variable connected to them.
State Variable has starting value, value as experimental data object and opcua_id
Parameter variables have guess value, value itself, and bounds for optimizer
Control variable same as parameter but has opcua_id

Variable list is ordered dict with Variables as values, and variable names as keys
Variables added here via add_variable(). Opc write and read handled here. 
State variables can be plotted via plot_states()

Model holds equations and all variables without any values. It makes difference only between State and Variables (THey represented as ParameterVariable). It is needed to add consistently equations.

Simulator does integrates a given model. It's intiated with Model, time_grid and Variable list. It doesn't care if variable is Parameter or Control variable. Matters only a fixed state: if variable is fixed it's value is used, if not, it's considered a casadi var for simulator and optimizer.
_reset_scaling() is important method that handles scaling of simulation
simulate() runs integrator at every time_grid point, and can also provide derivatives
generate_exp_data() transforms output of simulator to VariableList. In such a way simulated data can be plotted and analyzed.
Simulator holds variable list given on initialization.

Optimizer is intialized based only on model and variable_list. It's a general class that holds input variable list
scaling is either one number or array of values for scaling
_setup_simulator is used for filling variable lists
_setup_initialization sets guess for desicion variables and their bounds
_setup_scaling is determining scaling for desicion variables based on current guess
_objective generates objective calculation for optimization problem
_optimize last step that launches optimization solver

parameter estimation determines time_grid for itself, based on available state_variable data. It also fixes all controll variables.

OED uses time grid provided to it. It uses unfixed contoll variables as desicion variables and unfixed parameters to calculate sensitivities
_sensitivity_matrix is used to calculated Jy/dp
get_fim() is used to calculate FIM for a given set of variables

Requirements
============

python => 3
casadi https://github.com/casadi/casadi/wiki/InstallationInstructions

Note
====

This project has been set up using PyScaffold 3.2.2. For details and usage
information on PyScaffold see https://pyscaffold.org/.
