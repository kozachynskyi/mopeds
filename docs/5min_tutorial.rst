#################
5 Minute Tutorial
#################

Important Concepts
==================

Let's take a modified *BOD* example from "Nonlinear Regression Analysis and Its Applications" Bates, Watts 1988. NLE model with one function:

.. math::

    y = f(x,\theta) = \theta_1 (C - \exp(-\theta_2 x))

from point of view of |par_est|:

- :math:`y` - is an Algebraic Variable
- :math:`\theta_1` and :math:`\theta_2` - are Paramater Variables
- :math:`x` - is a Control Variable
- :math:`C` - is a Constant Variable

Let's create these variables::

    import par_est

    y = par_est.VariableAlgebraic("y", 8.3)
    x = par_est.VariableControl("x", 1)
    C = par_est.VariableConstant("C", 1)
    theta1 = par_est.VariableParameter("theta1", 20)
    theta2 = par_est.VariableParameter("theta2", 0.24)

Values that are used while creating variables have different meenings:

- Algebraic Variable is initiated with a *guess* that solver will use to solve an |NLE|
- Parameter and Control Variables will save provided value as a *guess* for optimizer and as a value for |NLE| solver
- Constant Variable will be subsituted by a provided value during the generation of the equation system

All variables have to be a part of a dictionary, instance of a VariableList. This class is used throughout the package and is usefull for plotting and result generation::

    variable_list = par_est.VariableList()
    variable_list.add_variable(y)
    variable_list.add_variable(x)
    variable_list.add_variable(C)
    variable_list.add_variable(theta1)
    variable_list.add_variable(theta2)

List of variables is further used to create a model instance, that holds infromation only about types of variables and how are they combined in equations::

    model = par_est.Model(variable_list)

Now we need to create an equation, using symbolic variables, that are stored for each Variable in attribute `.casadi_var`. Use |casadi| syntax  for complex functions, like `syn` or `log`. When creating equations be aware, that one should use symbolic variables from `model.varlist_all`::

    y = model.varlist_all["y"].casadi_var
    x = model.varlist_all["x"].casadi_var
    C = model.varlist_all["C"].casadi_var
    theta1 = model.varlist_all["theta1"].casadi_var
    theta2 = model.varlist_all["theta2"].casadi_var

    import casadi as ca
    equation = y - (theta1 * (1 - ca.exp(-theta2 * x)))

    model.add_equations_algebraic([equation])

Last line adds the equation to a model. From now on, one can "Simulate" or solve the equation system, because DoF is closed: all the variables are specified, and a guess for Algebraic Variable is provided::

    simulator = par_est.SimulatorNLE(model, variable_list)
    result = simulator.generate_exp_data()

Now let's look at a result, which is a VariableList with Algebraic Variables, that were calculated. Resulting value can be acessed via a `.dataframe` of the Variable of directly via `.value` property::

    print(result["y"].value)
    >>> [4.267442778668931]
    print(result["y"].dataframe)
    >>>                   y
    >>>1970-01-01  4.267443

In order to solve the example at different :math:`x`, one need to create a new Simulator, with the same `model`, but with different VariableList::

    variable_list["x"].value = 20
    simulator = par_est.SimulatorNLE(model, variable_list)
    result = simulator.generate_exp_data()

    print(result["y"].value)
    # >>> [19.8354050590196]

.. warning::
    You cannot change a value of a Constant Variable after a model is created!

Optimization
============

Now let's create a Parameter Estimation problem. First, one needs to provide experimental data. Let's say, that at :math:`x = 1` and :math:`x = 7`, :math:`y = 8.3` and :math:`y = 19.8` respectively. In order to infrom optimzer about this data sets we need to create a separate VariableList for each experimental point and assign values::

    import copy

    var_list_1 = copy.deepcopy(variable_list)
    var_list_1["x"].value = 1
    var_list_1["y"].value = 8.3

    var_list_2 = copy.deepcopy(variable_list)
    var_list_2["x"].value = 7
    var_list_2["y"].value = 19.8

    experimental_data = [var_list_1, var_list_2]

One needs to use `copy.deepcopy`. Generated experimental points can be further used for optimization. We need to specify, which Parameter to be optimized, by setting `.fixed` Flag, and providing guess and bound for optimizer. Result of an optimizer is a dictionary, where key `[x]` holds results::

    for var_list in experimental_data:
        var_list["theta1"].fixed = False
        var_list["theta1"].guess = 40
        var_list["theta1"].lower_bound = 0
        var_list["theta1"].upper_bound = 40

    pe = par_est.ParameterEstimationNLE(model, experimental_data)
    result = pe.optimize()
    print(result["x"])
    # >>> 25.2727
