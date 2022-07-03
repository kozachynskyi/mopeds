#################
Available options
#################

NLE Simulator
=============

You create a simulator while calling::

    simulator = par_est.SimulatorNLE(model, variable_list)

This line will actually mean the following::

    solver_settings = {
        "nlpsol": "ipopt",
        "verbose": False,
        "print_in": False,
        "print_out": False,
        "expand": True,
        "nlpsol_options": {
            "ipopt.hessian_approximation": "limited-memory",
            "ipopt.max_iter": 300,
            "ipopt.print_level": 0,
            "print_time": False,
        },
    }

    simulator = par_est.SimulatorNLE(
        model,
        variable_list,
        solver_settings=solver_settings,
        solver_name="rootfinder",
        use_bounds=True,
    )

By default, simulator is using a `rootfinder` algorithm of the casadi, with `IPOPT` as a NLP solver. `use_bounds=True` means, that `constraints` option of the rootifnder will be set, based on the bounds of the Algebraic Variables that user has chosen. To read more about all available options check `official API <https://web.casadi.org/python-api/#rootfinding>`_.

.. warning::
    Use of `ipopt` as a solver_name of the SimulatorNLE is not fully supported

.. warning::
    Use of `fast_newton` and `newton` rootfinder options of CasADi is also not supported

IPOPT Settings
--------------

If IPOPT fails -> it's 75% chance that you cannot fix it by tweaking the settings, but rather by either:

- Changing bounds of the Algebraic Variables
- Setting a better guess for NLE, or even writing a routine that will calculate guess before rootfinder finds roots
- Choosing better bounds for the Optimizer

However there is a set of settings, that is more robust and slower, than the default ones. It requires a use of `MA57 <https://www.hsl.rl.ac.uk/catalogue/ma57.html>`_ solver from HSL Library::

    solver_settings = {
        "nlpsol": "ipopt",
        "verbose": False,
        "print_in": False,
        "print_out": False,
        "expand": True,
        "nlpsol_options": {
            "ipopt.hessian_approximation": "exact",
            "ipopt.linear_solver": "ma57",
            "ipopt.ma57_automatic_scaling": "yes",
            "ipopt.max_iter": 300,
            "ipopt.print_level": 0,
            "print_time": False,
        },
    }

If you need to supply other `IPOPT options <https://coin-or.github.io/Ipopt/OPTIONS.html>`_, just put them all in the `solver_settings` dictionary as shown above. But really, you don't need other options, your model is just bad =)

NLP Settings
------------

You can also supply more options from `here <https://web.casadi.org/python-api/#nlp>`_ under `nlp_sol` options, or even change the solver::

    solver_settings = {
        "nlpsol": "qrsqp",
        "verbose": False,
        "print_in": False,
        "print_out": False,
        "expand": True,
        "nlpsol_options": {
            "print_iteration": False,
        },
    }

    simulator = par_est.SimulatorNLE(
        model,
        variable_list,
        solver_settings=solver_settings,
        solver_name="rootfinder",
        use_bounds=True,
    )


Parameter Estimation NLE
========================

When you create a Parameter Estimation instance like this::
    
    pe = par_est.ParameterEstimationNLE(model, experimental_data)

you actually do the following::

    solver_settings = {
        "nlpsol": "ipopt",
        "verbose": False,
        "print_in": False,
        "print_out": False,
        "expand": True,
        "nlpsol_options": {
            "ipopt.hessian_approximation": "limited-memory",
            "ipopt.max_iter": 300,
            "ipopt.print_level": 0,
            "print_time": False,
        },
    }

    pe = par_est.ParameterEstimationNLE(
        model,
        experimental_data,
        simulator_settings=solver_settings,
        simulator_name="rootfinder",
        use_simulator_bounds=True,
    )

    pe.solver_name = "ipopt"
    pe.solver_settings = {
        "verbose": False,
        "ipopt": {"max_iter": 300},
    }

Optimizer creates internal Simulators with default settings, which you can change, as described in :ref:`NLE Simulator`. If you want to change Optimizer settings, you can do that by setting it's attributes on the fly::

    pe = par_est.ParameterEstimationNLE(model, experimental_data)
    pe.solver_name = "qrsqp"
    pe.solver_settings = {
        "verbose": False,
        "max_iter": 300,
    }

Basically all the settings discussed above can be applied here, with only difference, that they are changed on the fly, after the object is created.

You should use HSL Library solvers in order to be sure in your optimization, but in general it should work. If you see `inf` or `NaN` errors during the optimization, and IPOPT starts to change the `alpha`, that probably means that Simulators, that are being solved with current step of the optimizer cannot converge - meaning you need to tweak access the converge of your Simulator at different inputs.
