# mopeds

<img align="right" src="https://git.tu-berlin.de/dbta/optimization/mopeds/-/raw/main/docs/logo.png" width="300px">

mopeds - **Mo**del based **P**arameter **E**stimation and **D**esign of Experiment**s** is a library wrapped around casadi to solve Simulation / Optimization problems based on steady state and dynamic models.

The system of equations can be a set of nonlinear equations (NLE), ordinary differential equations (ODE) or differential and algebraic equations (DAE) of index 1.

Parameter Estimation (PE) and Optimal Experimental Design (OED) are supported out of the box.

The project is looking for case studies and testers, so if you have any problems using the package or have any questions, do not hesitate to contact us.

Documentation is available [here](https://mopeds.readthedocs.io/en/latest/).


## Installation

`pip` Installation:

```
pip install mopeds
```

## Migration from par_est

If you used `par_est` before, in order to move to `mopeds` you need to replace the name, API did not change. In comparison to `par_est`, `mopeds` uses `casadi == 3.6.4`, so different results might be expected.
Considering creating a test to compare the numerical results while migrating.

## Citation

If you use MOPEDS in academic work, please cite:

Kozachynskyi, V., Illner, M., Esche, E., Repke, J.-U.
"The optimal experiment? Influence of solution strategies on model-based optimal experimental design"
Computers & Chemical Engineering, 2024.
https://doi.org/10.1016/j.compchemeng.2024.108746

## What's New?

### Upcoming

- Variables are automatically scaled, if lower and upper bound are provided, use_bounds for rootfinder NLE is deprecated
- Dynamic and steady state simulators have a consistent API now, sim.simulate() and sim.simulate_fast(), generate_exp_data is deprecated, used sim.simulate()[2] instead, to get a varlist.
- Model does not contain any variable lists, instead in just holds an order of variables and respecive casadi variables.
- Parameter Estimation and OED of NLE models supports direct optimization, before only the sequential optimization was used: for every experiment there was a rootfinder that found solution and provided gradient for optimizer.
- OED.optimize() ignore scale argument, use oed.objective_scaling instead
- VariableConstant allows multiple inputs, which are ignored. Used to easier switch from Independent Variable to constant
- Add tools.analyze_scaling() and pe.check_results_bounds to help with selection of scaling bounds
- Added linear example in mopeds.example
- Added "df_all" when calculating the objective and residual of the PE NLE
- Rework API of tools. Generate_varlist.. for NLE is now called generate_artificial_data..
- OED differently arranges jacobian, than before. Before it was sorted from top to bottom by measured variable, as in PE.jacobian. Now sorted simulation by simulation
- In PE, self.array_data has become a nlpsol_p parameter, meaning that PE data can be changed without reinitilization of the nlp solver, thus saving time
- Opimizer now support reusing the created .solver. It allows for fast repeated execution of the solver, e.g., for Monte Carlo simulations

### 0.10.3

- Added support and CI infrastructure for Python 3.13 and NumPy 2.x.

### 0.10.2

- Fix the numpy dependency to <2 and casadi<3.7 to avoid errors

### 0.10.1

- Fixed installation error in Windows python 3.11

### 0.10.0

- Rename par_est to mopeds and open-source the package

### 0.9.3.a1

- OED of dynamic models supports multiple different modes and strategies
- Added multiple regularization techniques

### 0.9.2

- fix WLS formulation  (remove division by 2)
- feature -> remove rounding of time_grid in Simulator

### 0.9.1

- fix bugs in identifiability analysis

### 0.9.0

- BREAKING: DAE simulators API change: from self.simulate() to self.simulate_sym()
- Added support for ACADOS ODE / DAE simulator
- Rework how PE for DAE and NLE works -> more simmilar code, easier to maintain
- Rework how Confidence Intervals of Parameters are calculated for multivariate measurements with different variance

### 0.8.0

NLE Simulator and Parameter Estimation were reworked, with focus on analysis of parameter variance-covariance matrix.
Parameter Estimation has different internals on how objective function is calculated, making it a bit faster and much more unrestandable.
Examples from Bates, Watts "Nonlinear Regression analysis and its applications" were imlemented and tested.

## Contributors

Many people have been involved in the development of this package, either by writing actual code, helping with the methods behind it, or simply using it and providing feedback and feature requests. Here are just a few names:

- Volodymyr Kozachynskyi
- Dario Staubach
- Martin Bubel
- Lorenz Hafner
- Mudassar Javed
- Torben Talis
- Joris Weigert
- Erik Esche
- Markus Illner
- Christian Hoffman
- Georg Brösigke
- Maria Stockman

and many, many others ...

## Development and Contributions

The primary repository of MOPEDS is hosted on TU Berlin GitLab:

https://git.tu-berlin.de/dbta/optimization/mopeds

A public GitHub mirror is available to improve visibility and simplify community contributions. Users are welcome to use either platform to:

- Report bugs
- Request features
- Ask questions
- Submit merge requests / pull requests

Changes submitted through GitHub are reviewed and synchronized with the primary GitLab repository.
Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## Acknowledgement

This work is funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - 56091768 and 466397921.
