# par_est

A library wrapped around casadi to solve Simulation / Optimization problems based on ODE and DAE models. Currently Parameter Estimation (PE) and Optimal Experimental Design (OED) are supported out of the box.

For a tutorial check https://www.user.tu-berlin.de/vovakozach/par_est/5min_tutorial.html

## Installation

- python3.8 or 3.9 is required
- install via pip or poetry
    - INSTALL: pip install git+https://git.tu-berlin.de/vovakozach/pe_oed_casadi
        - ITWM Repo: pip install git+https://gitlab.itwm.fraunhofer.de/bubel/parameter-estimation
    - UPDATE: pip install -U git+https://git.tu-berlin.de/vovakozach/pe_oed_casadi

### Installation ITWM

`pip` Installation:
```
pip install par_est --extra-index-url https://token:ceBtskLYpNqS6LR17NzP@gitlab.itwm.fraunhofer.de/api/v4/projects/3167/packages/pypi/simple
```

`poetry` Installation:
```
[[tool.poetry.source]]
name="par_est"
url="https://token:ceBtskLYpNqS6LR17NzP@gitlab.itwm.fraunhofer.de/api/v4/projects/3167/packages/pypi/simple"
```


## Development

- Clone this repo on your computer git clone https://git.tu-berlin.de/vovakozach/pe_oed_casadi
- Run `poetry install` (ensure that correct python version is installed ex. pyenv)
- Run tests via `pytest`, final tests should be run with `tox -r` command

### Old way

- You still can use `pip install -e .` in order to install package in development mode, but it's not recommended and will be deprecated

## TODO

Test_scaling in test_optimization shows that derivatives for first step with and without scaling are different, but it depends on a length of steps. Independent of ODE or DAE
Add tests for different integrators: ["idas", "collocation"]
Test and analyze hammersley generation in tools

Test_pendulum_dae fails if time_grid has more timesteps. Presicion of DAE solver should be checked
State variable starting value is redundant. It should be taken/stored from Experimental Data.
Control variable can hold only one control value for one experiment. No control change in time and additional experiments can be added.
Determine how to store time_grid - list or numpy array
in Optimizer if guess is 0 it's set to 1, which shouldn't hold all the time. Consider to determine scalingn based on lb and ub

Parameter estimation sequantially solves simulations. Multiprocessing doesn't work, because of picklign problem. Multithreading can be tested, but not sure it works. I shouldn't concentrate on that.

OED has builtin variances, it should be taken from state variables

## Development tips

- To mark something as unfinished use WIP (work in progress)
- To Debug casadi use: 
    - "print_in": True
    - "print_out": True
    - "verbose": True
    - "print_stats": True
- TO get list of function options: integrator.print_options()

## Structure

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

## Contributing

- Use poetry to work with package and contribute to it https://python-poetry.org/
- Use conventional commits to create commit messages https://www.conventionalcommits.org
    - git config --local commit.template .git-commit-message

## Note

This project has been set up using PyScaffold 3.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
