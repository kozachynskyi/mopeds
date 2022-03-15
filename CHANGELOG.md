## 0.5.2 (2022-03-15)

### Fix

- add importlib_metadata for python<3.8

### Feat

- remove setup.py
- add __version__
- added berty NLE example

## 0.5.1 (2022-03-09)

### Feat

- support for python 3.7 till 3.10
- NLE support unfixed variables in generate_exp_data

### Fix

- OptimizerNLE raise error when scale is used
- NLE optimizer ignores nan values in _objective()

## 0.5.0 (2022-02-19)

### Feat

- add NLE examples
- add NLE example to src
- add Parameter Estimation NLE
- add SImulatorNLE to __init__
- add varlist generator for NLE
- add Simulator NLE

### Fix

- remove set_starting_value, too coomplicated

### Refactor

- remove breakpoint

## 0.4.0 (2022-01-05)

### Feat

- set big markers only for minimal values
- added plotly backend option
- generate_exp_data() finds t0 for algebraic
- added _simulate_t0() method
- added plot_simulation() method for optimizers
- plot() methods support axis argument
- add map_objective method to optimzer
- add recalculate_algebraic argument to PE
- optimizer always recalculate algebraic vars
- remove "calculate_algebraic_experimental"
- added recalculate_algebraic option to simulator
- added integrator with algebraic recalculation
- analyze jacobian via PandasGUI, experimental
- added check_decision_bounds method
- generate_exp_data supports unfixed variables
- add ignore_plotting flag to variables
- added pandasqui
- added pandas, took changes from previous work
- calculation initial algebraic val with idas
- Make .factory call in Simulator optional
- MPC ca now use algebraic and state vars
- PE can calculate error of algebraic vars
- added singe shooting MPC
- added singe shooting MPC
- Remove pip install -e support, update versions

### Refactor

- added plot argument to plot_simulation
- used black
- used black
- use self.rootfinder, work only with DAE
- extracted _get_varlist_to_plot() method
- Fixed mypy and flake8 error. Used black
- Added mypy, bumped lock

### Fix

- fix scaling of rootfinder
- fix problems if only one var to plot
- fixed test after change in reinitialize_algebraic
- fix test after generate_exp_data() rework
- fix test after varlist.plot() rework
- rework variable plot using dataframe
- issue with recalculate_algebraic
- alow pip install -e, fix dataframe
- fixed initializtion of Constants in Model
- rounded time_grid of simulation to seconds
- all tests passed
- finished PE, working on tests
- Simulation is reformated, tests working
- fixed pytest usage of coverage
- updated poetry lock
- wrong .variance order
- fixed PE/MPC objective with algebraic var

## 0.3.0 (2021-05-16)

## moved_git (2021-01-03)

## 0.2.1 (2020-11-29)

## 0.2.0 (2020-11-29)

## 0.1.0 (2020-04-02)

## 0.0.5 (2019-12-13)

## 0.0.4 (2019-12-12)

## 0.0.3 (2019-11-27)

## 0.0.2 (2019-11-26)
