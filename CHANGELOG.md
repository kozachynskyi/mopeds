## 0.9.0 (2023-04-23)

### Refactor

- reroll back example file changes
- refactor the parameter analysis, return marginal CI
- **PE**: move paramter_analysis and yao to PE_base
- **SimDAE**: simplify simulate_t0
- **SimDAE**: simplify simulate_dae_recalcluate_algebraic
- **Sim**: reposition methods to set default similator settings
- modify .gitignore file
- **PE**: rearange some methods
- add example for NLE cstr
- move some methods around
- **PENLE**: remove unused atribute
- **PE**: move some methods
- **PE**: remove redundancy
- **PE**: use setup_simulator_mapping
- **PE**: create PE_Base class and move methods around

### Fix

- **PE**: typo. how parameter covariance matrix is calculated
- **PE**: fix WLS objective function , added 1/2
- rework paramter variance for multiresponce data
- acados import typo2
- typo acados import
- add fix to bug that happens on OptigodMarkI
- remove typo
- remove hessian-approximation option from ipopt
- generalize experiments scaling for PE
- **PENLE**: fix potential bug with sorting of variances

### Feat

- **PE**: add option to not show plot while doing parameter_analysis
- add _debug functions
- **PE**: add property DOF to calculate NumMeas-NumPar
- add identifiability analysis to ODE/DAE
- WIP fair function for PE
- **PE**: acados simulator C code generated only for first simulator
- add acados dae solver support
- **Model**: add .name attribute
- include casados_integrator with modifications
- **Sim**: add simulate_unfixed method for ODE/DAE
- add log scaling to sampling, andd tests
- **Variables**: PiecewiceConstant.value setting is supported
- **PE**: deprecate "use_idas_constraints"
- add example time derivative
- speed up VariableList.dataframe
- **Sim**: add change_independent_variables()
- **PE**: move to modern api of PENLE
- **PEDAE**: remove objective_alg()
- **SIM_ODE**: rename self.simulate() to self.simulate_sym()

## 0.8.3 (2022-12-09)

### Feat

- **SimDAE**: add self.mapping_*_variables dictionary
- **Varlist**: add show argument to plot()
- **PENLE**: calclulate_objective_residuals returns sim values

### Refactor

- **Variable**: allow lists to specify time_grid

### Fix

- **PENLE**: but at recursion of indetifiability analysis
- **PENLE**: bug when parameters are not fixed as in varlist_decision

## 0.8.2 (2022-11-02)

### Fix

- make python3.8 work, fix __future__ annotations
- remove np.int depreciation warning

## 0.8.1 (2022-10-27)

### Fix

- **PENLE**: calculate_objective_and_residual fixed scaling bug

## 0.8.0 (2022-08-24)

### Feat

- **NLE**: add option to GenerateDataNle to supply responce names
- rounds print output of bates examples
- adds a documentation for calculate_inference_bounds in optimization.py
- inference band for NLE PE using exp OR artificial data
- added examples for calculate_inference_bound based on Bates et al.
- added inference examples
- added puromycin examples from Nonlinear regression analysis (M. Bates)
- **PENLE**: rework how calculate_sensitivity() works
- **PEnle**: rename calculate_ols_value to calculate_objective_and_residual
- **PeNLE**: new ols and wls objective calculations (moredimensional)
- **PENLE**: add index_measurements_in_sim list
- **NLE**: data generator now returns true_parameters
- **PENLE**: add new way to data_array and data_mask calculation
- **SimNLE**: add simulate() method with option to select variables
- **SimNLE**: add change_independent_variables method
- **SimNLE**: add mapping of algebraic and independent variables
- **SimNLE**: add contains_unfixed argument to Simulator
- **utilities**: rework how generate_varlist_for_optimizer NLE works
- pe.calculate_sensitivity -> return meas_covaraince
- return all residuals, not only nonzero ones

### Refactor

- **examples**: add dostring
- fix typo preturbate to perturbate
- fix not working examples
- combine examples, remove redundancy
- remove redundant type conversion to array
- hide import scipy.stats in method
- remove redundancy and assert
- rename perturbation_mode
- avoid redundancy
- more readable code
- remove unused imports
- code formatting
- **NLE**: move all examples from bates to separate file
- **ExampleNLE**: add example file for puromycin
- pre-commit
- **PENLE**: rework calculate_ols_value
- **examples**: add simple_mixer NLE example
- add example from seminar 2022-08-02
- clean up pe.parameter_analysis

### Fix

- **MPC**: raise NotImplementedError (optimizer is not supported)
- fix OLS dictionary
- **NLEinference**: remove redundant seed attribute
- **Vars**: VariableConstant.casadi_var returns self.value
- **PENLE**: parameter_analysis plotting and calculations new API
- **PENLE**: parameter identifiability yao new API
- **PENLE**: fix how parameter_dict_to_list works, make it robust
- **PENLE**: raise error, when optimize(obj_func=) str is not supported
- **PEnle**: fix bug with _reset_scaling
- **SimNLE**: generate_exp_data -> var.variance is propagated further
- **PENLE**: array_data_new nan replaced with 0
- **PENLE**: add raise error
- (parameter_analysis) -> correctly scaled parameter covariance

## 0.7.4 (2022-07-23)

## 0.7.3 (2022-07-23)

### Refactor

- remove python 3.7 support, let isort skip __init__

## 0.7.2 (2022-07-19)

### Refactor

- add .env folder to .gitignore
- fix version to 0.7.1, add ITWM deploy Readme
- deleted unnecessary file
- remove redundant print statements
- create sim.simulate_sym_unfixed() method

### Fix

- use FIM scaled with meas var for eingevalue ranking
- fixed naming issues f <-> y and m <-> model in 5min_tutorial.rst
- resolved an error in examples/nle/cstr_collocation_sim.py
- wrong variables names while reporting optimize()

### Feat

- example 5-min-tut tested and implemented
- add parameter_identifiability*yao and *eigenvalue
- add peNLE.calculate_sensitivity_and_fim() method
- add pe.calculate_ols_value() to analyze objective
- add dict with variables and values as output of optimize()

## 0.7.1 (2022-07-04)

### Feat

- add SimulaotClass to PENLE
- WIP - fixed covaraince calculation for PENLE
- multistart sorts result by objective function
- optimizer prints result as dict after finishing
- add tqdm to signalize multistart progress
- add do_once() method for monkey-patching
- add tqdm as dependency

### Refactor

- add warning using "ipopt" with NLE
- do not delete tmp_file after close
- move import matplotlib inside methods
- add new line for better output
- remove breakpoins from examples

### Fix

- removed decorator @notypecheck

## 0.7.0 (2022-06-19)

### Feat

- value.setter alternative to dataframe.iloc[0]
- add .show() method to show dataframe to HTML
- remove pandasgui dependency, add Jinja2

### Refactor

- add type annotations
- remove Variable.get_data_opcua() method
- import from par_est directly
- remove unused code

### Fix

- fix DeprecationWarning pandas.Index.get_loc

## 0.6.2 (2022-06-14)

### Feat

- WIP - calculation of parameter variance

## 0.6.1 (2022-06-14)

### Feat

- for PE_NLE add OLS and WLS objective function
- add jac_calculation and variance extract NLE
- PE for NLE with variable controls
- add examples from Literatures
- Latin Hypercube sampling is default, pyDOE
- add SimulatorNLE.jacobian
- add sphinx requirement
- add scaling to NLEOptimization
- add bounds for rootfinder, make bounds optional

### Refactor

- fix type of list_simulators
- use __future__ for type anotations
- _objective() returns list
- rework how Simulatore-PE mapping works
- scaling SimulatorNLE from PENLE
- black8

### Fix

- fix example change API
- set_bounds() works for Control Variables

## 0.6.0 (2022-03-21)

### Feat

- update dependency versions, drop Python3.7

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
