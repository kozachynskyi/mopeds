## 0.10.0 (2023-12-07)

### Feat

- rework oed jacobian to work with new casadi syntax
- update dependencies, allow python 3.12
- move GN working PE to examples
- fix regularization tests, and setup ridge regularization
- move WIP tihkhnov tuning example in separate file
- prepare for merge, fix tests
- WIP impleementaton of pe with algebraic variables
- add warning message while using ipopt as nle solver during PE
- allow ipopt as solver of NLE equations system

### Fix

- support python >3.9
- update factory() syntax of casadi
- tools bug inttoduced with PiecewiseConstantCOntrols
- undo changes, prepare for merge
- fixed typo

### Refactor

- rename src folder
- rework poetry to support 3.12
- remove pyDOE artifacts
- turn off casados_integrator
- remove Jinja2 and pyDOE as dependencies
- rework to avoid deprecation warnings of casadi and numpy
- move example file

## 0.9.3a1 (2023-10-11)

### Feat

- add test for tikhonov
- add tikhonov WIP
- add parameter scaling of jacobian for regularization
- tools generate_varlist supports perturbation
- reworked regularization techniques

### Fix

- fix lopez regularization

## 0.9.3a0 (2023-09-24)

### Feat

- gitlab ci turn off, add poetry config
- OEDNLE supports simulator class as input
- add cstr NLE example in par_est
- OED Nle takes into acount previous measurements
- add OED for NLE WIP
- oed simulator supports ca.MX time_grid

### Fix

- add FixedSampling setting for OED to __init__
- oed.generate_data fixed Piecewise Controls
- pe.parameter_analysis catches error if hessian not calculated

### Refactor

- OED move _optimize to BaseOED

## 0.9.2.alpha1 (2023-08-31)

### Fix

- add FixedGridSamping as default mode for OED

## 0.9.2.alpha (2023-08-31)

### Feat

- VariablePiecewiseConstant reworked lb, ub and guess
- VariablePiecewiseConstant supports ca.MX as time_grid
- oed supports Adaptive based sampling with piecewise controls
- oed setups equality constraint after initialization
- rework oed.generate_exp_data
- add oed.separate_and_check_controls()
- raise error of supplied time_sp or weight_ control is wrong
- VariablePiecewiseConstant testing time_grid from ca.MX
- categorize OED solution strategies Optimal, Adaptive etc
- add deluca2016
- varlist.plot() fills NA before plotting, to plot Controls
- sim.generate_exp_data() returns PiecewiseConstantControls
- **OED**: add option to add custom objective function
- support for measruement weights
- add initial version of hoang2014 from quaglio
- **OED**: genereate_exp_data works with flexible sampling time
- wip add reference to "parent" variable in VaribleControl
- PieceswiseConstantControls are plotted in varlist
- **OED**: add methods to create varlist for PE, and simulate Sim
- **Sim**: make raised error more informative
- **OED**: add simulate method to OED
- **OED**: support for OED when t_sp is given, and t_sw is selected
- **OED**: support for initial values of state variables
- **OED**: jac returned is not scaled by parameters
- **OED**: add parameter scaling to objective function
- add OEDsettings, refactor OED to be more readable
- first working version of OED with flexible time_grid
- **OED**: don't allow unfixed PiecewiseControl vars with horizon
- remove redundant code
- **OED**: make time_grid_contorl an optional argument
- add D criteria
- utilities create_grid returns meschgrid for plotting
- and yeast model, working on quaglio
- add A criteria with FD

### Fix

- OED without provided settings works as expected
- oed.generate_data ignores redundant time controls
- OED generate_exp_data typo
- **OED**: generate_exp_data if time0 of state variable is decision var
- **OED**: typo in Objective_D
- typo in OED
- typos in OED
- add _parameter_scaling to Callback
- API change create_grid
- typo OED
- fix yeast model
- typo

### Refactor

- remove unused import
- modify error message
- temp_backup
- move tests around
- add correct typing
- **OED**: remove old self.optimize()
- **Sim**: change iteration over self.time_relative via range
- **Sim**: remove unused self.integrator
- **PE**: jacobian scaling is clearly formulated
- separate oed in opitmization_oed
- separate Simulator in NLE and Dynamic

## 0.9.2 (2023-05-02)

### Fix

- revert false caclulations of wls

## 0.9.1 (2023-04-23)

### Fix

- **PE**: identifiability use correct FIM, order of parameters irrelevant

## 0.9.0 (2023-04-23)

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

## 0.8.3 (2022-12-09)

### Feat

- **SimDAE**: add self.mapping_*_variables dictionary
- **Varlist**: add show argument to plot()
- **PENLE**: calclulate_objective_residuals returns sim values

### Fix

- **PENLE**: but at recursion of indetifiability analysis
- **PENLE**: bug when parameters are not fixed as in varlist_decision

### Refactor

- **Variable**: allow lists to specify time_grid

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

## 0.7.4 (2022-07-23)

## 0.7.3 (2022-07-23)

### Refactor

- remove python 3.7 support, let isort skip __init__

## 0.7.2 (2022-07-19)

### Feat

- example 5-min-tut tested and implemented
- add parameter_identifiability*yao and *eigenvalue
- add peNLE.calculate_sensitivity_and_fim() method
- add pe.calculate_ols_value() to analyze objective
- add dict with variables and values as output of optimize()

### Fix

- use FIM scaled with meas var for eingevalue ranking
- fixed naming issues f <-> y and m <-> model in 5min_tutorial.rst
- resolved an error in examples/nle/cstr_collocation_sim.py
- wrong variables names while reporting optimize()

### Refactor

- add .env folder to .gitignore
- fix version to 0.7.1, add ITWM deploy Readme
- deleted unnecessary file
- remove redundant print statements
- create sim.simulate_sym_unfixed() method

## 0.7.1 (2022-07-04)

### Feat

- add SimulaotClass to PENLE
- WIP - fixed covaraince calculation for PENLE
- multistart sorts result by objective function
- optimizer prints result as dict after finishing
- add tqdm to signalize multistart progress
- add do_once() method for monkey-patching
- add tqdm as dependency

### Fix

- removed decorator @notypecheck

### Refactor

- add warning using "ipopt" with NLE
- do not delete tmp_file after close
- move import matplotlib inside methods
- add new line for better output
- remove breakpoins from examples

## 0.7.0 (2022-06-19)

### Feat

- value.setter alternative to dataframe.iloc[0]
- add .show() method to show dataframe to HTML
- remove pandasgui dependency, add Jinja2

### Fix

- fix DeprecationWarning pandas.Index.get_loc

### Refactor

- add type annotations
- remove Variable.get_data_opcua() method
- import from par_est directly
- remove unused code

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

### Fix

- fix example change API
- set_bounds() works for Control Variables

### Refactor

- fix type of list_simulators
- use __future__ for type anotations
- _objective() returns list
- rework how Simulatore-PE mapping works
- scaling SimulatorNLE from PENLE
- black8

## 0.6.0 (2022-03-21)

### Feat

- update dependency versions, drop Python3.7

## 0.5.2 (2022-03-15)

### Feat

- remove setup.py
- add __version__
- added berty NLE example

### Fix

- add importlib_metadata for python<3.8

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
- Remove pip install -e support, update versions

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

### Refactor

- added plot argument to plot_simulation
- used black
- used black
- use self.rootfinder, work only with DAE
- extracted _get_varlist_to_plot() method
- Fixed mypy and flake8 error. Used black
- Added mypy, bumped lock

## 0.3.0 (2021-05-16)

## 0.2.1 (2020-11-29)

## 0.2.0 (2020-11-29)

## 0.1.0 (2020-04-02)

## 0.0.5 (2019-12-13)

## 0.0.4 (2019-12-12)

## 0.0.3 (2019-11-27)

## 0.0.2 (2019-11-26)
