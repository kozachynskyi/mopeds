# Release Notes

## 0.11.1 (2026-06-27)

This release focuses on improving the quality and maintainability of `mopeds`.

## Highlights

- Added support for NumPy 2.5.
- Removed unused and broken legacy code.
- Added software citation support, including a Zenodo DOI and publication list.

## 0.11.0 (2026-06-21)

This release modernizes MOPEDS by updating its core dependencies and development tooling.

### Breaking Changes

- MOPEDS now requires CasADi 3.7.
- Support for CasADi 3.6 has been removed.

### Development Infrastructure

- Migrated from Poetry to uv for dependency management, development environments, and package builds.

## 0.10.3 (2026-06-20)

This release focuses on compatibility updates and development infrastructure improvements.

### Highlights

- Added support for Python 3.13.
- Added support for NumPy 2.x.
- Added automated CI testing for supported Python versions.
- Migrated project tooling from Poetry to uv.
- Added `newton` and `fast_newton` solvers for nonlinear equation systems (NLE).

### Other Improvements

- Improved plotting support in `ErrorAnalyzer`.
- Improved package compatibility checks and dependency constraints.

## 0.10.2 (2025-04-14)

This release introduces major improvements to parameter estimation, optimal experimental design, scaling, and model diagnostics.

### Highlights

#### Direct Optimization for NLE Models

Parameter Estimation (PE) and Optimal Experimental Design (OED) for nonlinear equation systems now support direct optimization in addition to the traditional rootfinder-based formulation.

#### Automatic Scaling Framework

A new scaling framework simplifies model setup and improves numerical robustness. Variables can now be scaled automatically based on their bounds.

#### Simplified Model Architecture

Models no longer store variable lists internally. Instead, they only maintain the variable ordering and associated CasADi symbols, resulting in a cleaner and more maintainable architecture.

#### Faster Repeated Parameter Estimation

Experimental data is now supplied as an NLP parameter, allowing repeated parameter estimation runs without rebuilding the optimization problem. Optimizers can also reuse solver instances, significantly reducing overhead for workflows such as Monte Carlo studies.

#### Error Analysis and Diagnostics

Added the new `ErrorAnalyzer` framework with support for:

- Prediction analysis
- Outlier detection
- Covariance analysis
- Dynamic models
- Covariance ellipse visualization

Additional helper functions were introduced:

- `tools.analyze_scaling_nle()`
- `pe.check_decision_bounds()`

#### API Improvements

- Unified simulator API for dynamic and steady-state models.
- `generate_varlist...()` for NLE workflows was renamed to `generate_artificial_data()`.
- `VariableConstant` now accepts multiple ignored inputs to simplify switching between constants and independent variables.
- `PE.calculate_objective_and_residuals()` now provides `df_all` with all simulated variables.

## 0.10.1 (2023-12-13)

### Highlights

- Added GitLab CI and automated testing infrastructure.
- Improved support for Python 3.12.
- Made constant variables immutable after creation.

## 0.10.0 (2023-12-09)

### Highlights

- Renamed `par_est` to `mopeds`.
- Open-sourced the project.

## 0.9.3a1 (2023-10-11)

### Highlights

- Added regularization techniques for parameter estimation.
- Added parameter scaling for regularized optimization.
- Expanded OED functionality for nonlinear equation systems.

## 0.9.0 (2023-04-23)

### Highlights

- Added ACADOS support for ODE and DAE simulation.
- Added identifiability analysis for dynamic models.
- Reworked parameter estimation internals for improved maintainability and consistency.

### Breaking Changes

- DAE simulator API changed from `simulate()` to `simulate_sym()`.

## 0.8.0 (2022-08-24)

### Highlights

- Major redesign of nonlinear equation system (NLE) simulation and parameter estimation.
- Added inference and uncertainty analysis tools.
- Improved parameter variance and covariance calculations.
- Added examples based on Bates & Watts, *Nonlinear Regression Analysis and Its Applications*.
