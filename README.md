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

See:

- [Release Notes](RELEASE_NOTES.md) for user-facing release summaries and migration notes.
- [Changelog](CHANGELOG.md) for the complete automatically generated change history.


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
