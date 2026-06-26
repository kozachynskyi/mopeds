# mopeds

<img align="right" src="https://git.tu-berlin.de/dbta/optimization/mopeds/-/raw/main/docs/logo.png" width="300px">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20939503.svg)](https://doi.org/10.5281/zenodo.20939503)

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

## Citation

If you use MOPEDS in your research, please cite the Zenodo archive or one of the accompanying publications.

> Kozachynskyi, V. (2026). *mopeds – Model based Parameter Estimation and Design of Experiments* (Version 0.11.0). Zenodo. DOI: [10.5281/zenodo.20939504](https://doi.org/10.5281/zenodo.20939504)

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
