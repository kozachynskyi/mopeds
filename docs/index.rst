############################################
par_est: Parameter Estimation of DAE and NLE
############################################

Introduction
============

|par_est| is a tool developed with a following goal: create a top-level wrapper for |casadi|, to make modelling, simulation and optimization more simple. At current state of the package, package has following features:

- Abstraction levels: Variable, VariableList, Model, Simulation and Optimization
- Supported models: |ode|, |dae| and |nle|
- Supported optimization routines: parameter estimation, optimal experimental design

Requirements
------------

Python 3.8 or 3.9


Installation
------------

* Ask for acess to `git_par_est`_
* Simple installation::

    $ pip install git+https://git.tu-berlin.de/vovakozach/pe_oed_casadi

* Install as developer::

    $ git clone https://git.tu-berlin.de/vovakozach/pe_oed_casadi
    $ cd pe_oed_casadi
    $ portry install OR pip install .



TOC
===

.. toctree::
   :maxdepth: 2

   5min_tutorial.rst
   available_settings.rst
   changelog.rst

Glossary
========

.. glossary::
    NLE
        Systems of nonlinear algebraic equations
    DAE
        Differential-algebraic system of equations
    ODE
        Ordinary differential equations
