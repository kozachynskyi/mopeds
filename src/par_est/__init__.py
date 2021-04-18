# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from .variables import (
    Variable,
    VariableParameter,
    VariableAlgebraic,
    VariableState,
    VariableConstant,
    VariableControlPiecewiseConstant,
    VariableControl,
    VariableList,
    ExperimentData,
)
from .variables import BadVariableError
from .tools import MXPickler
from .model import Model
from .simulation import Simulator
from .optimization import (
    Optimizer,
    ParameterEstimation,
    OptimalExperimentalDesign,
)
