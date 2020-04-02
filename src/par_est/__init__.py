# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = 'unknown'
finally:
    del get_distribution, DistributionNotFound

from .variables import Variable, Parameter_variable, Algebraic_variable, State_variable, Control_variable, VariableList, Experimental_Data
from .model import Model
from .simulation import Simulator
from .optimization import ParameterEstimation, OptimalExperimentalDesign
