import importlib.metadata

__version__ = importlib.metadata.version("mopeds")

from .utilities import MXPickler, show_html_from_dataframe
from .variables import (
    get_options,
    set_options,
    options,
    Variable,
    VariableParameter,
    VariableAlgebraic,
    VariableState,
    VariableConstant,
    VariableControlPiecewiseConstant,
    VariableControl,
    VariableList,
    _consistent_scaling_decorator,
)
from .variables import BadVariableError
from .variables import ORIGIN_TS
from .model import Model
from .simulation_dynamic import Simulator
from .simulation_nle import SimulatorNLE
from .optimization import (
    Optimizer,
    ParameterEstimation,
    ParameterEstimationNLE,
    ParameterEstimationNLE_control,
)
from .optimization_oed import (
    OptimalExperimentalDesign,
    OED_objective,
    OptimalSampling,
    AdaptiveOptimalSampling,
    AdaptiveSampling,
    FixedGridSampling,
    OptimalExperimentalDesign_NLE,
)

import mopeds.examples
import mopeds.tools
