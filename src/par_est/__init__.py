try:
    import importlib.metadata as importlib_metadata
except ModuleNotFoundError:
    import importlib_metadata

__version__ = importlib_metadata.version(__name__)

from .utilities import MXPickler, show_html_from_dataframe
from .variables import (
    Variable,
    VariableParameter,
    VariableAlgebraic,
    VariableState,
    VariableConstant,
    VariableControlPiecewiseConstant,
    VariableControl,
    VariableList,
)
from .variables import BadVariableError
from .variables import ORIGIN_TS
from .model import Model
from .simulation import Simulator, SimulatorNLE
from .optimization import (
    Optimizer,
    ParameterEstimation,
    OptimalExperimentalDesign,
    ParameterEstimationNLE,
    ParameterEstimationNLE_control,
)
from .mpc import ModelPredictiveControl

import par_est.examples
