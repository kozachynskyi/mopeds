import importlib.metadata

__version__ = importlib.metadata.version("par_est")

try:
    import acados_template
    import par_est.casados_integrator

    _ACADOS_SUPPORT = True
except ImportError:
    _ACADOS_SUPPORT = False

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
from .simulation_dynamic import Simulator
from .simulation_nle import SimulatorNLE
from .optimization import (
    Optimizer,
    ParameterEstimation,
    ParameterEstimationNLE,
    ParameterEstimationNLE_control,
)
from .optimization_oed import OptimalExperimentalDesign, OEDsettings
from .mpc import ModelPredictiveControl

import par_est.examples
import par_est.tools
