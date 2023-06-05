from __future__ import annotations

import copy
import logging
from abc import abstractmethod
from collections.abc import Callable
from itertools import combinations
from typing import Sequence
from warnings import warn
from dataclasses import dataclass

import casadi as ca
import numpy as np
import pandas as pd
from scipy import linalg
from tqdm import tqdm

from par_est import (
    Model,
    Simulator,
    VariableControl,
    VariableControlPiecewiseConstant,
    VariableList,
    VariableParameter,
    VariableState,
    _ACADOS_SUPPORT,
    Optimizer
)

if _ACADOS_SUPPORT:
    from par_est import casados_integrator

@dataclass
class OEDsettings:
    """Class to setup complex OED optimization problems"""
    max_time_experiment: float = 1
    min_sampling_delay: float = 0.1
    num_sampling_times: int = 3

class OED_objective(ca.Callback):
    def __init__(self, name, jac, opts={}):
        opts["enable_jacobian"] = False
        opts["enable_forward"] = False
        opts["enable_reverse"] = False
        opts["enable_fd"] = True

        ca.Callback.__init__(self)
        self.nin = jac.shape
        self.construct(name, opts)

    def get_n_in(self): return 1
    def get_n_out(self): return 2

    def eval(self, args):
        raise NotImplementedError

    def get_sparsity_in(self,i):
        return ca.Sparsity.dense(*self.nin)

    def get_sparsity_out(self,i):
        if i == 0:
            return ca.Sparsity.dense(1)
        elif i == 1:
            return ca.Sparsity.dense(*self.nin)


class CriteriaA(OED_objective):
    def eval(self, args):
        jac = args[0]
        jac_scaled = args[0] * self._parameter_scaling
        obj = np.trace(np.linalg.inv(jac_scaled.T @ jac_scaled))

        return obj, jac

class CriteriaD(OED_objective):
    def eval(self, args):
        jac = args[0]
        jac_scaled = args[0] * self._parameter_scaling
        obj = np.linalg.det(np.linalg.inv(jac_scaled.T @ jac_scaled))

        return obj, jac

class OED_base(Optimizer):
    def select_objective_function(self, objective_function_name: str):
        if objective_function_name == "A":
            self._objective = self._objective_A
        elif objective_function_name == "A_fd":
            self._objective = self._objective_A_fd
        elif objective_function_name == "D":
            self._objective = self._objective_D_fd
        else:
            raise NotImplementedError(
                f"Objective function '{objective_function_name}' is not supported"
            )
        return self._objective

    def _objective_A(self):
        """A criteria"""
        jac = self.jacobian_scaled_mx
        jac_scaled = jac * self._parameter_scaling
        obj = ca.trace(ca.inv(jac_scaled.T @ jac_scaled))

        return obj, jac

    def _objective_A_fd(self):
        self._objective_func = CriteriaA("A", self.jacobian_scaled_mx)
        func_eval = self._objective_func(self.jacobian_scaled_mx)

        return func_eval[0], func_eval[1]

    def _objective_D_fd(self):
        self._objective_func = CriteriaD("D", self.jacobian_scaled_mx)
        func_eval = self._objective_func(self.jacobian_scaled_mx)

        return func_eval[0], func_eval[1]

    def optimize(self, scale:float = 1, objective_function="A"):
        """Function to select optimization function"""
        self.select_objective_function(objective_function)

        return self._optimize(scale)

    @property
    def _parameter_scaling(self):
        parameter_scaling = ca.repmat(self.parameter_values, 1, self.jacobian_scaled_mx.shape[0]).T
        return parameter_scaling

    def change_parameter_values(self):
        """Change parameter values in simulator"""

    def calculate_objective_and_jacobian(
        self,
        controls: dict[str, float],
        objective_function: str = "A",
    ) -> dict[str, float | np.ndarray]:
        self._setup_scaling(False)

        obj_f = self.select_objective_function(objective_function)()

        decision_variables = self.varlist_decision.get_casadi_variables()
        casadi_function = ca.Function(
            "objective",
            [decision_variables],
            [obj_f[0], obj_f[1]],
            ["x"],
            ["f", "jac"],
        )

        selected_parameters = self.variables_dict_to_list(controls)
        res = casadi_function(x=selected_parameters)
        result_np = {
            "f": float(res["f"]),
            "jac": res["jac"].toarray(),
        }

        return result_np

    def _setup_varlist_decision(self):
        parameter_values = []
        inverted_variances = []
        self.names_of_measurements = []

        for variable_name in self.model.varlist_all.keys():
            try:
                var = self.list_input_varlist[0][variable_name]
            except KeyError:
                continue

            if isinstance(var, VariableControl):
                if var.fixed is False:
                    if isinstance(var, VariableControlPiecewiseConstant):
                        if not len(var.variable_list) == 1:
                            raise NotImplementedError("Piecewise constant controls with time grid are not supported")
                        for var_control in var.variable_list.values():
                            if not var_control.fixed:
                                self.varlist_decision.add_variable(var_control)
                    else:
                        self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableParameter):
                if var.fixed is False:
                    self.varlist_parameter.add_variable(var)
                    parameter_values.append(var.value[0])

            elif isinstance(var, VariableState):
                if var.name in self.list_measureable_variables:
                    inverted_variances.append(1 / var.variance)
                    self.names_of_measurements.append(var.name)
                # if var.fixed is False:
                #     inverted_variances.append(1 / var.variance)
                #     self.names_of_measurements.append(var.name)

        for time_var in self.varlist_timegrid.values():
            self.varlist_decision.add_variable(time_var)

        if len(self.varlist_parameter) == 0:
            raise ValueError("All parameters are fixed, OED is not possible")
        self.array_inverted_variances: np.ndarray = np.array(inverted_variances)
        self.array_inverted_std = np.sqrt(inverted_variances)

        self.parameter_values = np.array(parameter_values)

    def generate_jacobian_function(self) -> None:
        """Combines simulate_sym() functions from simulator, and creates MX structure, that is used
        further in objective_function calculation"""
        parameter_variables = self.varlist_parameter.get_casadi_variables()

        if isinstance(self.list_simulators[0], Simulator):
            res_dict_name = "xf"
        else:
            raise NotImplementedError

        res_simulation = self.list_simulators[0].simulate_sym()[res_dict_name].T

        parameter_variables = self.varlist_parameter.get_casadi_variables()
        decision_variables = self.varlist_decision.get_casadi_variables()

        all_selected_measurements = res_simulation.get(
            False, ca.Slice(), self.index_measurements_in_sim
        )

        all_selected_measurements = all_selected_measurements.get(
            False, self.index_time_grid, ca.Slice()
        )

        self.simulate_all_function = ca.Function(
            "sim_all", [self.varlist_parameter.get_casadi_variables(), self.varlist_decision.get_casadi_variables()], [all_selected_measurements]
        )
        self.simulate_all_mx = self.simulate_all_function(parameter_variables, decision_variables)

        jacobian = {}
        jacobian_scaled = {}

        for index_measurement, meas_name in enumerate(self.names_of_measurements):
            jac_meas_mx = ca.jacobian(
                self.simulate_all_mx[:, index_measurement], parameter_variables
            )
            jac_meas_function = ca.Function(
                "jac_meas", [parameter_variables, decision_variables], [jac_meas_mx]
            )
            jac_meas_mx = jac_meas_function(self.parameter_values, decision_variables)

            jac_meas_scaled_mx = (
                jac_meas_mx * self.array_inverted_std[index_measurement]
            )
            jacobian[meas_name] = jac_meas_mx
            jacobian_scaled[meas_name] = jac_meas_scaled_mx

        jac_array = ca.vcat(list(jacobian.values()))
        jac_array_scaled = ca.vcat(list(jacobian_scaled.values()))

        self.jacobian_mx = jac_array
        self.jacobian_scaled_mx = jac_array_scaled


class OptimalExperimentalDesign(OED_base):
    def __init__(
        self,
        model: Model,
        variable_list: list[VariableList],
        time_grid_measurements: np.ndarray | OEDsettings,
        simulator_name: str = "idas",
        simulator_settings: dict = None,
        *,
        reinitialize_algebraic: bool = False,
        measurable_variables: list[str] | None = None,
    ) -> None:
        super().__init__(model, variable_list, simulator_name, simulator_settings)
        self.varlist_timegrid = VariableList()

        # User specified time_grid is used for initilizaiton of Simulators
        if isinstance(time_grid_measurements, OEDsettings):
            self._initialize_from_settings(time_grid_measurements)
        else:
            if not time_grid_measurements[0] == 0:
                raise ValueError("Time grid should start with 0")
            self.time_grid_measurements = np.sort(time_grid_measurements)

        if measurable_variables is None:
            self.list_measureable_variables = list(self.model.varlist_state.keys())
        else:
            self.list_measureable_variables = []
            # Do this so variable names are sorted as expected
            for var_name in self.model.varlist_state.keys():
                if var_name in measurable_variables:
                    self.list_measureable_variables.append(var_name)

        self._setup_simulator()
        self._setup_initialization()

        self.solver_name: str = "ipopt"
        self.solver_settings: dict = {
            "verbose": False,
            # "monitor": ["nlp_grad_f", "nlp_f"],
            "ipopt": {
                "max_iter": 300,
                "hessian_approximation": "limited-memory"
                # "print_level": 6,
            },
        }
        if reinitialize_algebraic:
            for sim in self.list_simulators:
                sim.calculate_algebraic_initials(apply_intials=True)

    def _initialize_from_settings(self, settings:OEDsettings):
        self.max_time_experiment = settings.max_time_experiment
        self.min_sampling_delay = settings.min_sampling_delay

        initial_guess = np.linspace(settings.min_sampling_delay, settings.max_time_experiment, settings.num_sampling_times)

        for i, guess in enumerate(initial_guess):
            new_var = VariableControl("time_sp" + str(i), guess, self.min_sampling_delay, self.max_time_experiment)
            self.varlist_timegrid.add_variable(new_var)
        time_grid_measurements = self.varlist_timegrid.get_casadi_variables()
        self.time_grid_measurements = ca.vcat([0, time_grid_measurements])

    def _setup_timegrid(self):
        # Simulator time_grid might have time_steps, at which "measueremnt" is not done,
        # and jacobian calculation for FIM is not needed. This selection is done via self.index_time_grid
        time_grid_simulator = self.list_simulators[0].time_grid_relative

        self.index_time_grid = []

        if isinstance(time_grid_simulator, ca.MX):
            if time_grid_simulator.shape == self.time_grid_measurements.shape:
                self.index_time_grid = [time for time in range(time_grid_simulator.shape[0] - 1)]
            else:
                raise NotImplementedError
        else:
            # Starting from [1:] because 0 is always first element in time_grid
            for time_index, time in enumerate(time_grid_simulator[1:]):
                if time in self.time_grid_measurements:
                    self.index_time_grid.append(time_index)

    def _setup_simulator(self, *, use_idas_constraints: bool = False) -> None:
        """Initializes simulator class. Parameter variables are fixed, and an index of an unfixed
        parameter is saved in self.select_independent list.
        This list is used during the calculation of the objective, to ignore jacobian of fixed parameters.
        self.index_all_states is used additionaly to self.select_independent list to get required jacobian.
        """
        self._setup_varlist_decision()

        self.index_measurements_in_sim = []

        self.list_simulators: list[Simulator] = [
            Simulator(
                self.model,
                self.time_grid_measurements,
                self.list_input_varlist[0],
                self.simulator_name,
                self.simulator_settings,
                simulate_jac=False,
            )
        ]

        for name in self.names_of_measurements:
            index = self.list_simulators[0].mapping_state_variables[name]
            self.index_measurements_in_sim.append(index)

        self._setup_timegrid()

        self.mapping_simulator_decisions: list[dict[int, int]] = [self.list_simulators[0].mapping_independent_variables]
        self.generate_jacobian_function()

    def _setup_equality_constraints(self):
        casadi_vars = self.varlist_timegrid.get_casadi_variables()
        g = []
        self.lower_bound_g = []
        self.upper_bound_g = []

        for i in range(len(self.varlist_timegrid) - 1):
            g.append(casadi_vars[i+1] - casadi_vars[i])
            self.lower_bound_g.append(self.min_sampling_delay)
            self.upper_bound_g.append(self.max_time_experiment)

        self.equality_constraints = ca.vcat(g)

    def _optimize(self, scale: float) -> dict[str, ca.DM | ca.MX]:
        """Runs optimizer, uses scaling if needed. Returned values is scaled back.
        Scaling should be done before setting a solver and solver settings."""
        self._setup_equality_constraints()

        self.solver: ca.Function = ca.nlpsol(
            "solver",
            self.solver_name,
            {
                "x": self.varlist_decision.get_casadi_variables(),
                "f": self._objective()[0] * scale,
                "g": self.equality_constraints,
            },
            self.solver_settings,
        )

        res_solver = self.solver(
            x0=self.guess,
            lbx=self.lower_bound,
            ubx=self.upper_bound,
            lbg=self.lower_bound_g,
            ubg=self.upper_bound_g,
        )

        res_solver["x"] = res_solver["x"]

        res_dict = {}
        for solution, var_name in zip(
            res_solver["x"].toarray(), list(self.varlist_decision.keys())
        ):
            res_dict[var_name] = float(solution[0])

        res_solver["x_dict"] = res_dict
        self.reset_acados()

        return res_solver
