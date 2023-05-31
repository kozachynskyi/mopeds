from __future__ import annotations

import copy
import logging
from abc import abstractmethod
from collections.abc import Callable
from itertools import combinations
from typing import Sequence
from warnings import warn

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


class OED_base(Optimizer):
    def _objective_A(self):
        """A criteria"""
        raise NotImplementedError

    def optimize(self, scale=False, objective_function="A"):
        """Function to select optimization function"""
        # Scaling works unpredictably. It was shown during creation of VariableControlPiecewiseConstant
        if scale:
            raise NotImplementedError

        if objective_function == "A":
            self._objective = self._objective_A
            raise NotImplementedError(
                f"Objective function '{objective_function}' is not supported"
            )

        return self._optimize(scale)

    def change_parameter_values(self):
        """Change parameter values in simulator"""

    def calculate_objective_and_fim(
        self,
        controls: dict[str, float],
        objective_function: str = "A",
    ) -> dict[str, float | np.ndarray]:
        self._setup_scaling(False)
        if objective_function == "A":
            obj_f = self._objective_A()
        else:
            raise NotImplementedError

        decision_variables = self.varlist_decision.get_casadi_variables()
        casadi_function = ca.Function(
            "objective",
            [decision_variables],
            [obj_f[0], obj_f[1], self.simulate_all_mx],
            ["x"],
            ["f", "residuals", "y"],
        )

        selected_parameters = self.parameters_dict_to_list(parameters)
        res = casadi_function(x=selected_parameters)
        result_np = {
            "f": float(res["f"]),
            "residuals": res["residuals"].toarray(),
            "y": res["y"].toarray(),
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
                if not var.fixed:
                    if isinstance(var, VariableControlPiecewiseConstant):
                        for var_control in var.variable_list.values():
                            if not var_control.fixed:
                                self.varlist_decision.add_variable(var_control)
                    else:
                        self.varlist_decision.add_variable(var)
            elif isinstance(var, VariableParameter):
                if var.fixed is False:
                    self.varlist_parameters.add_variable(var)
                    parameter_values.append(var.value[0])

            elif isinstance(var, VariableState):
                if var.fixed is False:
                    inverted_variances.append(1 / var.variance)
                    self.names_of_measurements.append(var.name)

        self.array_inverted_variances: np.ndarray = np.array(inverted_variances)

        self.parameter_values: list[float] = parameter_values

    def generate_fim_function(self) -> None:
        """Combines simulate_sym() functions from simulator, and creates MX structure, that is used
        further in objective_function calculation"""
        parameter_variables = self.varlist_parameters.get_casadi_variables()

        if isinstance(self.list_simulators[0], Simulator):
            res_dict_name = "xf"
        else:
            raise NotImplementedError

        list_simulation_T = []

        for simulator in self.list_simulators:
            res_simulation = simulator.simulate_sym()

            list_simulation_T.append(res_simulation[res_dict_name].T)

        free_variables = ca.vcat([
                self.varlist_parameters.get_casadi_variables(),
                self.varlist_decision.get_casadi_variables(),
                ]
                )

        all_selected_measurements = ca.vcat(list_simulation_T).get(
            False, ca.Slice(), self.index_measurements_in_sim
        )
        self.simulate_all_function = ca.Function(
            "sim_all", [free_variables], [all_selected_measurements]
        )
        self.simulate_all_mx = self.simulate_all_function(free_variables)

        jacobian = {}
        jacobian_scaled = {}

        for index_measurement, meas_name in enumerate(self.names_of_measurements):
            jac_meas_mx = ca.jacobian(
                self.simulate_all_mx[:, index_measurement], parameter_variables
            )
            jac_meas_function = ca.Function(
                "jac_meas", [parameter_variables], [jac_meas_mx]
            )
            jac_meas_dm = jac_meas_function(self.parameter_values)
            jac_meas_scaled_dm = (
                jac_meas_dm * self.array_inverted_std[:, index_measurement]
            )
            jacobian[meas_name] = jac_meas_dm
            jacobian_scaled[meas_name] = jac_meas_scaled_dm

        jac_array = np.concatenate(list(jacobian.values()))
        jac_array_scaled = np.concatenate(list(jacobian_scaled.values()))
        raise NotImplementedError


class OptimalExperimentalDesign(OED_base):
    def __init__(
        self,
        model: Model,
        variable_list: list[VariableList],
        time_grid_relative: np.ndarray,
        simulator_name: str = "idas",
        simulator_settings: dict = None,
        *,
        reinitialize_algebraic: bool = False,
    ) -> None:
        super().__init__(model, variable_list, simulator_name, simulator_settings)

        # User specified time_grid is used for initilizaiton of Simulators
        self.time_grid_original: np.ndarray = time_grid_relative

        self._setup_simulator()
        self._setup_initialization()

        # User specified time grid might not include every time_stamp of VariableControlPiecewiseConstant
        # So the corrected time_grid of Simulator is used in optimization
        self.time_grid_modified = self.list_simulators[0].time_grid_relative
        if not np.array_equal(self.time_grid_original, self.time_grid_modified):
            raise NotImplementedError

        self.solver_name: str = "ipopt"
        self.solver_settings: dict = {
            "verbose": False,
            # "monitor": ["nlp_grad_f", "nlp_f"],
            "ipopt": {
                "max_iter": 100,
                # "print_level": 6,
            },
        }
        if reinitialize_algebraic:
            for sim in self.list_simulators:
                sim.calculate_algebraic_initials(apply_intials=True)

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
                self.time_grid_original,
                self.list_input_varlist[0],
                self.simulator_name,
                self.simulator_settings,
                simulate_jac=False,
            )
        ]

        for name in self.names_of_measurements:
            index = self.list_simulators[0].mapping_state_variables[name]
            self.index_measurements_in_sim.append(index)

        self.mapping_simulator_decisions: list[dict[int, int]] = [self.list_simulators[0].mapping_independent_variables]
        self.generate_fim_function()

    def _objective(
        self, analyze: bool = False, values: list[float] | None = None
    ) -> tuple[ca.MX | ca.DM, ca.MX | ca.DM]:
        """Calculates an A OED criteria, beacuse casadi cannot do other.
        "analyze" Flag and values are used for debugging.

        Args:
            analyze: used for debugging, calculates objective and covariance.
            values: this values are used as desicionb variables for calculating of the objective.
        """
        covariance_full = None
        # -1 ignores time point zero in self.time_grid
        num_time = len(self.time_grid_modified) - 1
        # +1 account for tau variable in Simulator class
        num_param = len(self.model.varlist_independent) + 1
        num_state = len(self.model.varlist_state)

        result_simulation = self.list_simulators[0].simulate_jac()
        result_jacobian = result_simulation["jac_xf_p"]

        # Used only for debugging
        if analyze is True:
            jacobian_full = None
            covariance_all = []
            objective_all = []

            self._setup_scaling(False)
            evaluate = ca.Function(
                "eval_fim",
                [self.varlist_decision.get_casadi_variables()],
                [result_jacobian],
            )
            evaluate_sim = ca.Function(
                "eval_fim",
                [self.varlist_decision.get_casadi_variables()],
                [result_simulation["xf"]],
            )
            if values is None:
                result_jacobian = evaluate(self.guess)
                result_sim = evaluate_sim(self.guess)
            else:
                result_jacobian = evaluate(values)
                result_sim = evaluate_sim(values)  # noqa: F841

        # Simulation returns jacobian that has to be split, to get jac at each time point.
        # list_jacobian_at_timepoint contains a list of that jacobians.
        split_vector = np.linspace(0, num_time, num_time + 1, dtype=int) * num_param
        list_jacobian_at_timepoint = ca.horzsplit(result_jacobian, split_vector)

        # Jacobian is also scaled based on parameter values.
        parameter_scaling = ca.repmat(ca.DM(self.parameter_values).T, num_state, 1)

        # For jacobian at each timepoint, select only unfixed parameters.
        # Afterwards scale the jacobiand and calculate a parameter covariabce matrix.
        # Finally, sum covariances at every time point.
        for jacobian in list_jacobian_at_timepoint:
            jacobian_selected = jacobian.get(
                False, self.index_all_states, self.select_independent
            )
            jacobian_selected = jacobian_selected * parameter_scaling

            cov_at_timepoint = (
                jacobian_selected.T
                @ np.diag(self.inverted_variances)
                @ jacobian_selected
            )

            if analyze:
                covariance_all.append(cov_at_timepoint)
                objective_all.append(ca.trace(ca.inv(cov_at_timepoint)))

            if covariance_full is None:
                covariance_full = cov_at_timepoint
            else:
                covariance_full = covariance_full + cov_at_timepoint

        if analyze:
            for jacobian in list_jacobian_at_timepoint:
                jacobian_selected = jacobian.get(
                    False, self.index_all_states, self.select_independent
                )
                jacobian_selected = jacobian_selected * parameter_scaling

                if jacobian_full is None:
                    jacobian_full = jacobian_selected
                else:
                    jacobian_full = ca.vertcat(jacobian_full, jacobian_selected)

        error = ca.trace(ca.inv(covariance_full))

        if analyze:
            return error, covariance_full, jacobian_full, covariance_all, objective_all  # type: ignore

        return error, covariance_full
