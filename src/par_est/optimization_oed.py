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
        jac_scaled = args[0]
        obj = np.trace(np.linalg.inv(jac_scaled.T @ jac_scaled))

        return obj, jac_scaled

class CriteriaD(OED_objective):
    def eval(self, args):
        jac_scaled = args[0]
        obj = np.linalg.det(np.linalg.inv(jac_scaled.T @ jac_scaled))

        return obj, jac_scaled


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
        jac_scaled = self.jacobian_scaled_mx
        obj = ca.trace(ca.inv(jac_scaled.T @ jac_scaled))

        return obj, jac_scaled

    def _objective_A_fd(self):
        self._objective_func = CriteriaA("A", self.jacobian_scaled_mx)
        func_eval = self._objective_func(self.jacobian_scaled_mx)

        return func_eval[0], func_eval[1]

    def _objective_D_fd(self):
        self._objective_func = CriteriaD("D", self.jacobian_scaled_mx)
        func_eval = self._objective_func(self.jacobian_scaled_mx)

        return func_eval[0], func_eval[1]

    def optimize(self, scale=False, objective_function="A"):
        """Function to select optimization function"""
        # Scaling works unpredictably. It was shown during creation of VariableControlPiecewiseConstant
        if scale:
            raise NotImplementedError

        self.select_objective_function(objective_function)

        return self._optimize(scale)

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
                if not var.fixed:
                    if isinstance(var, VariableControlPiecewiseConstant):
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

        if len(self.varlist_parameter) == 0:
            raise ValueError("All parameters are fixed, OED is not possible")
        self.array_inverted_variances: np.ndarray = np.array(inverted_variances)
        self.array_inverted_std = np.sqrt(inverted_variances)

        self.parameter_values: list[float] = parameter_values

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
        time_grid_measurements: np.ndarray,
        time_grid_control_switch: np.ndarray | None = None,
        simulator_name: str = "idas",
        simulator_settings: dict = None,
        *,
        reinitialize_algebraic: bool = False,
        measurable_variables: list[str] | None = None,
    ) -> None:
        super().__init__(model, variable_list, simulator_name, simulator_settings)

        # User specified time_grid is used for initilizaiton of Simulators
        self._setup_timegrid(time_grid_measurements, time_grid_control_switch)

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

        # User specified time grid might not include every time_stamp of VariableControlPiecewiseConstant
        # So the corrected time_grid of Simulator is used in optimization
        self.time_grid_modified = self.list_simulators[0].time_grid_relative
        if not np.array_equal(self.time_grid_original, self.time_grid_modified):
            raise NotImplementedError

        self.solver_name: str = "ipopt"
        self.solver_settings: dict = {
            "verbose": False,
            # "enable_fd": True,
            # "enable_jacobian": False,
            # "enable_forward": False,
            # "enable_reverse": False,
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

    def _setup_timegrid(self, time_measurements, time_control_switch):
        if time_control_switch is None:
            time_control_switch = []
        time_grid = np.unique(list(time_measurements) + list(time_control_switch))
        self.time_grid_original: np.ndarray = time_grid


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
        self.generate_jacobian_function()

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
