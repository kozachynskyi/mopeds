import mopeds
import casadi as ca
import numpy as np
import copy

# a = mopeds.VariableAlgebraic("a", -11, 3e1, 5e1)
# print(a.scale_from_original(a.lower_bound))
# print(a.scale_from_original(a.upper_bound))
# print(var := a.scale_from_original(35))
# print(a.scale_to_original(var))
# breakpoint()

# def scale(lb, value, ub):
#     v = 1  / (ub - lb)
#     r = 2 - (ub / (ub - lb))
#     return (((np.array([lb, value, ub]) - r) / v) , ((np.array([lb, value, ub]) * v) + r))

# print(scale(-1,10,21))
# print(scale(-1,1.5, 21))
# # print(scale(1e4,1e5,2e9))
# breakpoint()

class PE(mopeds.ParameterEstimationNLE):
    def _objective_ols(self):
        """Objective function is a trace(Z.T * Z), where Z is a residual matrix with shape:
        numRows -> amount of supplied experiments, numCol -> amount of variables that have measurements
        If experiments do not supply a measurement for one of the measurements, self.array_data_mask will
        have 0 as the respective element of the martix, otherwise 1"""
        residuals = (
            (self.simulate_all_mx - self.array_data)
            * self.array_data_mask
            * np.sqrt(self.experiments_scale)
        )
        objective = ca.sumsqr(residuals)

        return objective, residuals

# Gill2008 Practical opitmization isbn: 978-0-12-283952-8
def example8dot7(dae=False) -> tuple[
    mopeds.VariableList, mopeds.Model, list[mopeds.VariableList]
]:
    variable_list = mopeds.variables.VariableList()  # Preallocate variable_list

    # variable_list.add_variable(mopeds.VariableAlgebraic("F", 0.008, 0.008, 0.027))
    variable_list.add_variable(mopeds.VariableAlgebraic("F", 0.000008, 0.000008, 0.000027))
    # variable_list.add_variable(mopeds.VariableAlgebraic("F", 1, 0.9, 1.1))

    # variable_list.add_variable(mopeds.VariableAlgebraic("F", 0.008, -ca.inf, ca.inf))

    variable_list.add_variable(mopeds.VariableControl("x1", 2, 0.1, 10))
    variable_list.add_variable(mopeds.VariableControl("x2", 1e-2, 1e-2, 1e-1))

    variable_list.add_variable(mopeds.VariableParameter("b", 3, 1, 7))
    variable_list.add_variable(mopeds.VariableParameter("c", 3, 1, 7))

    if dae:
        variable_list.add_variable(mopeds.VariableState("tt", 0, 0, 0.0004))

    m = mopeds.Model(variable_list)  # adding all variables to the model

    F = m.varlist_all["F"].casadi_var  # noqa: E501
    x1 = m.varlist_all["x1"].casadi_var  # noqa: E501
    x2 = m.varlist_all["x2"].casadi_var  # noqa: E501
    b = m.varlist_all["b"].casadi_var  # noqa: E501
    c = m.varlist_all["c"].casadi_var  # noqa: E501

    f = variable_list["F"]

    v = 1
    r = 0

    eq1 = (F / v) - r - ((x1**b) * (x2**c))


    m.add_equations_algebraic([eq1])  # adding the equations to model

    if dae:
        m.add_equations_differential([F])

    return variable_list, m

def nle_scaling():
    vl, m = example8dot7()
    sim = mopeds.SimulatorNLE(m, vl)

    exp, truepar = mopeds.tools.generate_varlist_with_data_NLE(m, vl, {"x1":[2,3,3]}, perturbate=False)
    for i in exp:
        print(i.dataframe)

    exp[0]["b"].fixed = False
    exp[0]["c"].fixed = False
    exp[0]["b"].guess = 3.1
    exp[0]["c"].guess = 3.1

    pe = PE(m, exp)
    v = pe.list_simulators[0]
    par = {"b":3.1, "c":3.1} 
    vv = pe.calculate_objective_and_residual({"b":3.1, "c":3.1})
    breakpoint()


    # pe.solver_settings["ipopt"]["max_iter"] = 4
    # pe.prepare_nle()
    res = pe.optimize(scale=False, objective_function="ols")
    print(res)
    breakpoint()

def nle_jac():
    vl, m = example8dot7()
    sim = mopeds.SimulatorNLE(m, vl)

    exp, truepar = mopeds.tools.generate_varlist_with_data_NLE(m, vl, {"x1":[2,3,3]}, perturbate=False)
    for i in exp:
        print(i.dataframe)

    exp[0]["b"].fixed = False
    exp[0]["c"].fixed = False
    exp[0]["b"].guess = 3.1
    exp[0]["c"].guess = 3.1
    par = {"b":3.0, "c":3.0} 

    mopeds.variables.VARIABLE_SCALING = False
    pe = mopeds.ParameterEstimationNLE(m, exp)
    v = pe.calculate_sensitivity_and_fim(par)
    mopeds.variables.VARIABLE_SCALING = True
    pe = mopeds.ParameterEstimationNLE(m, exp)
    vv = pe.calculate_sensitivity_and_fim(par)
    print(vv["jac_full"] / v["jac_full"])
    breakpoint()


    # pe.solver_settings["ipopt"]["max_iter"] = 4
    # pe.prepare_nle()
    res = pe.optimize(scale=False, objective_function="ols")
    print(res)
    breakpoint()


if __name__ == "__main__":
    nle_jac()

    vl, m = example8dot7(True)
    grid = np.linspace(0,100,5)
    sim = mopeds.Simulator(m, grid, vl)
    res1 = sim.simulate(algebraic=True)[2]
    sim.change_independent_variables({"x1": 3})
    res2 = sim.simulate(algebraic=True)[2]

    exp = []
    for i, res in enumerate([res1, res2]):
        exp_i = copy.deepcopy(vl)
        if i == 1:
            exp_i["x1"].value = 3
        exp_i["tt"].dataframe = res["tt"].dataframe
        print(res.dataframe)

        exp_i["b"].fixed = False
        exp_i["c"].fixed = False
        exp_i["b"].guess = 3.1
        exp_i["c"].guess = 3.1
        exp.append(exp_i)

    # pe = PE(m, exp)
    pe = mopeds.ParameterEstimation(m, exp)
    v = pe.list_simulators[0]
    vv = pe.calculate_objective_and_residual({"b":3.0, "c":3.0})
    breakpoint()


    # pe.solver_settings["ipopt"]["max_iter"] = 4
    # pe.prepare_nle()
    res = pe.optimize(scale=False, objective_function="ols")
    print(res)
