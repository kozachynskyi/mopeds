import mopeds
import casadi as ca

try:
    import sympy
except ImportError as exc:
    pass

import numpy as np
import copy

mopeds.set_options(variable_scaling=False)


# def get_jacobian_from_sympy():
#     symbols = ["p", "T", "a", "b", "c"]
#     p, T, a, b, c = sympy.symbols(symbols)
#     vars = [p, T, a, b, c]
#     eq = p - 10 ** (a - b / (T + c))

#     for v in vars:
#         print(sympy.simplify(sympy.diff(eq, v)))
# print(get_jacobian_from_sympy())
# breakpoint()


class ArhenuisJac(ca.Callback):
    def __init__(self, opts={}):
        ca.Callback.__init__(self)
        self.construct("jac_f", opts)

    def get_n_in(self):
        return 2

    def get_n_out(self):
        return 1

    def get_sparsity_in(self, i):
        if i == 0:  # nominal input
            return ca.Sparsity.dense(5, 1)
        elif i == 1:  # nominal output
            return ca.Sparsity(1, 1)

    def get_sparsity_out(self, i):
        return ca.Sparsity.dense(1, 5)

    def eval(self, arg):
        p, T, a, b, c = ca.vertsplit(arg[0])
        res = ca.DM(1, 5)
        res[0, 0] = 1
        res[0, 1] = -(10 ** (a - b / (T + c))) * b * np.log(10) / (T + c) ** 2
        res[0, 2] = -(10 ** (a - b / (T + c))) * np.log(10)
        res[0, 3] = 10 ** (a - b / (T + c)) * np.log(10) / (T + c)
        res[0, 4] = -(10 ** (a - b / (T + c))) * b * np.log(10) / (T + c) ** 2

        return [res]


# print(ca.log(10))
fun_jac = ArhenuisJac()


class Arhenius(ca.Callback):
    def __init__(self, name, opts={"enable_fd": False}):
        # opts = {"enable_fd":True}
        ca.Callback.__init__(self)
        self.construct(name, opts)

    def get_n_in(self):
        return 1

    def get_n_out(self):
        return 1

    def get_sparsity_in(self, i):
        return ca.Sparsity.dense(5, 1)

    def get_sparsity_out(self, i):
        return ca.Sparsity.dense(1, 1)

    def eval(self, arg):
        p, T, a, b, c = ca.vertsplit(arg[0])
        res = p - 10 ** (a - b / (T + c))
        return [res]

    def has_jacobian(self):
        return False

    def get_jacobian(self, name, inames, onames, opts):

        # You are required to keep a reference alive to the returned Callback object
        self.jac_callback = fun_jac
        return self.jac_callback


fun_arh = Arhenius("c1")


def raults():
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4&Type=ANTOINE&Plot=on
    variable_list = mopeds.VariableList()
    # fmt:off
    variable_list.add_variable(mopeds.VariableState("e0_x_c1", 0.1, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableState("test", 0, 0, 1))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P_LV_o_c1", 1.4144774629101626, 0.1, 10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P_LV_o_c2", 0.3584844184269499, 0.1, 10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_c1", 3.5595, 3, 4))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_B_c1", 643.748))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_C_c1", -198.043))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_c2", 4.92531, 4, 6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_B_c2", 1432.526))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_C_c2", -61.819))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_c2", 0.5, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c1", 0.78888, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c2", 0.21112, 0, 1))  # noqa: E501

    # variable_list.add_variable(mopeds.VariableControlPiecewiseConstant("e0_x_c1", 0.1, 0, 1))  # noqa: E501

    variable_list.add_variable(mopeds.VariableControl("e0_T", 346.4149, 100, 500))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P", 1.013, 0.1, 6))  # noqa: E501

    m = mopeds.Model(variable_list)

    e0_P = m.varlist_all["e0_P"].casadi_var  # noqa: E501
    e0_x_c1 = m.varlist_all["e0_x_c1"].casadi_var  # noqa: E501
    e0_P_LV_o_c1 = m.varlist_all["e0_P_LV_o_c1"].casadi_var  # noqa: E501
    e0_P_LV_o_c2 = m.varlist_all["e0_P_LV_o_c2"].casadi_var  # noqa: E501
    e0_A_c1 = m.varlist_all["e0_A_c1"].casadi_var  # noqa: E501
    e0_A_c2 = m.varlist_all["e0_A_c2"].casadi_var  # noqa: E501
    e0_B_c1 = m.varlist_all["e0_B_c1"].casadi_var  # noqa: E501
    e0_B_c2 = m.varlist_all["e0_B_c2"].casadi_var  # noqa: E501
    e0_C_c1 = m.varlist_all["e0_C_c1"].casadi_var  # noqa: E501
    e0_C_c2 = m.varlist_all["e0_C_c2"].casadi_var  # noqa: E501
    e0_T = m.varlist_all["e0_T"].casadi_var  # noqa: E501
    e0_x_c2 = m.varlist_all["e0_x_c2"].casadi_var  # noqa: E501
    e0_y_c1 = m.varlist_all["e0_y_c1"].casadi_var  # noqa: E501
    e0_y_c2 = m.varlist_all["e0_y_c2"].casadi_var  # noqa: E501

    eq_s = 0.1
    eq_s2 = 1
    EQ_alg1 = ((e0_y_c1*e0_P)-((e0_x_c1*e0_P_LV_o_c1)))  # noqa: E501,E226
    EQ_alg2 = ((e0_y_c2*e0_P)-((e0_x_c2*e0_P_LV_o_c2)))  # noqa: E501,E226
    EQ_alg3 = (1.0-(((e0_x_c1+e0_x_c2))))  # noqa: E501,E226
    EQ_alg4 = (1.0-(((e0_y_c1+e0_y_c2))))  # noqa: E501,E226

    if False:
        EQ_alg5 = fun_arh(ca.vcat([e0_P_LV_o_c1, e0_T, e0_A_c1, e0_B_c1, e0_C_c1]))
        EQ_alg6 = fun_arh(ca.vcat([e0_P_LV_o_c2, e0_T, e0_A_c2, e0_B_c2, e0_C_c2]))
    else:
        EQ_alg5 = (e0_P_LV_o_c1-10**(e0_A_c1 - e0_B_c1 / (e0_T + e0_C_c1)))  # noqa: E501,E226
        EQ_alg6 = (e0_P_LV_o_c2-10**(e0_A_c2 - e0_B_c2 / (e0_T + e0_C_c2)))  # noqa: E501,E226

    list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6]  # noqa: E501
    # fmt:on

    m.add_equations_differential([eq_s, eq_s2])
    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


vl, m = raults()
grid = np.linspace(0, 8, 20)

# times = np.linspace(1, 30, 12)
# data = np.linspace(0.1, 0.9, 12)
# vl["e0_x_c1"].expand_horizon(times, data)
sim = mopeds.Simulator(m, grid, vl, integrator_settings={"expand": True})
res = sim.simulate(algebraic=True)[2]
print(res.dataframe)


exp_data = copy.deepcopy(vl)
rng = np.random.default_rng(0)

names = ["test"]  # , "e0_x_c2", "e0_P_LV_o_c1"]
# names = ["test", "e0_y_c1"]
# for var_name in names:
#     exp_data[var_name].dataframe = res[var_name].dataframe
#     # exp_data["e0_T"].dataframe = exp_data["e0_T"].dataframe.to_numpy() + rng.normal([0]*exp_data.dataframe.shape[0], 1)
#     v = exp_data[var_name].dataframe
#     # exp_data["e0_T"].dataframe = (v.T  + rng.normal([0]*exp_data.dataframe.shape[0], 1)).T
for var in res.values():
    exp_data[var.name].dataframe = res[var.name].dataframe


exp_data["e0_A_c1"].fixed = False
# exp_data["e0_A_c2"].fixed = False
exp_data["e0_A_c1"].guess = 3.5
# exp_data["e0_A_c2"].guess = 4

pe = mopeds.ParameterEstimation(
    m, [exp_data], simulator_settings={"expand": True, "enable_fd": False}
)
# pe.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"
# v = pe.calculate_sensitivity_and_fim({"e0_A_c1": 3.5})
# vv = pe.calculate_objective_and_residual({"e0_A_c1": 3.1})
# vvv = pe.calculate_objective_and_residual({"e0_A_c1": 3.5595})
# breakpoint()
res = pe.optimize()
print(res["x_dict"])
