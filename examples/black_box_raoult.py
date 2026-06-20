import mopeds
import casadi as ca

try:
    import sympy
except ImportError as exc:
    raise ImportError(
        "SymPy is required to run examples/black_box_dynamic.py. "
        "Install it with `pip install sympy`."
    ) from exc

import numpy as np

def get_jacobian_from_sympy():
    symbols = [ "p", "T", "a", "b", "c"]
    p, T, a, b, c = sympy.symbols(symbols)
    vars = [p, T, a, b, c]
    eq = (p-10**(a - b / (T + c)))

    for v in vars:
        print(sympy.simplify(sympy.diff(eq, v)))

# print(get_jacobian_from_sympy())
# breakpoint()

class ArhenuisJac(ca.Callback):
    def __init__(self, opts={}):
        ca.Callback.__init__(self)
        self.construct("jac_f", opts)

    def get_n_in(self): return 2
    def get_n_out(self): return 1

    def get_sparsity_in(self,i):
        if i==0: # nominal input
          return ca.Sparsity.dense(5,1)
        elif i==1: # nominal output
          return ca.Sparsity(1,1)

    def get_sparsity_out(self,i):
        return ca.Sparsity.dense(1,5)

    def eval(self, arg):
        p, T, a, b, c = ca.vertsplit(arg[0])
        res = ca.DM(1,5)
        res[0,0] = 1
        res[0,1] = -10**(a - b/(T + c))*b*np.log(10)/(T + c)**2
        res[0,2] = -10**(a - b/(T + c))*np.log(10)
        res[0,3] = 10**(a - b/(T + c))*np.log(10)/(T + c)
        res[0,4] = -10**(a - b/(T + c))*b*np.log(10)/(T + c)**2

        return [res]

# print(ca.log(10))
fun_jac = ArhenuisJac()
class Arhenius(ca.Callback):
    def __init__(self, name, opts={}):
        # opts = {"enable_fd":True}
        ca.Callback.__init__(self)
        self.construct(name, opts)

    def get_n_in(self):
        return 1

    def get_n_out(self):
        return 1

    def get_sparsity_in(self,i):
        return ca.Sparsity.dense(5,1)

    def get_sparsity_out(self,i):
        return ca.Sparsity.dense(1,1)

    def eval(self, arg):
        p, T, a, b, c = ca.vertsplit(arg[0])
        res = (p-10**(a - b / (T + c)))
        return [res]

    def has_jacobian(self): return True
    def get_jacobian(self,name,inames,onames,opts):

        # You are required to keep a reference alive to the returned Callback object
        self.jac_callback = fun_jac
        return self.jac_callback


fun_arh = Arhenius("c1")

def raults():
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4&Type=ANTOINE&Plot=on
    variable_list = mopeds.VariableList()
    # fmt:off
    variable_list.add_variable(mopeds.VariableParameter("e0_A_c1", 3.5595, 3, 4))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_B_c1", 643.748))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_C_c1", -198.043))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_A_c2", 4.92531, 4, 6))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_B_c2", 1432.526))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_C_c2", -61.819))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P_LV_o_c1", 1.4144774629101626, 0.1, 10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_P_LV_o_c2", 0.3584844184269499, 0.1, 10))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_c2", 0.5, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c1", 0.78888, 0, 1))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c2", 0.21112, 0, 1))  # noqa: E501

    variable_list.add_variable(mopeds.VariableControl("e0_x_c1", 0.5, 0, 1))  # noqa: E501

    if True:
        variable_list.add_variable(mopeds.VariableControl("e0_T", 346.4149, 100, 500))  # noqa: E501
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_P", 1.013, 0.1, 6))  # noqa: E501
    else:
        variable_list.add_variable(mopeds.VariableAlgebraic("e0_T", 346.4149, 100, 500))  # noqa: E501
        variable_list.add_variable(mopeds.VariableControl("e0_P", 1.013, 0.1, 6))  # noqa: E501

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

    EQ_alg1 = ((e0_y_c1*e0_P)-((e0_x_c1*e0_P_LV_o_c1)))  # noqa: E501,E226
    EQ_alg2 = ((e0_y_c2*e0_P)-((e0_x_c2*e0_P_LV_o_c2)))  # noqa: E501,E226
    EQ_alg3 = (1.0-(((e0_x_c1+e0_x_c2))))  # noqa: E501,E226
    EQ_alg4 = (1.0-(((e0_y_c1+e0_y_c2))))  # noqa: E501,E226

    if True:
        EQ_alg5 = fun_arh(ca.vcat([e0_P_LV_o_c1, e0_T, e0_A_c1, e0_B_c1, e0_C_c1]))
        EQ_alg6 = fun_arh(ca.vcat([e0_P_LV_o_c2, e0_T, e0_A_c2, e0_B_c2, e0_C_c2]))
    else:
        EQ_alg5 = (e0_P_LV_o_c1-10**(e0_A_c1 - e0_B_c1 / (e0_T + e0_C_c1)))  # noqa: E501,E226
        EQ_alg6 = (e0_P_LV_o_c2-10**(e0_A_c2 - e0_B_c2 / (e0_T + e0_C_c2)))  # noqa: E501,E226

    list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6]  # noqa: E501
    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m

vl, m = raults()
# sim = mopeds.SimulatorNLE(m, vl)

vl_list, true_params, _ = mopeds.tools.generate_varlist_with_data_NLE(m, vl, {"e0_x_c1": [0.1, 0.9, 10]}, perturbate=False, measurement_names=["e0_y_c1"])

for vl_i in vl_list:
    vl_i["e0_A_c1"].fixed = False
    vl_i["e0_A_c2"].fixed = True
    vl_i["e0_A_c1"].guess = 3.5
    vl_i["e0_A_c2"].guess = 4

sim = mopeds.SimulatorNLE(m, vl_i)
print(sim.simulate(unfixed_variables=true_params)[2].dataframe)

pe = mopeds.ParameterEstimationNLE(m, vl_list)
pe.solver_settings["ipopt"]["hessian_approximation"] = "limited-memory"
res = pe.optimize(direct_optimization=True)
print(pe.check_result_bounds(res))
print(res["x_dict"])
