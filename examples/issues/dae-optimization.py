try:
    import mopeds
    import numpy as np

    def raults():
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on
        # https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175&Mask=4&Type=ANTOINE&Plot=on
        variable_list = mopeds.VariableList()
        # fmt:off
        # e0_x_c1 = ca.MX.sym("e0_x_c1")
        variable_list.add_variable(mopeds.VariableState("e0_x_c1", 0.1, 0, 1))  # noqa: E501
        variable_list.add_variable(mopeds.VariableState("var_state2", 0, 0, 1))  # noqa: E501

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
        eq_s2 = e0_y_c1
        EQ_alg1 = ((e0_y_c1*e0_P)-((e0_x_c1*e0_P_LV_o_c1)))  # noqa: E501,E226
        EQ_alg2 = ((e0_y_c2*e0_P)-((e0_x_c2*e0_P_LV_o_c2)))  # noqa: E501,E226
        EQ_alg3 = (1.0-(((e0_x_c1+e0_x_c2))))  # noqa: E501,E226
        EQ_alg4 = (1.0-(((e0_y_c1+e0_y_c2))))  # noqa: E501,E226

        EQ_alg5 = (e0_P_LV_o_c1-10**(e0_A_c1 - e0_B_c1 / (e0_T + e0_C_c1)))  # noqa: E501,E226
        EQ_alg6 = (e0_P_LV_o_c2-10**(e0_A_c2 - e0_B_c2 / (e0_T + e0_C_c2)))  # noqa: E501,E226

        list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6]  # noqa: E501
        # fmt:on

        m.add_equations_differential([eq_s, eq_s2])
        m.add_equations_algebraic(list_algebraic_equations)

        return variable_list, m
    vl, m = raults()
    grid = np.linspace(0, 8, 2)
    sim = mopeds.Simulator(m, grid, vl)
    res_mopeds = sim.simulate(algebraic=True)[2]
except Exception:
    pass


import casadi as ca
import numpy as np
import copy


ISSUE_1_trigger = False
ISSUE_2_trigger = False

def integrate_over_time_grid(time_grid, integrator, y_0, z_0, p_0):
    x0 = y_0
    z0 = z_0
    res_states = []
    res_algebraic = []
    previos_t = 0
    num_steps = time_grid[1:].shape[0] - 1

    for time_index in range(num_steps + 1):
        time_step = time_grid[time_index + 1]
        res_integration = integrator(
            x0=x0,
            z0=z0,
            p=ca.vertcat(time_step - previos_t, *p_0),
        )

        previos_t = time_step
        x0 = res_integration["xf"][:, 1]
        z0 = res_integration["zf"][:, 1]

        if time_index == 0:
            res_states.append(res_integration["xf"][:, 0])
            res_algebraic.append(res_integration["zf"][:, 0])

        res_states.append(x0)
        res_algebraic.append(z0)

    res_states = ca.hcat(res_states)
    res_algebraic = ca.hcat(res_algebraic)

    res = ca.vcat([res_states, res_algebraic])
    return res

# VariableState
e0_x_c1 = ca.MX.sym("e0_x_c1")
var_state2 = ca.MX.sym("var_state2")
VariableState_0 = [0.1, 0.0]
VariableState_mx = [e0_x_c1, var_state2]
# VariableAlgebraic
e0_P_LV_o_c1 = ca.MX.sym("e0_P_LV_o_c1")
e0_P_LV_o_c2 = ca.MX.sym("e0_P_LV_o_c2")
e0_x_c2 = ca.MX.sym("e0_x_c2")
e0_y_c1 = ca.MX.sym("e0_y_c1")
e0_y_c2 = ca.MX.sym("e0_y_c2")
e0_P = ca.MX.sym("e0_P")
VariableAlgebraic_0 = [1.4144774629101626, 0.3584844184269499, 0.5, 0.78888, 0.21112, 1.013]
VariableAlgebraic_mx = [e0_P_LV_o_c1, e0_P_LV_o_c2, e0_x_c2, e0_y_c1, e0_y_c2, e0_P]
# VariableParameter
e0_A_c1 = ca.MX.sym("e0_A_c1")
e0_B_c1 = ca.MX.sym("e0_B_c1")
e0_C_c1 = ca.MX.sym("e0_C_c1")
e0_A_c2 = ca.MX.sym("e0_A_c2")
e0_B_c2 = ca.MX.sym("e0_B_c2")
e0_C_c2 = ca.MX.sym("e0_C_c2")
e0_T = ca.MX.sym("e0_T")
VariableParameter_0 = [3.5595, 643.748, -198.043, 4.92531, 1432.526, -61.819, 346.4149]
VariableParameter_mx = [e0_A_c1, e0_B_c1, e0_C_c1, e0_A_c2, e0_B_c2, e0_C_c2, e0_T]

EQ_diff1 = 0.1
if ISSUE_1_trigger:
    EQ_diff2 = 1
else:
    EQ_diff2 = e0_y_c1
EQ_alg1 = ((e0_y_c1*e0_P)-((e0_x_c1*e0_P_LV_o_c1)))  # noqa: E501,E226
EQ_alg2 = ((e0_y_c2*e0_P)-((e0_x_c2*e0_P_LV_o_c2)))  # noqa: E501,E226
EQ_alg3 = (1.0-(((e0_x_c1+e0_x_c2))))  # noqa: E501,E226
EQ_alg4 = (1.0-(((e0_y_c1+e0_y_c2))))  # noqa: E501,E226

EQ_alg5 = (e0_P_LV_o_c1-10**(e0_A_c1 - e0_B_c1 / (e0_T + e0_C_c1)))  # noqa: E501,E226
EQ_alg6 = (e0_P_LV_o_c2-10**(e0_A_c2 - e0_B_c2 / (e0_T + e0_C_c2)))  # noqa: E501,E226

list_diff_equations = [EQ_diff1, EQ_diff2]  # noqa: E501
list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6]  # noqa: E501

tau = ca.MX.sym("tau")
dae_system = {
    "x": ca.vcat(VariableState_mx),
    "p": ca.vertcat(tau, *VariableParameter_mx),
    "ode": ca.vcat(list_diff_equations) * tau,
    "alg": ca.vcat(list_algebraic_equations),
    "z": ca.vcat(VariableAlgebraic_mx),
}

integrator = ca.integrator( "integrator", "idas", dae_system, 0, [0, 1], {})

res = integrate_over_time_grid(grid, integrator, VariableState_0, VariableAlgebraic_0, VariableParameter_0)

unfixed_symbols = copy.deepcopy(VariableParameter_0)
unfixed_symbols[0] = e0_A_c1
res_sym = integrate_over_time_grid(grid, integrator, VariableState_0, VariableAlgebraic_0, unfixed_symbols)

if ISSUE_2_trigger:  # estimated based on "e0_y_c1" algebraic variable
    index = 5
else:  # estimate based on "var_state2" variable
    index = 1

exp_data = res[index, :]
function = ca.Function("f", [e0_A_c1], [res_sym])

obj_MX = ca.sumsqr(exp_data - function(e0_A_c1)[index, :])
obj_function = ca.Function("obj", [e0_A_c1], [obj_MX])

nlpsol_dict = { "x": e0_A_c1, "f": obj_MX }

solver = ca.nlpsol("par_est", "ipopt", nlpsol_dict, {"monitor": "nlp_grad_f", "ipopt": {"max_iter": 1, "print_level": 0}})

nlpsol_args = { "x0": 3.4, "lbx": 3, "ubx": 4}

res_pe = solver.call(nlpsol_args)
print(f"Expected 3.5595, estimated: {res_pe['x']}")

final_res = function(3.5595)
res_obj = obj_function(3.5595)
obj_jac = obj_function.jacobian()
jac = obj_jac(3.4, 0)
print(jac)

print(res_obj)

# Test if model correctly impelemnted, do not submit to issue
assert np.isclose(res_mopeds.dataframe.to_numpy(), res.T).all()
