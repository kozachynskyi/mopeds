import copy
import casadi as ca
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 24

import mopeds
mopeds.set_options(variable_scaling=False)


def initialize_problem():  # noqa: C901

    variable_list = mopeds.VariableList()
    # fmt:off

    variable_list.add_variable(mopeds.VariableConstant("e0_A_c1", 4.42448))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_A_c2", 4.20772))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_c1", 1312.253))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_B_c2", 1233.129))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_c1", -32.445))  # noqa: E501
    variable_list.add_variable(mopeds.VariableConstant("e0_C_c2", -40.953))  # noqa: E501

    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_gamma_c1", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_greek_gamma_c2", 1.0, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_LV_o_c1", 1.0, 0.5, 10.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_p_LV_o_c2", 1.0, 0.5, 10.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_x_c2", 0.4, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c1", 0.7, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_y_c2", 0.3, 0.0, 1.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableAlgebraic("e0_T", 330.0, 273.0, 373.0))  # noqa: E501

    variable_list.add_variable(mopeds.VariableParameter("e0_greek_lambda_c1", -0.8404, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableParameter("e0_greek_lambda_c2", -0.561, -1.0E9, 1.0E9))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_p", 1.0, 0.5, 10.0))  # noqa: E501
    variable_list.add_variable(mopeds.VariableControl("e0_x_c1", 0.6, 0.0, 1.0))  # noqa: E501


    m = mopeds.Model(variable_list)

    e0_greek_lambda_c1 = m.varlist_all["e0_greek_lambda_c1"].casadi_var  # noqa: E501
    e0_greek_lambda_c2 = m.varlist_all["e0_greek_lambda_c2"].casadi_var  # noqa: E501
    e0_p = m.varlist_all["e0_p"].casadi_var  # noqa: E501
    e0_x_c1 = m.varlist_all["e0_x_c1"].casadi_var  # noqa: E501
    e0_A_c1 = m.varlist_all["e0_A_c1"].casadi_var  # noqa: E501
    e0_A_c2 = m.varlist_all["e0_A_c2"].casadi_var  # noqa: E501
    e0_B_c1 = m.varlist_all["e0_B_c1"].casadi_var  # noqa: E501
    e0_B_c2 = m.varlist_all["e0_B_c2"].casadi_var  # noqa: E501
    e0_C_c1 = m.varlist_all["e0_C_c1"].casadi_var  # noqa: E501
    e0_C_c2 = m.varlist_all["e0_C_c2"].casadi_var  # noqa: E501
    e0_greek_gamma_c1 = m.varlist_all["e0_greek_gamma_c1"].casadi_var  # noqa: E501
    e0_greek_gamma_c2 = m.varlist_all["e0_greek_gamma_c2"].casadi_var  # noqa: E501
    e0_p_LV_o_c1 = m.varlist_all["e0_p_LV_o_c1"].casadi_var  # noqa: E501
    e0_p_LV_o_c2 = m.varlist_all["e0_p_LV_o_c2"].casadi_var  # noqa: E501
    e0_x_c2 = m.varlist_all["e0_x_c2"].casadi_var  # noqa: E501
    e0_y_c1 = m.varlist_all["e0_y_c1"].casadi_var  # noqa: E501
    e0_y_c2 = m.varlist_all["e0_y_c2"].casadi_var  # noqa: E501
    e0_T = m.varlist_all["e0_T"].casadi_var  # noqa: E501

    EQ_alg1 = ((e0_y_c1*e0_p)-((e0_x_c1*(e0_greek_gamma_c1*e0_p_LV_o_c1))))  # noqa: E501,E226
    EQ_alg2 = ((e0_y_c2*e0_p)-((e0_x_c2*(e0_greek_gamma_c2*e0_p_LV_o_c2))))  # noqa: E501,E226
    EQ_alg3 = (((e0_x_c1+e0_x_c2))-(1.0))  # noqa: E501,E226
    EQ_alg4 = (((e0_y_c1+e0_y_c2))-(1.0))  # noqa: E501,E226
    EQ_alg5 = (e0_p_LV_o_c1-(((10.0))**(1.0*(e0_A_c1-(e0_B_c1/(e0_T+e0_C_c1))))))  # noqa: E501,E226
    EQ_alg6 = (e0_p_LV_o_c2-(((10.0))**(1.0*(e0_A_c2-(e0_B_c2/(e0_T+e0_C_c2))))))  # noqa: E501,E226
    EQ_alg7 = (e0_greek_gamma_c1-(ca.exp((((e0_greek_lambda_c1+(2.0*(((e0_greek_lambda_c2-e0_greek_lambda_c1))*e0_x_c1))))*((e0_x_c2))**(1.0*2.0)))))  # noqa: E501,E226
    EQ_alg8 = (e0_greek_gamma_c2-(ca.exp((((e0_greek_lambda_c2+(2.0*(((e0_greek_lambda_c1-e0_greek_lambda_c2))*e0_x_c2))))*((e0_x_c1))**(1.0*2.0)))))  # noqa: E501,E226

    list_algebraic_equations = [EQ_alg1, EQ_alg2, EQ_alg3, EQ_alg4, EQ_alg5, EQ_alg6, EQ_alg7, EQ_alg8, ]  # noqa: E501

    # fmt:on

    m.add_equations_algebraic(list_algebraic_equations)

    return variable_list, m


if __name__ == "__main__":

    variable_list, m = initialize_problem()
    true_params = {}
    true_params["e0_greek_lambda_c1"] = variable_list["e0_greek_lambda_c1"].value[0]
    true_params["e0_greek_lambda_c2"] = variable_list["e0_greek_lambda_c2"].value[0]

    # Create simulation Object
    sim = mopeds.SimulatorNLE(m, variable_list)
    res = sim.simulate()
    # print(res[2].dataframe[["e0_T", "e0_x_c2", "e0_y_c2"]])

    def plot_vle_diagram(sim):
        """Plots the VLE diagram for fixed pressure"""
        y = []
        T = []
        x =  np.linspace(0.0, 1, 100)
        for xx in x:
            sim.change_independent_variables({"e0_x_c1": xx})
            res = sim.simulate()[2]
            y.append(res["e0_y_c1"].value[0])
            T.append(res["e0_T"].value[0])
        fig, ax = plt.subplots(1,1)
        plt.plot(x, T)
        plt.plot(y, T)
        plt.ylabel("e0_T")
        plt.xlabel("e0_x_c1 / e0_y_c1")
        return fig, ax

    if False:
        plot_vle_diagram(sim)
        plt.show()

    MEAS_STD = 0.005
    def generate_artificial_data(sim):
        """Generate artificial data by simulating model and perturbating results"""
        list_exp_data = []
        y = []
        T = []
        x =  np.round(np.linspace(0.01, 0.99, 5), 4)
        rng = np.random.default_rng(0)
        y_exp = []
        for xx in x:
            sim.change_independent_variables({"e0_x_c1": xx})
            res = sim.simulate()[2]
            y = res["e0_y_c1"].value[0]
            y_random = round(rng.normal(y, MEAS_STD), 4)
            y_exp.append(y_random)

            vl_i = copy.deepcopy(variable_list)
            vl_i["e0_y_c1"].value = y_random
            vl_i["e0_x_c1"].value = xx
            list_exp_data.append(vl_i)

            T.append(round(res["e0_T"].value[0], 2))

        print(f"y_data = {repr(y_exp)}")
        print(f"x_data = {repr(list(x))}")
        print(f"T_data = {repr(list(T))}")

        return y_exp, x, T, list_exp_data


    # Generate artificial data
    generate_artificial_data(sim)

    y_data = [0.0059, 0.226, 0.5626, 0.833, 0.9924]
    x_data = [0.01, 0.255, 0.5, 0.745, 0.99]
    T_data = [334.16, 337.1, 336.57, 333.12, 329.19]

    if False:
        fig, ax = plot_vle_diagram(sim)
        ax.scatter(y_data, T_data, c="r", s=100)
        plt.show()

    # Add new experiments with OED
    if True:
        variable_list_oed = copy.deepcopy(variable_list)
        variable_list_oed["e0_x_c1"].fixed = False
        variable_list_oed["e0_greek_lambda_c1"].fixed = False
        variable_list_oed["e0_greek_lambda_c2"].fixed = False
        previous_measurements = [{"e0_x_c1": i} for i in x_data]
        new_x_data = copy.deepcopy(x_data)
        for i in range(2):
            rng = np.random.default_rng(0)
            oed = mopeds.OptimalExperimentalDesign_NLE(m, [variable_list_oed], measurable_variables=["e0_y_c1"], previous_measurements=previous_measurements)
            res = oed.optimize()
            new_exp = res["x_dict"]
            new_exp["e0_x_c1"] = round(new_exp["e0_x_c1"], 4)
            sim.change_independent_variables(new_exp)
            res = sim.simulate()[2]
            y = res["e0_y_c1"].value[0]
            T = res["e0_T"].value[0]
            y_random = round(rng.normal(y, MEAS_STD), 4)
            y_data.append(y_random)
            x_data.append(new_exp["e0_x_c1"])
            T_data.append(T)
            previous_measurements.append(new_exp)


    # y_data = [0.0065, 0.0803, 0.2002, 0.335, 0.48, 0.6333, 0.7675, 0.8649, 0.9272, 0.9825]
    # x_data = [0.01, 0.1189, 0.2278, 0.3367, 0.4456, 0.5544, 0.6633, 0.7722, 0.8811, 0.99]

    def generate_exp_varlist(y_data, x_data):
        """Prepare exp data to parameter estimation"""
        list_exp_varlist = []
        for xx, yy in zip(x_data, y_data):
            vl_i = copy.deepcopy(variable_list)
            vl_i["e0_x_c1"].value = xx
            vl_i["e0_y_c1"].value = yy
            vl_i["e0_y_c1"].variance = MEAS_STD**2
            list_exp_varlist.append(vl_i)
        return list_exp_varlist

    list_exp_varlist = generate_exp_varlist(y_data, x_data)
    list_exp_varlist[0]["e0_greek_lambda_c1"].fixed = False
    list_exp_varlist[0]["e0_greek_lambda_c2"].fixed = False

    pe = mopeds.ParameterEstimationNLE(m, list_exp_varlist)
    res = pe.optimize(direct_optimization=True)
    x_dict = res["x_dict"]

    for param_name in ["e0_greek_lambda_c1", "e0_greek_lambda_c2"]:
        print(f"{param_name} estimated: ", round(x_dict[param_name], 4), "; true value: ", true_params[param_name])
        print("".join(["-"]*50))

    def create_pe_prediction(m):
        x =  np.linspace(0.0, 1.0, 40)
        list_exp_varlist = []
        for xx in x:
            vl_i = copy.deepcopy(variable_list)
            vl_i["e0_x_c1"].value = xx
            vl_i["e0_y_c1"].variance = MEAS_STD**2
            vl_i["e0_T"].variance = 0.03**2
            vl_i["e0_y_c1"].value = 1
            vl_i["e0_T"].value = 1
            vl_i["e0_greek_lambda_c1"].fixed = False
            vl_i["e0_greek_lambda_c2"].fixed = False
            list_exp_varlist.append(vl_i)
        pe = mopeds.ParameterEstimationNLE(m, list_exp_varlist)
        return pe, x

    # Plot after estimation with uncertainty
    if True:
        cov_linearized = pe.calculate_sensitivity_and_fim_fast(x_dict)[2]
        pe_prediction, x_prediction = create_pe_prediction(m)
        jac_prediction = pe_prediction.calculate_sensitivity_and_fim_fast(x_dict)[0]
        y_prediction = pe_prediction.calculate_objective_and_residual(x_dict)["df_all"]
        prediction_std = np.sqrt(np.diag(jac_prediction @ cov_linearized @ jac_prediction.T)).reshape(pe_prediction.array_data.shape, order="F")

        x = y_prediction["e0_y_c1"]
        y = y_prediction["e0_T"]
        plt.plot(x, y)
        plt.plot(x_prediction, y)
        plt.scatter(y_data, T_data, c="r", s=100)
        plt.errorbar(x, y, xerr=prediction_std[:,0],yerr=prediction_std[:,1])
        plt.errorbar(x_prediction, y, yerr=prediction_std[:,1])
        plt.ylabel("e0_T")
        plt.xlabel("e0_x_c1 / e0_y_c1")
        plt.show()




    if False:
        variable_list["e0_y_c1"].variance = MEAS_STD**2
        prediction_grid = {"e0_x_c1": [0, 1, 20]}
        meas_grid = {"e0_x_c1": [0.01, 0.99, 15],}
        analyzer = mopeds.tools.ErrorAnalyzer(variable_list, m, prediction_grid, meas_grid, ["e0_greek_lambda_c1", "e0_greek_lambda_c2"], ["e0_y_c1"])
        analyzer.parameter_covariance_mc(plot=False, num_samples=100)
        analyzer.plot_parameter_covariance_ellipse(normalize_parameters=False)
        plt.show()



