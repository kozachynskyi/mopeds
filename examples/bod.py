import par_est

VAR_LIST, MODEL, EXP_DATA = par_est.examples.isomerization_model()
# VAR_LIST, MODEL, EXP_DATA = par_est.examples.bod_model()


def parameter_estimation(parameters=None):
    pe = par_est.ParameterEstimationNLE(MODEL, EXP_DATA)
    pe.solver_settings["ipopt"]["linear_solver"] = "ma27"
    if parameters:
        res = pe.parameter_analysis(parameters)
    else:
        res = pe.optimize()
    print(res)


# real_p = [19.143, 0.5311]
real_p = {"theta1": 35.92, "theta2": 0.0708, "theta3": 0.377, "theta4": 0.167}
parameter_estimation()
parameter_estimation(real_p)
