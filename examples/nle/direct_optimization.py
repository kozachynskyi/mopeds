import mopeds
import casadi as ca


def playground():
    mopeds.set_options(variable_scaling=False)
    vl, m = mopeds.examples.linear_example()

    exp_vl, true_par, _ = mopeds.tools.generate_artificial_data_nle(
        m, vl, [{"u": 1, "v": 3}], perturbate=False
    )

    for var in exp_vl[0].values():
        if isinstance(var, mopeds.VariableParameter):
            var.fixed = False
            true_par[var.name] = var.value[0]

    pe = mopeds.ParameterEstimationNLE(m, exp_vl)
    v = pe.calculate_sensitivity_and_fim(true_par)
    print("PE jac \n", v["jac_full"])

    sim = mopeds.SimulatorNLE(m, vl)
    print("Res sim\n", sim.simulate()[0])

    print("Calculate jac sim\n", sim.calculate_jac()["jac_x_p"])

    jac = sim.function.jacobian()

    d = {"x0": [16, 2, 9], "p": sim._independent_variables}
    print("MAPPING\n", sim.mapping_independent_variables)

    print("Manual eval\n", eval := sim.function_v.call(d))
    breakpoint()
    jacjac = jac.call(d)
    # print("Manual jac\n",jacjac)
    # print( -jacjac["jac_y_p"][:, 2:])#/ jacjac["jac_y_x0"])
    print(ca.inv(jacjac["jac_x_x0"]) @ -jacjac["jac_x_p"])
    print(ca.solve(jacjac["jac_x_x0"], -jacjac["jac_x_p"]))
    breakpoint()

    # print( jacjac["jac_y_x0"] / jacjac["jac_y_p"])


if __name__ == "__main__":
    # playground()
    mopeds.set_options(variable_scaling=False)
    vl, m = mopeds.examples.linear_example()

    for var in vl.values():
        if isinstance(var, mopeds.VariableControl):
            var.fixed = False
        if isinstance(var, mopeds.VariableParameter):
            var.fixed = False
    vl["a"].fixed = True
    # vl["b"].fixed = False
    # vl["y"].variance = 0.1
    # vl["z"].variance = 0.3
    # vl["q"].variance = 0.2


    meas_vars = None
    previous_meas = [
        {"v": 1.5, "u": 3.5},
        {"v": 2.5, "u": 4.5},
        {"v": 3.5, "u": 3.5},
        {"v": 4.5, "u": 3.5},
    ]
    oed = mopeds.OptimalExperimentalDesign_NLE(
        m,
        [vl],
        measurable_variables=meas_vars,
        previous_measurements=previous_meas,
    )
    oed._setup_direct_optimization("OED")
    # oed.generate_jacobian_function_direct()
    oed.solver_settings["ipopt"]["max_iter"] = 100
    # oed.solver_settings["ipopt"]["tol"] = 1e-12
    # oed.solver_settings["ipopt"]["hessian_approximation"] = "exact"
    res1 = oed.optimize(direct_optimization=False)
    res2 = oed.optimize(direct_optimization=True)
    # v1 = {"u": 2.000000014718749, "v": 3.0154388971361943}
    # v2 = {"u": 2.0000000199999075, "v": 3.422127538838273}
    # vv1 = oed.calculate_objective_and_jacobian(v1)
    # vv2= oed.calculate_objective_and_jacobian(v2)
    print(res1["x_dict"])
    print(res2["x_dict"])
    print(oed.calculate_objective_and_jacobian(res1["x_dict"])["f"])
    print(oed.calculate_objective_and_jacobian(res2["x_dict"])["f"])

    breakpoint()
