import numpy as np
import casadi as ca
import copy
import par_est


class SimSaskia(par_est.SimulatorNLE):
    def model_analysis(self):
        J_mx = ca.jacobian(self.model.equations_algebraic, self.model.varlist_algebraic.get_casadi_variables())
        J = ca.Function("n", [ca.vcat([self.model.varlist_algebraic.get_casadi_variables(), self.model.varlist_independent.get_casadi_variables()])], [J_mx])
        J0 = J(ca.vcat([self.call_arg["x0"], self.call_arg["p"]]))

        if not np.all(np.isnan(J0)):
            print("Jacobian is zero at provided initial values")
            index_nans = np.argwhere(np.isnan(J0))

            var_names = list(self.model.varlist_algebraic.keys())
            equation_formulas = self.model.equations_algebraic

            previous_equation_index = None
            for equation_index, variable_index in index_nans:
                print(f"Variable {var_names[variable_index]}")
                if previous_equation_index is None or not previous_equation_index == equation_index:
                    print(f"{equation_formulas[equation_index]}")
                previous_equation_index = equation_index

        nb, rowperm, colperm, rowblock, colblock, coarse_rowblock, coarse_colblock = J0.sparsity().btf()
        breakpoint()

        solution = []
        solution_dict = {}
        solved_mx = []
        solver_settings = copy.deepcopy(self.solver_settings)
        try:
            solver_settings.pop("constraints")
        except:
            pass

        variables_algebraic_list = list(self.model.varlist_algebraic.values())


        for block_index in range(0, nb):
            if block_index == 0:
                parameters_mx = self.model.varlist_independent.get_casadi_variables()
                parameters_values = self.call_arg["p"]
            else:
                parameters_mx = ca.vertcat(parameters_mx, block_variables)
                parameters_values = ca.vertcat(self.call_arg["p"], *solution)

            equation_indexes = rowperm[rowblock[block_index]: rowblock[block_index+1]]
            var_indexes = colperm[colblock[block_index]: colblock[block_index+1]]
            block_equations = self.model.equations_algebraic.get(False, equation_indexes)
            vars_in_equations = ca.symvar(block_equations)
            # block_variables = ca.vcat(ca.symvar(self.model.varlist_algebraic.get_casadi_variables().get(False, var_indexes, [0])))
            # block_variables = self.model.varlist_algebraic.get_casadi_variables().get(False, var_indexes, [0])


            block_variables = []
            guess = []
            for var_index in var_indexes:
                if variables_algebraic_list[var_index].casadi_var == "e0_D_L_j1_ii3_jj1":
                    breakpoint()
                block_variables.append(variables_algebraic_list[var_index].casadi_var)
                guess.append(self.call_arg["x0"][var_index])

            guess = ca.DM(guess)
            block_variables = ca.vcat(block_variables)



            # for var in vars_in_equations:
            #     if var.name() in self.mapping_independent_variables:

            f = ca.rootfinder("block_solver", "nlpsol", {"x": block_variables, "p": parameters_mx, "g": block_equations}, solver_settings)

            f_sol = f.call({"x0": guess, "p": parameters_values})["x"]
            solution.append(f_sol)

            for index in range(f_sol.shape[0]):
                solution_dict[block_variables[index].name()] = float(f_sol[index])
            # print(block_variables)
            print(solution_dict.keys())
            # print(solution)


            if block_index == 0:
                solved_mx = ca.vcat(vars_in_equations)
            else:
                solved_mx = ca.vcat([solved_mx, *vars_in_equations])

            # print(block_equations, block_variables, vars_in_equations)
            # breakpoint()

VAR_LIST, MODEL, EXP_DATA = par_est.examples.isomerization_model()
sim = SimSaskia(MODEL, VAR_LIST)
sim.model_analysis()
