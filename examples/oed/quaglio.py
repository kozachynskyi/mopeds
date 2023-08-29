import numpy as np
import copy
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from tools import unfix_parameters
import casadi as ca

import par_est

class CriteriaD_espie1989(par_est.OED_objective):
    def eval(self, args):
        jac = args[0]
        jac_scaled = args[0] * self._parameter_scaling
        obj = np.linalg.det(jac_scaled.T @ jac_scaled)
        return obj, jac

class CriteriaD_asprey2002(par_est.OED_objective):
    def eval(self, args):
        jac = args[0]
        jac_scaled = args[0] * self._parameter_scaling
        obj = np.linalg.det(jac_scaled.T @ jac_scaled)
        return obj, jac

class CriteriaD_quaglio2018(par_est.OED_objective):
    def eval(self, args):
        jac = args[0]
        jac_scaled = args[0]# * self._parameter_scaling
        obj = np.linalg.det(jac_scaled.T @ jac_scaled)
        return obj, jac

class OED_magnusson2015(par_est.OptimalExperimentalDesign):
    def _objective_A(self):
        """A criteria"""
        jac = self.jacobian_mx
        jac_scaled = jac
        obj = ca.trace(ca.inv(jac_scaled.T @ jac_scaled))

        return obj, jac


def espie1989():
    varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True)

    par_initial = {
        "theta1": 0.3,
        "theta2": 0.25,
        "theta3": 0.56,
        "theta4": 0.02,
    }

    for par_name, par_value in par_initial.items():
        varlist[par_name].value = par_value

    # varlist["x1"].variance = 0.2**2
    # varlist["x2"].variance = 0.2**2

    varlist["x1"].variance = 1
    varlist["x2"].variance = 1

    varlist["x1"].value = 1.0
    varlist["x2"].value = 0.01

    varlist["u1"].value = 0.2
    varlist["u1"].fixed = True
    varlist["u2"].ignore_plotting = False

    varlist["u2"].lower_bound = 5
    varlist["u2"].upper_bound = 35

    mode = "constant"
    mode = "determinant"
    # mode = "eigenvalue"
    # mode = "optimize_solution"
    # mode = "optimize_dict"

    time_grid = np.linspace(0, 72, int((72/0.75 + 1)))
    oed_settings = par_est.OEDsettings(num_control_switches=1)

    if mode == "constant":
        # Constant input, Figure 2 and 3
        varlist["u2"].value = 35
        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_espie1989)["f"]
    elif mode == "determinant":
        varlist["u2"].value = 35
        varlist["u2"].expand_horizon([37.5, 50], [5, 35])
        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_espie1989)["f"]
    elif mode == "eigenvalue":
        varlist["u2"].value = 5
        varlist["u2"].expand_horizon([12, 23.25], [12, 35])
    elif mode == "optimize_solution":
        varlist["u2"].value = 35
        varlist["u2"].expand_horizon([0.75, 1.5, 2.25, 3.0, 3.75, 4.5, 5.25, 6.0, 6.75, 7.5, 8.25, 9.0, 9.75, 10.5, 11.25, 12.0, 12.75, 13.5, 14.25, 15.0, 15.75, 16.5, 17.25, 18.0, 18.75, 19.5, 20.25, 21.0, 21.75, 22.5, 23.25, 24.0, 24.75, 25.5, 26.25, 27.0, 27.75, 28.5, 29.25, 30.0, 30.75, 31.5, 32.25, 33.0, 33.75, 34.5, 35.25, 36.0, 36.75, 37.5, 38.25, 39.0, 39.75, 40.5, 41.25, 42.0, 42.75, 43.5, 44.25, 45.0, 45.75, 46.5, 47.25, 48.0, 48.75, 49.5, 50.25, 51.0, 51.75, 52.5, 53.25, 54.0, 54.75, 55.5, 56.25, 57.0, 57.75, 58.5, 59.25, 60.0, 60.75, 61.5, 62.25, 63.0, 63.75, 64.5, 65.25, 66.0, 66.75, 67.5, 68.25, 69.0, 69.75, 70.5, 71.25], [35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 34.41,  5.  ,  5.  ,  5.  ,  5.  ,  5.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  , 35.  ])
        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, "D")["f"]
    elif mode == "optimize_dict":
        solution = {'u2_t0': 34.99988643472994, 'u2_t1': 34.99968940498816, 'u2_t2': 34.99961631831301, 'u2_t3': 34.999624642771735, 'u2_t4': 34.99966619272635, 'u2_t5': 34.99971097586224, 'u2_t6': 34.999745491330195, 'u2_t7': 34.999766106679175, 'u2_t8': 34.99977354109989, 'u2_t9': 34.999769528072406, 'u2_t10': 34.99975612983923, 'u2_t11': 34.99973541412478, 'u2_t12': 34.99970942891819, 'u2_t13': 34.999680672585505, 'u2_t14': 34.999651267461395, 'u2_t15': 34.99962308695431, 'u2_t16': 34.999597750484185, 'u2_t17': 34.99957655155194, 'u2_t18': 34.99956014858205, 'u2_t19': 34.999548970866044, 'u2_t20': 34.99954279351537, 'u2_t21': 34.99954158963804, 'u2_t22': 34.99954483384351, 'u2_t23': 34.999552030280235, 'u2_t24': 34.99956262051535, 'u2_t25': 34.999575997238665, 'u2_t26': 34.999591647137834, 'u2_t27': 34.99960910583337, 'u2_t28': 34.99962790252663, 'u2_t29': 34.99964771056676, 'u2_t30': 34.999668192143595, 'u2_t31': 34.99968910894469, 'u2_t32': 34.999710214808054, 'u2_t33': 34.99973134657434, 'u2_t34': 34.99975232668174, 'u2_t35': 34.99977302209309, 'u2_t36': 34.99979328450823, 'u2_t37': 34.99981298221178, 'u2_t38': 34.99983195553349, 'u2_t39': 34.99985005580783, 'u2_t40': 34.99986712832524, 'u2_t41': 34.9998830381301, 'u2_t42': 34.999897688233375, 'u2_t43': 34.99991102228453, 'u2_t44': 34.99992301328389, 'u2_t45': 34.99993368050121, 'u2_t46': 34.99994324337886, 'u2_t47': 34.99995203851732, 'u2_t48': 34.999960365794266, 'u2_t49': 34.999968177119584, 'u2_t50': 34.999974883294534, 'u2_t51': 34.99997839284283, 'u2_t52': 34.99997554112406, 'u2_t53': 34.99994920573331, 'u2_t54': 34.99993792303316, 'u2_t55': 34.99993553743843, 'u2_t56': 34.999933945833625, 'u2_t57': 34.999933661929134, 'u2_t58': 34.99993465557629, 'u2_t59': 34.99993584956759, 'u2_t60': 34.999936813379044, 'u2_t61': 34.99993750563603, 'u2_t62': 34.99993809019119, 'u2_t63': 34.99993878387512, 'u2_t64': 34.999939668023636, 'u2_t65': 34.99994079802956, 'u2_t66': 34.99994224596091, 'u2_t67': 34.99994397555426, 'u2_t68': 34.99994588065947, 'u2_t69': 34.999947925768446, 'u2_t70': 34.999950043538604, 'u2_t71': 34.99995211103686, 'u2_t72': 34.999954073006194, 'u2_t73': 34.999955835576046, 'u2_t74': 34.999957300853815, 'u2_t75': 34.99995839376943, 'u2_t76': 34.99995902174444, 'u2_t77': 34.99995909782687, 'u2_t78': 34.99995866248742, 'u2_t79': 34.99995814183527, 'u2_t80': 34.999956747992385, 'u2_t81': 34.999928834963, 'u2_t82': 34.99986840445507, 'u2_t83': 34.40587381794541, 'u2_t84': 5.000031067126143, 'u2_t85': 5.000021163004472, 'u2_t86': 5.000017264379194, 'u2_t87': 5.000013019313238, 'u2_t88': 5.000017953746731, 'u2_t89': 34.99996246439409, 'u2_t90': 34.99996581467536, 'u2_t91': 34.99996727250451, 'u2_t92': 34.9999621589868, 'u2_t93': 34.99996252301505, 'u2_t94': 34.99994455108265, 'u2_t95': 34.99983025127333}
        varlist_oed = unfix_parameters(varlist)
        varlist_oed["u2"].fixed = False
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid, oed_settings)
        obj = oed.calculate_objective_and_jacobian(solution, "D")["f"]
    else:
        raise NotImplementedError
    print(mode)
    print(obj / 1e13)
    breakpoint()

    if False:
        sim = par_est.Simulator(m_monod, time_grid, varlist)
        sim.generate_exp_data().plot()

    unfix_variables = ["theta1", "theta2", "theta3","theta4", "u2"]
    varlist_oed = copy.deepcopy(varlist)

    for par_name in unfix_variables:
        varlist_oed[par_name].fixed = False

    oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid, oed_settings)
    oed.solver_settings["ipopt"]["max_iter"] = 20
    oed.solver_settings["ipopt"]["linear_solver"] = "ma57"
    oed.optimize(1, "D")

    print(oed.calculate_objective_and_jacobian(solution, "D")["f"])
    a = oed.generate_experimental_data(solution, par_initial).plot()
    breakpoint()

def asprey2002():
    varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True)

    par_initial = {
        "theta1": 0.5,
        "theta2": 0.5,
        "theta3": 0.5,
        "theta4": 0.5,
    }

    # varlist["x1"].variance = 0.2**2
    # varlist["x2"].variance = 0.2**2

    varlist["x1"].variance = 1
    varlist["x2"].variance = 1

    varlist["x1"].lower_bound = 1.0
    varlist["x1"].upper_bound = 10.0

    varlist["x2"].value = 0.1
    varlist["x2"].fixed = True

    varlist["u1"].lower_bound = 0.05
    varlist["u1"].upper_bound = 0.2
    varlist["u1"].ignore_plotting = False

    varlist["u2"].lower_bound = 5
    varlist["u2"].upper_bound = 35
    varlist["u2"].ignore_plotting = False

    for par_name, par_value in par_initial.items():
        varlist[par_name].value = par_value

    mode = "initial"
    # mode = "determinant"
    # mode = "eigenvalue"
    # mode = "optimize_solution"
    # mode = "optimize_dict"


    oed_settings = par_est.OEDsettings(num_control_switches=1)

    if mode == "initial":
        time_grid = [0,2,4,6,8,10,12,14,16,18,20]
        varlist["x1"].value = 5.5

        varlist["u1"].value = 0.12

        varlist["u2"].value = 15
        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_asprey2002)

    elif mode == "determinant":
        time_grid = [0, 21.2, 22.2, 23.2, 24.2, 25.2, 26.2, 27.2, 28.2, 29.2, 30.2]
        varlist["x1"].value = 8.53

        varlist["u1"].value = 0.2
        varlist["u1"].expand_horizon([5.4], [0.05])

        varlist["u2"].value = 35
        varlist["u2"].expand_horizon([20, 25.2], [22.8, 15.])

        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_asprey2002)
    else:
        raise NotImplementedError
    print(mode)
    print(obj["f"])
    breakpoint()

    if True:
        time_grid = np.linspace(0, 32, 1000)
        sim = par_est.Simulator(m_monod, time_grid, varlist)
        sim.generate_exp_data().plot()

    unfix_variables = ["theta1", "theta2", "theta3","theta4", "u2"]
    varlist_oed = copy.deepcopy(varlist)

    for par_name in unfix_variables:
        varlist_oed[par_name].fixed = False

    oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid, oed_settings)
    oed.solver_settings["ipopt"]["max_iter"] = 20
    oed.solver_settings["ipopt"]["linear_solver"] = "ma57"
    oed.optimize(1, "D")

    print(oed.calculate_objective_and_jacobian(solution, "D")["f"])
    a = oed.generate_experimental_data(solution, par_initial).plot()
    breakpoint()

def hoang2013():
    varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True, u1_piecewise_linear=True)

    # par_initial = {
    #     "theta1": 0.5,
    #     "theta2": 0.5,
    #     "theta3": 0.5,
    #     "theta4": 0.5,
    # }
    par_initial = {'theta1': 0.5,
                 'theta2': 0.06,
                 'theta3': 0.5,
                 'theta4': 0.25}

    # varlist["x1"].variance = 0.2**2
    # varlist["x2"].variance = 0.2**2

    # varlist["x1"].variance = 1
    # varlist["x2"].variance = 1

    varlist["x1"].lower_bound = 1.0
    varlist["x1"].upper_bound = 10.0

    varlist["x2"].value = 0.1
    varlist["x2"].fixed = True

    varlist["u1"].lower_bound = 0.05
    varlist["u1"].upper_bound = 0.5
    varlist["u1"].ignore_plotting = False

    varlist["u2"].lower_bound = 0.2
    varlist["u2"].upper_bound = 35
    varlist["u2"].ignore_plotting = False

    for par_name, par_value in par_initial.items():
        varlist[par_name].value = par_value

    # mode = "initial"
    mode = "determinant1"
    # mode = "eigenvalue"
    # mode = "optimize_solution"
    # mode = "optimize_dict"
    # mode = "pe"


    oed_settings = par_est.OEDsettings(num_control_switches=1)

    if mode == "initial":
        time_grid = [0,3,6,9,12]
        varlist["x1"].value = 2
        varlist["x2"].value = 0.1

        varlist["u1_dot"].value = 0.1
        varlist["u1"].value = 0

        varlist["u2"].value = 15

        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_asprey2002)

    elif mode == "determinant1":
        time_grid = [0,3,6,9,12]
        varlist["x1"].value = 2
        varlist["x2"].value = 0.1

        varlist["u1_dot"].value = 0.1
        varlist["u1"].value = -0.005
        varlist["u1"].expand_horizon([1,2,3,6,11], [0.0623683249324275, -0.0430126378844328, 0.000762496762496763, 0.000459587089724076, -0.0167749287749288])

        varlist["u2"].value = 15
        # varlist["u2"].expand_horizon([11], [20])

        varlist_oed = unfix_parameters(varlist)
        oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_asprey2002)
    elif mode == "pe":
        varlist["x1"].value = 2
        varlist["x2"].value = 0.1

        varlist["u1_dot"].value = 0.1
        varlist["u1"].value = 0

        varlist["u2"].value = 15

        result_pe = {'theta1': 0.5,
                     'theta2': 0.06,
                     'theta3': 0.5,
                     'theta4': 0.25}

        x1_data = [
                    [0,2],
                    [0.630508474576271,2.009977827051],
                    [2.61016949152542,2.06984478935698],
                    [4.35932203389831,2.09977827050998],
                    [5.9728813559322,2.12305986696231],
                    [8.4135593220339,2.13968957871397],
                    [10.4474576271186,2.12638580931264],
                    [11.9796610169492,2.12971175166297]
                ]
        x2_data = [
                [0,0.1], 
                [0.942372881355932,0.157427937915743],
                [2.09491525423729,0.144124168514412],
                [3.5864406779661,0.144124168514412],
                [4.98305084745763,0.137472283813747],
                [6.52881355932203,0.127494456762749],
                [7.84406779661017,0.130820399113082],
                [9.66101694915254,0.14079822616408],
                [11.9932203389831,0.137472283813747]
                ]

        varlist["theta2"].fixed = False
        varlist["theta4"].fixed = False

        for var, data in zip([varlist["x1"], varlist["x2"]], [x1_data, x2_data]):
            data_numpy = np.array(data)
            var.set_dataframe_from_value_and_time(data_numpy[:,1], data_numpy[:,0])

        pe = par_est.ParameterEstimation(m_monod, [varlist])
        print(pe.calculate_objective_and_residual(result_pe))
        # a = pe.optimize()
        breakpoint()

    else:
        raise NotImplementedError
    print(mode)
    print(obj["f"])

    if True:
        time_grid = np.linspace(0, 12, 100)
        sim = par_est.Simulator(m_monod, time_grid, varlist)
        a = sim.generate_exp_data()
        breakpoint()
        a.plot()

    unfix_variables = ["theta1", "theta2", "theta3","theta4", "u2"]
    varlist_oed = copy.deepcopy(varlist)

    for par_name in unfix_variables:
        varlist_oed[par_name].fixed = False

    oed = par_est.OptimalExperimentalDesign(m_monod, [varlist_oed], time_grid, oed_settings)
    oed.solver_settings["ipopt"]["max_iter"] = 20
    oed.solver_settings["ipopt"]["linear_solver"] = "ma57"
    oed.optimize(1, "D")

    print(oed.calculate_objective_and_jacobian(solution, "D")["f"])
    a = oed.generate_experimental_data(solution, par_initial).plot()
    breakpoint()

def magnusson2015():
    varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True, ode=True)

    par_initial = {
        "theta1": 0.1,
        "theta2": 0.1,
        "theta3": 0.1,
        "theta4": 0.1,
    }

    varlist["x1"].variance = 5**2
    varlist["x2"].variance = 5**2

    # varlist["x1"].variance = 1
    # varlist["x2"].variance = 1

    varlist["x1"].value = 7

    varlist["x2"].value = 0

    varlist["u1"].lower_bound = 0.05
    varlist["u1"].upper_bound = 0.5
    varlist["u1"].ignore_plotting = False

    varlist["u2"].lower_bound = 0.2
    varlist["u2"].upper_bound = 35
    varlist["u2"].ignore_plotting = False

    for par_name, par_value in par_initial.items():
        varlist[par_name].value = par_value

    mode = "trace_modelica"
    # mode = "trace_matlab"
    # mode = "trace_gproms"

    time_grid = np.linspace(0, 20, 6)
    varlist["x1"].value = 7
    varlist["x2"].value = 0

    if mode == "trace_modelica":
        varlist["u1"].value = 0.05
        varlist["u1"].expand_horizon([4,8,12], [0.08, 0.1, 0.05])

        varlist["u2"].value = 5
        varlist["u2"].expand_horizon([4,12,16], [35, 5, 20])

        varlist_oed = unfix_parameters(varlist)
        oed = OED_magnusson2015(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, "A")

    elif mode == "trace_matlab":
        varlist["u1"].value = 0.05
        varlist["u1"].expand_horizon([4,8,12,16], [0.08, 0.12, 0.05, 0.14])

        varlist["u2"].value = 15

        varlist_oed = unfix_parameters(varlist)
        oed = OED_magnusson2015(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, "A")

    elif mode == "trace_gproms":
        varlist["u1"].value = 0.05
        varlist["u1"].expand_horizon([8,12,16], [0.06, 0.18, 0.05])

        varlist["u2"].value = 5
        varlist["u2"].expand_horizon([4,8,12,16], [23, 28,30,5])

        varlist_oed = unfix_parameters(varlist)
        oed = OED_magnusson2015(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, "A")
    else:
        raise NotImplementedError
    print(mode)
    print(obj["f"])

    if True:
        time_grid = np.linspace(0, 20, 100)
        sim = par_est.Simulator(m_monod, time_grid, varlist)
        a = sim.generate_exp_data()
        a.plot()

def quaglio2018():
    varlist, m_monod, exp_data = par_est.examples.yeast_growth("monod", piecewise=True)

    par_initial = {
        "theta1": 0.531,
        "theta2": 7.854,
        "theta3": 0.474,
        "theta4": 0.019,
    }

    varlist["x1"].variance = 0.01
    varlist["x2"].variance = 0.05

    varlist["u1"].ignore_plotting = False

    varlist["u2"].ignore_plotting = False

    for par_name, par_value in par_initial.items():
        varlist[par_name].value = par_value

    mode = "first_d"

    time_grid = np.array([0, 5, 10, 15, 20])

    if mode == "first_d":
        varlist["x1"].value = 5
        varlist["x2"].value = 0.01

        varlist["u1"].value = 0.2
        varlist["u2"].value = 35

        varlist_oed = unfix_parameters(varlist)
        oed = OED_magnusson2015(m_monod, [varlist_oed], time_grid)
        obj = oed.calculate_objective_and_jacobian({}, CriteriaD_quaglio2018)

    else:
        raise NotImplementedError
    print(mode)
    print(obj["f"]/1e16)

    if True:
        time_grid = np.linspace(0, 20, 100)
        sim = par_est.Simulator(m_monod, time_grid, varlist)
        a = sim.generate_exp_data()
        a.plot()


if __name__ == "__main__":
    # espie1989()
    # asprey2002()
    # hoang2013()
    # magnusson2015()
    quaglio2018()
