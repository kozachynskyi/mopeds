from casadi import *
from numpy import *
from pylab import *
import par_est
from par_est.simulation import SimulatorNLE
from par_est.optimization import ParameterEstimationNLE
import copy

eps = SX.sym("eps")
mu = SX.sym("mu")
alpha = SX.sym("alpha")
k = SX.sym("k")
sigma = SX.sym("sigma")
params = [eps, mu, alpha, k, sigma]

#! Variables
a = SX.sym("a")
gamma = SX.sym("gamma")

#! Equations
res0 = mu * a + 1.0 / 2 * k * a * sin(gamma)
res1 = -sigma * a + 3.0 / 4 * alpha * a ** 3 + k * a * cos(gamma)

#! Numerical values
sigma_ = 0.1
alpha_ = 0.1
k_ = 0.2
params_ = [0.1, 0.1, alpha_, k_, sigma_]

#! We create a Function instance
f = Function("f", [vertcat(a, gamma), vertcat(*params)], [vertcat(res0, res1)])
opts = {}
opts["strategy"] = "linesearch"
opts["abstol"] = 1e-14

# $ Require $a > 0$ and $\gamma < 0$
opts["constraints"] = [2, -2]
s = rootfinder("s", "kinsol", f, opts)
x_ = s([1,-1], params_)
print("Solution = ", x_)

variable_list = par_est.VariableList()
variable_list.add_variable(par_est.VariableState("a", 1.0))
variable_list.add_variable(par_est.VariableState("gamma", -1.0))

variable_list.add_variable(par_est.VariableParameter("eps", 0.1, 0.05 , 0.2))
variable_list.add_variable(par_est.VariableParameter("mu", 0.1, 0.05 , 0.2))
variable_list.add_variable(par_est.VariableParameter("sigma", 0.1, 0.05 , 0.2))
variable_list.add_variable(par_est.VariableParameter("alpha", 0.1, 0.05 , 0.2))
variable_list.add_variable(par_est.VariableParameter("k", 0.2, 0.15 , 0.25))

for var in variable_list.values():
    var.guess = var.lower_bound

m = par_est.Model(variable_list)

res0 = m.varlist_all["mu"].casadi_var * m.varlist_all["a"].casadi_var + 1.0 / 2 * m.varlist_all["k"].casadi_var * m.varlist_all["a"].casadi_var * sin( m.varlist_all["gamma"].casadi_var)
res1 = - m.varlist_all["sigma"].casadi_var * m.varlist_all["a"].casadi_var+ 3.0 / 4 * m.varlist_all["alpha"].casadi_var * m.varlist_all["a"].casadi_var ** 3 + m.varlist_all["k"].casadi_var* m.varlist_all["a"].casadi_var* cos(m.varlist_all["gamma"].casadi_var)

m.add_equations_differential([res0, res1])

var_fixed = copy.deepcopy(variable_list)
for var in var_fixed.values():
    var.fixed = True


sim_fixed = SimulatorNLE(m, var_fixed)
var_list_exp = sim_fixed.generate_exp_data()

variable_list_optimizer = copy.deepcopy(variable_list)
for key, var in var_list_exp.items():
    variable_list_optimizer[key] = var


for var in variable_list_optimizer.values():
    var.fixed = True

variable_list_optimizer["k"].fixed = False

print(var_list_exp["a"].value.value)

pe = ParameterEstimationNLE(m, [variable_list_optimizer, variable_list_optimizer])

res = pe.optimize(False)
print(res)

# res = sim.generate_exp_data()
