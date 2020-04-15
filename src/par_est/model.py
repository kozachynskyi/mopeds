import casadi as ca
from par_est import VariableAlgebraic, VariableState, Variable, VariableList


class Model(object):

    def __init__(self, variable_list):
        self.varlist_state = VariableList()
        self.varlist_algebraic = VariableList()
        self.varlist_independent = VariableList()
        self.varlist_all = VariableList()
        self.equations_differential = None
        self.equations_algebraic = None

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    self.varlist_state.add_variable(VariableState(var.name))
                elif isinstance(var, VariableAlgebraic):
                    self.varlist_algebraic.add_variable(VariableAlgebraic(var.name))
                else:
                    self.varlist_independent.add_variable(type(var)(var.name))
            else:
                raise (TypeError)

        self.varlist_all.update(self.varlist_state)
        self.varlist_all.update(self.varlist_algebraic)
        self.varlist_all.update(self.varlist_independent)

    def add_equations_differential(self, equations):
        if self.equations_differential is None:
            self.equations_differential = ca.vcat(equations)
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    def add_equations_algebraic(self, equations):
        if self.equations_algebraic is None:
            self.equations_algebraic = ca.vcat(equations)
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)
