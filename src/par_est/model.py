import casadi as ca
from par_est import (
    VariableAlgebraic,
    VariableState,
    VariableParameter,
    VariableControl,
    VariableConstant,
    Variable,
    VariableList,
)


class Model(object):
    """ Model class is used to get lists of variables used in this model,
    create equations and determine if model is DAE or ODE.
    """

    def __init__(self, variable_list):
        self.varlist_state = VariableList()
        self.varlist_algebraic = VariableList()
        # Includes Parameters and Controls
        self.varlist_independent = VariableList()
        self._varlist_constant = VariableList()
        self.varlist_all = VariableList()
        self.equations_differential = None
        self.equations_algebraic = None
        self.DAE = False

        for var in variable_list.values():
            if isinstance(var, Variable):
                if isinstance(var, VariableState):
                    self.varlist_state.add_variable(VariableState(var.name))
                elif isinstance(var, VariableAlgebraic):
                    self.varlist_algebraic.add_variable(VariableAlgebraic(var.name))
                elif isinstance(var, VariableParameter) or isinstance(
                    var, VariableControl
                ):
                    self.varlist_independent.add_variable(type(var)(var.name))
                elif isinstance(var, VariableConstant):
                    self._varlist_constant.add_variable(
                        VariableConstant(var.name, var.value)
                    )
                else:
                    raise VariableTypeError(var.name)
            else:
                raise VariableTypeError(var.name)

        self.varlist_all.update(self.varlist_state)
        self.varlist_all.update(self.varlist_algebraic)
        self.varlist_all.update(self.varlist_independent)
        self.varlist_all.update(self._varlist_constant)

    def add_equations_differential(self, equations):
        if self.equations_differential is None:
            self.equations_differential = ca.vcat(equations)
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)

    def add_equations_algebraic(self, equations):
        if self.equations_algebraic is None:
            self.equations_algebraic = ca.vcat(equations)
            self.DAE = True
        else:
            # Adding additional equations is not implemented
            raise (NotImplementedError)


class VariableTypeError(Exception):
    def __init__(self, name):
        message = f"Not a supported par_est_casadi variable class! Wrong variable with name: {name}"
        super().__init__(message)
