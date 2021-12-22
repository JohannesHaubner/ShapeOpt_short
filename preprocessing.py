from dolfin import *
from dolfin_adjoint import *
import numpy as np


class Preprocessing:
    def __init__(self, FunctionSpace):
        self.V = FunctionSpace

    def dof_to_rho(self, x):
        func = Function(V)
        func.vector()[:] = x
        return func

    def dof_to_rho_chainrule(self, djy, option=1):
        if option == 2:
            return djy[:]
        else:
            return djy.vector()[:]




