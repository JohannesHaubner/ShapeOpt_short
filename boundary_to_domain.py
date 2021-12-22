from dolfin import *
from dolfin_adjoint import *

import numpy as np


def boundary_to_domain(b, W, B, dof_map):
    ##
    ## TODO check if b is in the correct function space
    # dof map
    w = Function(W)
    w.vector()[dof_map] = b.vector().get_local()
    return w