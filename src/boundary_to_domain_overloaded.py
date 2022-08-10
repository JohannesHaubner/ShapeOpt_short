from dolfin import *
from dolfin_adjoint import *

from pyadjoint import Block
from pyadjoint.overloaded_function import overload_function

from src.boundary_to_domain import boundary_to_domain

backend_boundary_to_domain = boundary_to_domain

class B2DBlock(Block):
    def __init__(self, func, W, B, dof_map, **kwargs):
        super(B2DBlock, self).__init__()
        self.kwargs = kwargs
        self.add_dependency(func)
        self.W = W
        self.B = B
        self.dof_map = dof_map

    def __str__(self):
        return 'B2DBlock'

    def evaluate_adj_component(self, inputs, adj_inputs, block_variable, idx, prepared=None):
        # dof map
        
        adj_input = adj_inputs[0]
        adj = adj_input[self.dof_map]

        adjv = Function(self.B)
        adjv.vector()[:] = adj
        return adjv.vector()

    def recompute_component(self, inputs, block_variable, idx, prepared):
        return backend_boundary_to_domain(inputs[0], self.W, self.B, self.dof_map)




boundary_to_domain = overload_function(boundary_to_domain, B2DBlock)

