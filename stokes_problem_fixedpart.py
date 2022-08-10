from create_mesh import inflow_marker, outflow_marker, wall_marker, obstacle_marker, c_x, c_y, L, H
from dolfin import *
from dolfin_adjoint import *
from pyadjoint.overloaded_type import create_overloaded_object
from scipy import sparse

from src import ipopt_solver as ipopt
from src.preprocessing import Preprocessing
from src.update_mesh import update_mesh

from src.boundary_to_domain_overloaded import boundary_to_domain
import numpy as np
set_log_level(LogLevel.ERROR)

def control_to_deformation(c, reg):
    # We map the scalar valued control c on the design_mesh to a
    # vector valued quantity on the design mesh

    design_mesh = c.function_space().mesh()
    zero = Constant([0] * mesh.geometric_dimension())

    dxb = Measure("dx", domain=design_mesh)

    Bs = FunctionSpace(design_mesh, "CG", 1)
    B = VectorFunctionSpace(design_mesh, "CG", 1)

    # space dependent reg

    r = TrialFunction(Bs)
    psir = TestFunction(Bs)

    a = reg * inner(grad(r), grad(psir)) * dx + inner(r, psir) * dx
    l = reg * psir * dx

    bcs = DirichletBC(Bs, Constant(1.0), 'on_boundary')

    reg_ = Function(Bs)
    solve(a == l, reg_, bcs)

    # control to vector valued quantity on design b

    b = TrialFunction(B)
    psib = TestFunction(B)

    a = inner(reg_ * grad(b), grad(psib)) * dx + inner(b, psib) * dx
    l = inner(c * n, psib) * dxb

    bcs = DirichletBC(B, Constant((0.0, 0.0)), 'on_boundary')

    b = Function(B)
    solve(a == l, b, bcs)

    # We map b to a deformation field on the whole mesh

    W = VectorFunctionSpace(mesh, "CG", 1)
    w = TrialFunction(W)
    psiw = TestFunction(W)

    b = boundary_to_domain(b, W, B, dof_map)

    a = inner(grad(w) + grad(w).T, grad(psiw)) * dx + inner(w, psiw) * dx
    l = inner(b, psiw) * ds #xb  # dS(1) #dX
    bcs = []
    for marker in [inflow_marker, outflow_marker, wall_marker]:
        bcs.append(DirichletBC(W, zero, mf, marker))
    w = Function(W, name="mesh deformation")
    solve(a == l, w, bcs)
    return w

def inner_product(V):
    v = TestFunction(V)
    u = TrialFunction(V)
    M = assemble( v * u * dx)
    return sparse.csc_matrix(M.array())

def forward(mesh):
    # The next step is to set up :eq:`state`. We start by defining the
    # stable Taylor-Hood finite element space.

    V2 = VectorElement("CG", mesh.ufl_cell(), 2)
    S1 = FiniteElement("CG", mesh.ufl_cell(), 1)
    VQ = FunctionSpace(mesh, V2 * S1)

    # Then, we define the test and trial functions, as well as the variational form
    zero = Constant([0] * mesh.geometric_dimension())

    (u, p) = TrialFunctions(VQ)
    (v, q) = TestFunctions(VQ)
    a = inner(grad(u), grad(v)) * dx - div(u) * q * dx - div(v) * p * dx
    l = inner(zero, v) * dx

    # The Dirichlet boundary conditions on :math:`\Gamma` is defined as follows
    zero = Constant([0] * mesh.geometric_dimension())
    #g = Expression(("sin(pi*x[1]/H)", "0"), H = H, degree=2)
    g = Expression(("4/(H*H)*x[1]*(H-x[1])", "0"), H = H, degree=2)
    #x, y = SpatialCoordinate(mesh)
    #g = project(as_vector([4/(H*H)*y*(1-y), 0, 0]), VQ)
    #g1, g2 = g.split(deepcopy=True)

    bc_inlet = DirichletBC(VQ.sub(0), g, mf, inflow_marker)
    bc_obstacle = DirichletBC(VQ.sub(0), zero, mf, obstacle_marker)
    bc_walls = DirichletBC(VQ.sub(0), zero, mf, wall_marker)
    bcs = [bc_inlet, bc_obstacle, bc_walls]

    up = Function(VQ, name="Mixed State Solution")
    solve(a == l, up, bcs=bcs)
    u, p = up.split()

    return u, p

def evaluate_objective_state_independent(w):
    gammaP = 1e5
    etaP = 0.05

    def smoothmax(r, eps=1e-4):
        return conditional(gt(r, eps), r - eps / 2, conditional(lt(r, 0), 0, r ** 2 / (2 * eps)))

    J = assemble(0.5 * gammaP * smoothmax(etaP - det(Identity(mesh.geometric_dimension()) + grad(w))) ** 2 * dx)
    return J

def evaluate_objective(w, u, p):
    J = assemble(0.5*inner(grad(u), grad(u)) * dx)
    return J

def postprocess_mesh(w_opt, par):
    mesh = w_opt.function_space().mesh()

    W = VectorFunctionSpace(mesh, "CG", 1)
    w = TrialFunction(W)
    psiw = TestFunction(W)

    a = inner(grad(w), grad(psiw))*dx(mesh)
    l = par * inner(w_opt, psiw)*dx(mesh)

    zero = Constant([0] * mesh.geometric_dimension())
    bcs = DirichletBC(W, zero, "on_boundary")

    w = Function(W)
    solve(a == l, w, bcs)
    return w


if __name__ == "__main__":
    mesh = Mesh()
    with XDMFFile("mesh.xdmf") as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile("mf.xdmf") as infile:
        infile.read(mvc, "name_to_read")
        mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)

    # We compute the initial volume of the obstacle
    one = Constant(1)
    Vol0 = L * H - assemble(one * dx(domain=mesh))

    file = File("./solution2.pvd")
    file << mesh

    for reg in [1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
        # We create a Boundary-mesh and function space for our control :math:`c`
        mesh, mf, design_mesh, n, dof_map = update_mesh(mesh, mf, obstacle_marker)

        # dof maps between V and B
        B = VectorFunctionSpace(design_mesh, "CG", 1)
        W = VectorFunctionSpace(mesh, "CG", 1)

        set_working_tape(Tape())

        C = FunctionSpace(design_mesh, "CG", 1)
        c = Function(C, name="artificial Control for preconditioning")

        # compute the deformation induced by c
        w = control_to_deformation(c, reg)

        # evaluate state independent parts of objective
        J = evaluate_objective_state_independent(w)

        # We move the domain.
        ALE.move(mesh, w)

        # solve forward equation
        u, p = forward(mesh)

        # evaluate objective
        J += evaluate_objective(w, u, p)

        # We define the reduced functional, where :math:`c` is the design parameter# and use scipy to minimize the objective.
        m = Control(c)
        Jhat = [ReducedFunctional(J, m)]
        scaling_Jhat = [1.0]

        #perturbation = interpolate(Expression(("-A*x[0]"),
        #                                      A=0.005, degree=2), C)
        #results = taylor_to_dict(Jhat, c, perturbation)

        # Define constraints
        (x, y) = SpatialCoordinate(mesh)
        vol_constraint = L*H - assemble( Constant(1.0)*dx(domain=mesh)) - Vol0
        bc_x_constraint = ((L**2 * H / 2 - assemble(x * dx(domain=mesh))) / (Vol0)) - c_x
        bc_y_constraint = ((L * H**2 / 2 - assemble(y * dx(domain=mesh))) / (Vol0)) - c_y

        vc = ReducedFunctional(vol_constraint, m)
        bcx = ReducedFunctional(bc_x_constraint, m)
        bcy = ReducedFunctional(bc_y_constraint, m)

        constraints = [vc, bcx, bcy]
        scaling_constraints = [1.0, 1.0, 1.0]
        bounds = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
        preproc = Preprocessing(C)

        inner_product_matrix = inner_product(C)

        # problem
        problem = ipopt.IPOPTProblem(Jhat, scaling_Jhat, constraints, scaling_constraints, bounds, preproc,
                               inner_product_matrix, reg)
        solver = ipopt.IPOPTSolver(problem)

        parameters = {#'derivative_test': 'first-order',
                      'maximum_iterations': 100, 
                      'tol': 1e-5,
                      'point_perturbation_radius': 0.0}

        c_opt = solver.solve(c.vector()[:])

        # compute optimal deformation
        w_opt = control_to_deformation(preproc.dof_to_rho(c_opt), reg)

        # move mesh
        ALE.move(mesh, w_opt)
        file << w_opt

        # postprocess mesh
        par = 1e0 #1.0
        w_pp = postprocess_mesh(w_opt, par)
        ALE.move(mesh, w_pp)

        file << w_opt
