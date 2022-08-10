from create_mesh import inflow_marker, outflow_marker, wall_marker, obstacle_marker, c_x, c_y, L, H
from dolfin import *
from dolfin_adjoint import *
from pyadjoint.overloaded_type import create_overloaded_object
from scipy import sparse

from src import ipopt_solver as ipopt
from src.preprocessing import Preprocessing

from src.boundary_to_domain_overloaded import boundary_to_domain
import numpy as np
set_log_level(LogLevel.ERROR)

def control_to_deformation(c, reg):
    # We map the scalar valued control c on the design_mesh to a
    # vector valued quantity on the design mesh
    # control to vector valued quantity on design b

    design_mesh = c.function_space().mesh()
    zero = Constant([0] * mesh.geometric_dimension())

    dxb = Measure("dx", domain=design_mesh)

    B = VectorFunctionSpace(design_mesh, "CG", 1)
    b = TrialFunction(B)
    psib = TestFunction(B)

    a = reg * inner(grad(b), grad(psib)) * dx + inner(b, psib) * dx
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

def dof_mapping_scalar(mesh, bmesh, dmesh):
    C = FunctionSpace(dmesh, "CG", 1)
    Vb = FunctionSpace(bmesh, "CG", 1)
    W = FunctionSpace(mesh, "CG", 1)

    # Build cell mapping between sub and parent mesh
    cell_map = dmesh.topology().mapping()[bmesh.id()].cell_map()

    # Get cell dofmaps
    dofmap = C.dofmap()
    dofmap_full = Vb.dofmap()

    f = Function(C)

    # Transfer dofs
    GValues = np.zeros(np.size(f.vector().get_local()))
    for c in cells(dmesh):
        GValues[dofmap.cell_dofs(c.index())] = dofmap_full.cell_dofs(cell_map[c.index()])

    # from boundary mesh to mesh
    mapb = bmesh.entity_map(0)
    d2v = dof_to_vertex_map(Vb)
    v2d = vertex_to_dof_map(W)

    f = Function(Vb)

    Vb_to_W_map = np.zeros(np.size(f.vector().get_local()))
    for i in range(np.size(f.vector().get_local())):
        GVertID = Vertex(bmesh, d2v[i]).index()  # Local Vertex ID for given dof on boundary mesh
        PVertID = mapb[GVertID]  # Local Vertex ID of parent mesh
        PDof = v2d[PVertID]
        Vb_to_W_map[i] = int(PDof)

    dof_map = [Vb_to_W_map[int(c)] for c in GValues]
    return dof_map

def dof_mapping(mesh, bmesh, dmesh):
    dof_map_scalar = dof_mapping_scalar(mesh, bmesh, dmesh)
    W = VectorFunctionSpace(mesh, "CG", 1)
    D = VectorFunctionSpace(dmesh, "CG", 1)
    Ds = FunctionSpace(dmesh, "CG", 1)

    w = Function(W)
    d = Function(D)
    w.vector()[:] = range(np.size(w.vector()[:]))
    #d.vector()[:] = range(np.size(d.vector()[:]))
    w_is = w.split(deepcopy=True)
    d_is = d.split(deepcopy=True)
    ws = []
    for w_i in w_is:
        d = Function(Ds)
        d.vector()[:] = w_i.vector()[dof_map_scalar]
        ws.append(d)
    split_to_vec = FunctionAssigner(D, [wi.function_space() for wi in ws])
    wn = Function(D)
    split_to_vec.assign(wn, ws)
    return wn.vector()[:]


def update_mesh(mesh, mf):
    mesh_new = Mesh(mesh)
    mvc = MeshValueCollection("size_t", mesh_new, 1)
    mf_new = cpp.mesh.MeshFunctionSizet(mesh_new, mvc)
    mf_new.set_values(mf.array())
    mesh = mesh_new
    mf = mf_new
    bmesh = BoundaryMesh(mesh, "exterior")
    bmvc = MeshValueCollection("size_t", bmesh, 1)
    bboundaries = cpp.mesh.MeshFunctionSizet(bmesh, bmvc)
    dofs = bmesh.entity_map(1)
    bnum = bmesh.num_vertices()
    for i in range(bnum):
        bboundaries.set_value(i, mf[dofs[i]])
    design_mesh = MeshView.create(bboundaries, obstacle_marker)
    design_mesh = create_overloaded_object(design_mesh)

    dof_map = dof_mapping(mesh, bmesh, design_mesh)

    B = VectorFunctionSpace(design_mesh, "CG", 1)
    nm = FacetNormal(mesh)
    nm = proj_normal(nm, VectorFunctionSpace(mesh, "CG", 1))
    n = Function(VectorFunctionSpace(design_mesh, "CG", 1))
    n.vector()[:] = nm.vector()[dof_map]

    return mesh, mf, design_mesh, n, dof_map

def proj_normal(n, V):
    # Define variational problem for projection
    w = TestFunction(V)
    Pn = TrialFunction(V)
    a = inner(w, Pn)*ds
    L = dot(w, n)*ds
    A = assemble(a, keep_diagonal=True)
    A.ident_zeros()
    b = assemble(L)
    np = Function(V)
    solve(A, np.vector(), b)
    return np


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
        mesh, mf, design_mesh, n, dof_map = update_mesh(mesh, mf)

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

        #parameters = {#'derivative_test': 'first-order',
        #              'maximum_iterations': 100,
        #              'tol': 1e-5,
        #              'point_perturbation_radius': 0.0}

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
