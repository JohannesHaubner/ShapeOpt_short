from dolfin import *
from dolfin_adjoint import *
from pyadjoint import create_overloaded_object
import numpy as np
from src.project_normal import proj_normal

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


def update_mesh(mesh, mf, design_marker):
    """
    :param mesh: mesh
    :param mf: facet mesh with marked boundary parts
    :param design_marker: marker id for design boundary
    :return:
    """
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
    design_mesh = MeshView.create(bboundaries, design_marker)
    design_mesh = create_overloaded_object(design_mesh)

    dof_map = dof_mapping(mesh, bmesh, design_mesh)

    B = VectorFunctionSpace(design_mesh, "CG", 1)
    nm = FacetNormal(mesh)
    nm = proj_normal(nm, VectorFunctionSpace(mesh, "CG", 1))
    n = Function(VectorFunctionSpace(design_mesh, "CG", 1))
    n.vector()[:] = nm.vector()[dof_map]

    return mesh, mf, design_mesh, n, dof_map