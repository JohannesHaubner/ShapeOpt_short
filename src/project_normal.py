from dolfin import *

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