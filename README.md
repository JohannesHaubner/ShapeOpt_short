# ShapeOpt_short

Short version of implementation https://github.com/JohannesHaubner/ShapeOpt for

>J. Haubner, M. Ulbrich: Advanced Numerical Methods for Shape Optimal Design of Fluid-Structure Interaction Problems.

with an iterative procedure as used in 

>J. Haubner, M. Siebenborn, M. Ulbrich: A continuous perspective on shape optimization via domain transformations, SIAM Journal on Scientific Computing 43 (3), A1997-A2018.

Here:
- the code does not run in parallel
- only a Stokes example is implemented

## Usage/Examples

The Dockerfile (preliminary version) can be used by running:
```
docker build -t shapeopt .
docker run -it shapeopt
```
or
```
docker pull ghcr.io/johanneshaubner/shapeopt
docker run -ti -v ${PWD}:/root/shared -w /root/shared --entrypoint=/bin/bash --rm ghcr.io/johanneshaubner/shapeoopt
```

first run create_mesh.py, then stokes_problem.py

If not ran from Docker image:
Requires a recent master version of dolfin with MeshView support. Requires the changes propsed in https://bitbucket.org/fenics-project/dolfin/issues/1123/assemble-on-mixed-meshview-forms-returns.
