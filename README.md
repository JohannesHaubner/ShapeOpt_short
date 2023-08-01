# ShapeOpt_short

Short version of implementation https://github.com/JohannesHaubner/ShapeOpt 

Here:
- the code does not run in parallel
- only the Stokes example is implemented

## Usage/Examples

The Dockerfile (preliminary version) can be used by running:
```
docker build -t shapeopt .
docker run -it shapeopt
```
or
```
docker pull ghcr.io/johanneshaubner/shapeopt:latest
docker run -it ghcr.io/johanneshaubner/shapeopt:latest
```

first run create_mesh.py, then stokes_problem.py

If not ran from Docker image:
Requires a recent master version of dolfin with MeshView support. Requires the changes propsed in https://bitbucket.org/fenics-project/dolfin/issues/1123/assemble-on-mixed-meshview-forms-returns.
