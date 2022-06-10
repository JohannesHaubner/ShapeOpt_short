# ShapeOpt_short

Requires a recent master version of dolfin with MeshView support.
(If it does not work properly, it might require the changes propsed in https://bitbucket.org/fenics-project/dolfin/issues/1123/assemble-on-mixed-meshview-forms-returns )

The Dockerfile (preliminary version) can be used by running:
```
docker build -t shapeopt .
docker run -it shapeopt
```

first run create_mesh.py, then stokes_problem.py
