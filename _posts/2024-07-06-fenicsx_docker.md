---
layout: post
title: 'Setting up Fenicsx on Windows OS using Docker'
date: 2024-07-06
permalink: /posts/2024/07/fenicsx_docker/
tags:
---

<br><br>

## Installing Docker

1. Install Docker Desktop from <https://www.docker.com/products/docker-desktop/>
2. Once installed, you would have to manually launch Docker Desktop app to start the docker engine.
3. The Docker engine needs to be running whenever you want to use Docker.

## Installing Fenicsx and other requirements

1. 
```shell
docker run --init -ti -p 8888:8888 dolfinx/lab:stable
```
    This will run a Fenicsx provided image which will launch a jupyterlab on localhost where you can run your fenicsx codes.

2. To install additional packages for visualisation, you need to make a Dockerfile and build your own image derived from the above fenicsx image. Add below code in your Dockerfile.
```
# Use the existing dolfinx image as the base image
FROM dolfinx/lab:stable
# Set environment variables to prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
# Set environment variable for display
ENV DISPLAY=host.docker.internal:0
# Install additional system packages
RUN apt-get update && apt-get install -y \
    libglu1-mesa-dev freeglut3-dev mesa-common-dev \
    && rm -rf /var/lib/apt/lists/*
# Install Python packages using pip
RUN pip3 install \
    sympy \
    ipywidgets \
    ipygany \
    trame \
    trame-vuetify \
    trame-vtk  \
    trame_jupyter_extension
# Set the working directory
WORKDIR /
# Expose port for Jupyter Notebook
EXPOSE 8888
```

3. In order to launch Gmsh GUI inside your docker container, you need X11 Forwarding on your Windows. 
    - Download Vcxsrv from <https://sourceforge.net/projects/vcxsrv/> . 
    - Launch VcXsrv using the "Multiple windows" option.
    - Ensure that the display number matches (typically 0).
    - Enable the following options:"Disable access control" (to allow connections from Docker).


4. To build your image open Command Prompt and go to the folder where the above Dockerfile is saved. Run the below command with your image name.
```shell
docker build -t <name of your image> .
```


5. Now to maintain your fenicsx files in your host system you would need volume mounting. Make a folder where you want to maintain your file system. This folder will share files with the directory '/home/fenicsx/shared' which is inside your container. Make a batch file (.bat) and add the below command in it.
```shell
docker run -it --rm -v <path to your folder>:/home/fenicsx/shared/ -p 8888:8888 -e DISPLAY=host.docker.internal:0 <name of your image>
```

6. Now, whenever you want to work on Fenicsx, make sure the docker engine is turned on and double click on the batch file you made. It will launch the Jupyterlab with all the requirements. 

7. To test you installation launch jupyter notebook.

```shell
jupyter notebook
```

Make a new notebook and add below code performing linear elastic simulation of cantilever under UDL.

```python
#Import the necessary packages
import pyvista
from dolfinx import mesh, fem, plot, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import numpy as np
import gmsh
import sys
```

```python
#Function to make the geometry and mesh in Gmsh
def gmsh_beam(L,W,T):

gmsh.initialize()
gmsh.option.setNumber('General.Verbosity', 0)
gmsh.model.add("beam")
gmsh.model.occ.addBox(0,0,0,L,W,T)

gmsh.model.occ.synchronize()

edges = gmsh.model.getEntities(1)
faces = gmsh.model.getEntities(2)
volumes = gmsh.model.getEntities(3)

for i in edges:
length = gmsh.model.occ.getMass(1,i[1])
if length == L:
    gmsh.model.mesh.setTransfiniteCurve(i[1],100)
    
elif length == W:
    gmsh.model.mesh.setTransfiniteCurve(i[1],10)
    
else:
    gmsh.model.mesh.setTransfiniteCurve(i[1],10)

for i in faces:
gmsh.model.mesh.setTransfiniteSurface(i[1])
gmsh.model.mesh.setRecombine(i[0],i[1])

for i in volumes:
gmsh.model.mesh.setTransfiniteVolume(i[1])


gmsh.model.addPhysicalGroup(2,[1],101)
gmsh.model.setPhysicalName(2,101,"fixed_surface")
gmsh.model.addPhysicalGroup(2,[6],102)
gmsh.model.setPhysicalName(2,102,"pressure_surface")

gmsh.model.addPhysicalGroup(3,[1],103)
gmsh.model.setPhysicalName(3,103,"volume1")

gmsh.model.mesh.generate(3)

if '-nopopup' not in sys.argv:
gmsh.fltk.run()

domain, markers, facets = io.gmshio.model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, gdim=3)

gmsh.finalize()

return domain,markers,facets
```

```python
#Function to setup the FE case and solve
def fe_beam(L,W,T):

domain, markers, facets = gmsh_beam(L,W,T)

rho = 0;g = 1
E = 1e5;nu = 0.3
mu = E/2/(1+nu);lambda_ = E*nu/(1+nu)/(1-2*nu)

# Define strain and stress
def epsilon(u):
    return ufl.sym(ufl.grad(u))  # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)

def sigma(u):
    return lambda_ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2 * mu * epsilon(u)


V = fem.functionspace(domain, ("Lagrange", 1, (3,)))

bc_dofs = fem.locate_dofs_topological(V, 2, facets.find(101))
bcs = fem.dirichletbc(np.zeros((3,)), bc_dofs, V)


ds = ufl.Measure("ds", domain=domain,subdomain_data=facets)

n = ufl.FacetNormal(domain)
u = ufl.TrialFunction(V)
d = domain.geometry.dim
v = ufl.TestFunction(V)
f = fem.Constant(domain, default_scalar_type((0, 0, 0)))
T = fem.Constant(domain,-1.0e1)
a = ufl.inner(sigma(u), epsilon(v))*ufl.dx
L = ufl.dot(f, v)*ufl.dx + ufl.dot(T*n, v)*ds(102)

# Compute solution
problem = LinearProblem(a, L, bcs=[bcs], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

pyvista.start_xvfb()

# Create plotter and pyvista grid
p = pyvista.Plotter()
topology, cell_types, geometry = plot.vtk_mesh(V)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

# Attach vector values to grid and warp grid by vector
grid["u"] = uh.x.array.reshape((geometry.shape[0], 3))
actor_0 = p.add_mesh(grid, style="wireframe", color="k")
warped = grid.warp_by_vector("u", factor=0.1)
actor_1 = p.add_mesh(warped, show_edges=True)
p.show_axes()
if not pyvista.OFF_SCREEN:
    p.show()
    
print("FEM Max Deflection = " + str(abs(uh.vector.min()[1])))
print("Analytical Max Deflection = w*L^4/(8*E*I) = " + str(1.5))
```

```python
fe_beam(10,1,1)
```

You should get the below output. 

<iframe src="{{ '/assets/fenicsx_plots/beam.html' | relative_url }}" width="100%" height="600px" style="border:none;"></iframe>


## References

[1] <https://github.com/FEniCS/dolfinx>  
[2] <https://bleyerj.github.io/comet-fenicsx/>  
[3] <https://jsdokken.com/dolfinx-tutorial/>  