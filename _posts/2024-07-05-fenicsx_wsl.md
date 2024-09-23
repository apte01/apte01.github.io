---
layout: post
title: 'Setting up Fenicsx on Windows Subsystem for Linux (WSL)'
date: 2024-07-05
permalink: /posts/2024/07/fenicsx_wsl/
tags:
excerpt: 
---
<br><br>

## Installing and setting up WSL

1. If you are installing WSL for first time, the below command will install Ubuntu distribution.
```shell
wsl --install
```
2. If you already have WSL installed, it is better to create separate instance specifically for fenicsx. To clone the existing instance,export it to a .tar file and then import the same .tar file for the new instance.
```shell
wsl --export <name of your existing distribution> <file.tar>
```

```shell
wsl --import <name of your new instance> <InstallLocation> <file.tar>
```
3. To add a new user
```bash
adduser <username>
```
    To give sudo permission to the new user.
```bash
usermod -aG sudo <username>
```

4. To make the new user as the default user add the below lines to ./etc/wsl.conf file.
```
[user]
default=<username>
```

5. For more details, visit <https://learn.microsoft.com/en-us/windows/wsl/install>


## Installing Fenicsx and other requirements

1. Install Fenicsx
```bash
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt update
sudo apt install fenicsx
```

2. Install Gmsh for meshing
```bash
pip install gmsh
```

3. Install Jupyter
```bash
sudo apt install jupyter
```

4. To visualize in jupyter notebook
```bash
pip install --upgrade pyvista ipywidgets ipygany
sudo apt-get -y install libglu1 libxcursor-dev libxft2 libxinerama1 libfltk1.3-dev libfreetype6-dev libgl1-mesa-dev libgl1-mesa-glx xvfb
pip install trame
pip install trame-vuetify
pip install --upgrade trame-vtk
pip install "pandas==1.5.3"
pip install openpyxl
```


5. To test you installation launch jupyter notebook.
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