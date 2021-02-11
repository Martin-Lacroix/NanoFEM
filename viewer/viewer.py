import numpy as np
import gmsh

# %% Load Data

VM = np.loadtxt(r"..\output\stress.txt")
elem = np.loadtxt(r"..\output\elements.txt",delimiter=",")
u = np.loadtxt(r"..\output\displacement.txt",delimiter=",")
nXYZ = np.loadtxt(r"..\output\coordinates.txt",delimiter=",")

# Cubic domain size

xLen = [np.min(nXYZ[:,0]),np.max(nXYZ[:,0])]
yLen = [np.min(nXYZ[:,1]),np.max(nXYZ[:,1])]
zLen = [np.min(nXYZ[:,2]),np.max(nXYZ[:,2])]
eLen = VM.shape[0]
nLen = u.shape[0]

# Stores the data for Gmsh

VM = [[VM[i]] for i in range(eLen)]
ux = [[u[i,0]] for i in range(nLen)]
uy = [[u[i,1]] for i in range(nLen)]
uz = [[u[i,2]] for i in range(nLen)]
U = [[np.linalg.norm(u[i])] for i in range(nLen)]

# %% Domain Geometry

gmsh.initialize()
gmsh.model.add("Nascam")

# Model points

for i in range(2):
    for j in range(2):
        for k in range(2):
            gmsh.model.geo.addPoint(xLen[i],yLen[j],zLen[k])

# Model lines

gmsh.model.geo.addLine(1,2)
gmsh.model.geo.addLine(2,4)
gmsh.model.geo.addLine(4,3)

gmsh.model.geo.addLine(3,1)
gmsh.model.geo.addLine(1,5)
gmsh.model.geo.addLine(5,7)

gmsh.model.geo.addLine(7,3)
gmsh.model.geo.addLine(4,8)
gmsh.model.geo.addLine(8,6)

gmsh.model.geo.addLine(6,2)
gmsh.model.geo.addLine(7,8)
gmsh.model.geo.addLine(6,5)

# Model faces

gmsh.model.geo.addCurveLoop([1,2,3,4])
gmsh.model.geo.addCurveLoop([4,5,6,7])
gmsh.model.geo.addCurveLoop([2,8,9,10])
gmsh.model.geo.addCurveLoop([11,9,12,6])
gmsh.model.geo.addCurveLoop([3,-7,11,-8])
gmsh.model.geo.addCurveLoop([1,-10,12,-5])

for i in range(1,7):
    gmsh.model.geo.addPlaneSurface([i])

# Model volume

gmsh.model.geo.addSurfaceLoop([1,2,3,4,5,6])
gmsh.model.geo.addVolume([1])

# %% Makes the Mesh

coord = nXYZ.flatten()
nTag = np.arange(nXYZ.shape[0])+1
eTag = np.arange(elem.shape[0])+1

for i in range(elem.shape[0]):
    
    n = elem[i]
    elem[i] = [n[0],n[4],n[6],n[2],n[1],n[5],n[7],n[3]]
    
eNode = elem.flatten()+1

# Adds the nodes to the mesh

gmsh.model.geo.synchronize()
gmsh.model.mesh.addNodes(3,1,nTag,coord)
gmsh.model.mesh.addElements(3,1,[5],[eTag],[eNode])

# %% Plots the solution

gmsh.view.add("Displacement x",1)
gmsh.view.add("Displacement y",2)
gmsh.view.add("Displacement z",3)
gmsh.view.add("Total Displacement",4)
gmsh.view.add("Von Mises Stress",5)

# Adds the data to the model

gmsh.view.addModelData(1,0,"Nascam","NodeData",nTag,ux)
gmsh.view.addModelData(2,0,"Nascam","NodeData",nTag,uy)
gmsh.view.addModelData(3,0,"Nascam","NodeData",nTag,uz)
gmsh.view.addModelData(4,0,"Nascam","NodeData",nTag,U)
gmsh.view.addModelData(5,0,"Nascam","ElementData",eTag,VM)

# Writes the solution

gmsh.write(r"..\output\result.msh")
for i in range(1,6): gmsh.view.write(i,r"..\output\result.msh",append=True)