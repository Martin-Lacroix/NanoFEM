import numpy as np
import gmsh

# %% Load Data

s = np.loadtxt(r"..\output\stress.txt",delimiter=",")
elem = np.loadtxt(r"..\output\elements.txt",delimiter=",")
u = np.loadtxt(r"..\output\displacement.txt",delimiter=",")
nXYZ = np.loadtxt(r"..\output\coordinates.txt",delimiter=",")

s = np.atleast_2d(s)
elem = np.atleast_2d(elem)

# Cubic domain size

xLen = [np.min(nXYZ[:,0]),np.max(nXYZ[:,0])]
yLen = [np.min(nXYZ[:,1]),np.max(nXYZ[:,1])]
zLen = [np.min(nXYZ[:,2]),np.max(nXYZ[:,2])]
eLen = s.shape[0]
nLen = u.shape[0]
U = []
S = []

# Stores the displacement and stress field

for j in range(3): U.append([[u[i,j]] for i in range(nLen)])
for j in range(6): S.append([[s[i,j]] for i in range(eLen)])

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

for i in range(elem.shape[0]):
    
    n = elem[i]
    elem[i] = [n[0],n[4],n[6],n[2],n[1],n[5],n[7],n[3]]
    
coord = nXYZ.flatten()
eNode = elem.flatten()+1
nTag = np.arange(nXYZ.shape[0])+1
eTag = np.arange(elem.shape[0])+1

# Adds the nodes to the mesh

gmsh.model.geo.synchronize()
gmsh.model.mesh.addNodes(3,1,nTag,coord)
gmsh.model.mesh.addElements(3,1,[5],[eTag],[eNode])

# %% Plots the solution

gmsh.view.add("Displacement x",1)
gmsh.view.add("Displacement y",2)
gmsh.view.add("Displacement z",3)

# Stress field

gmsh.view.add("Stress xx",4)
gmsh.view.add("Stress yy",5)
gmsh.view.add("Stress zz",6)
gmsh.view.add("Stress xy",7)
gmsh.view.add("Stress zx",8)
gmsh.view.add("Stress yz",9)

# Adds the data to the model

for i in range(3): gmsh.view.addModelData(i+1,0,"Nascam","NodeData",nTag,U[i])
for i in range(6): gmsh.view.addModelData(i+4,0,"Nascam","ElementData",eTag,S[i])

# Creates the output files

gmsh.write(r"..\output\stress.msh")
gmsh.write(r"..\output\displacement.msh")

# Writes the solutions

for i in range(3): gmsh.view.write(i+1,r"..\output\displacement.msh",append=True)
for i in range(6): gmsh.view.write(i+4,r"..\output\stress.msh",append=True)