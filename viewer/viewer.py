import numpy as np
import gmsh

# %% Load Data

nXYZ = np.loadtxt(r"..\output\coordinates.txt",delimiter=",")
elem = np.atleast_2d(np.loadtxt(r"..\output\elements.txt",delimiter=","))
uList = np.atleast_2d(np.loadtxt(r"..\output\displacement.txt",delimiter=","))

eLen = elem.shape[0]
nLen = nXYZ.shape[0]

# Computes the cubic domain size

xLen = [np.min(nXYZ[:,0]),np.max(nXYZ[:,0])]
yLen = [np.min(nXYZ[:,1]),np.max(nXYZ[:,1])]
zLen = [np.min(nXYZ[:,2]),np.max(nXYZ[:,2])]
step = len(uList)

# Stores the displacement field

ux = [0]*step
uy = [0]*step
uz = [0]*step

for k in range(step):
    
    ux[k] = [[uList[k][i]] for i in range(nLen)]
    uy[k] = [[uList[k][i+nLen]] for i in range(nLen)]
    uz[k] = [[uList[k][i+2*nLen]] for i in range(nLen)]

# Stores the stress field

try:
    sig = []
    s = np.atleast_2d(np.loadtxt(r"..\output\stress.txt",delimiter=","))
    for j in range(6): sig.append([[s[i,j]] for i in range(eLen)])

except: pass


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
    
    if(elem.shape[1]==8): elem[i] = elem[i][[0,4,6,2,1,5,7,3]]
    else: elem[i] = elem[i][[0,18,20,2,6,24,26,8,9,1,3,19,21,11,23,5,15,7,25,17,10,12,4,22,14,16,13]]
    
    
coord = nXYZ.flatten()
eNode = elem.flatten()+1
nTag = np.arange(nXYZ.shape[0])+1
eTag = np.arange(elem.shape[0])+1

# Adds the nodes to the mesh

gmsh.model.geo.synchronize()
gmsh.model.mesh.addNodes(3,1,nTag,coord)

if(elem.shape[1]==8): gmsh.model.mesh.addElements(3,1,[5],[eTag],[eNode])
else: gmsh.model.mesh.addElements(3,1,[12],[eTag],[eNode])

# %% Plots the solution

gmsh.view.add("Displacement x",1)
gmsh.view.add("Displacement y",2)
gmsh.view.add("Displacement z",3)

# Stress field

gmsh.view.add("Stress xx",4)
gmsh.view.add("Stress yy",5)
gmsh.view.add("Stress zz",6)
gmsh.view.add("Stress yz",7)
gmsh.view.add("Stress zx",8)
gmsh.view.add("Stress xy",9)

# Writes the displacement data in the model

for i in range(step):
    
    gmsh.view.addModelData(1,i,"Nascam","NodeData",nTag,ux[i])
    gmsh.view.addModelData(2,i,"Nascam","NodeData",nTag,uy[i])
    gmsh.view.addModelData(3,i,"Nascam","NodeData",nTag,uz[i])

gmsh.write(r"..\output\displacement.msh")
for i in range(3): gmsh.view.write(i+1,r"..\output\displacement.msh",append=True)

# Writes the stress data in the model

if(len(sig)>0):
    
    gmsh.write(r"..\output\stress.msh")
    for i in range(6): gmsh.view.addModelData(i+4,0,"Nascam","ElementData",eTag,sig[i])
    for i in range(6): gmsh.view.write(i+4,r"..\output\stress.msh",append=True)