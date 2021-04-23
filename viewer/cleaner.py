import numpy as np
import gmsh

# %% Load Data

s = np.loadtxt(r"..\output\stress.txt")
nXYZ = np.loadtxt(r"..\output\node.txt",delimiter=",")
u = np.atleast_2d(np.loadtxt(r"..\output\disp.txt",delimiter=","))
elem = np.atleast_2d(np.loadtxt(r"..\output\elem.txt",delimiter=",")+0.1).astype(int)

# Computes the cubic domain size

xLen = [np.min(nXYZ[:,0]),np.max(nXYZ[:,0])]
yLen = [np.min(nXYZ[:,1]),np.max(nXYZ[:,1])]
zLen = [np.min(nXYZ[:,2]),np.max(nXYZ[:,2])]

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

E = np.zeros(elem.shape[0]).astype(bool)
N = np.zeros(nXYZ.shape[0]).astype(bool)

for i in range(elem.shape[0]):
    
    if(elem.shape[1]==8): elem[i] = elem[i][[0,4,6,2,1,5,7,3]]
    else: elem[i] = elem[i][[0,18,20,2,6,24,26,8,9,1,3,19,21,11,23,5,15,7,25,17,10,12,4,22,14,16,13]]

    if(abs(s[i])!=0):
        
        N[elem[i]] = True
        E[i] = True
        
# Deletes empty elements

nTag = np.arange(nXYZ.shape[0])+1
eTag = np.arange(elem.shape[0])+1

eTag = eTag[E]
nTag = nTag[N]

elem = elem[E]
nXYZ = nXYZ[N]

s = s[E]
u = u[N]

coord = nXYZ.flatten()
eNode = elem.flatten()+1

# Adds the nodes to the mesh

gmsh.model.geo.synchronize()
gmsh.model.mesh.addNodes(3,1,nTag,coord)

if(elem.shape[1]==8): gmsh.model.mesh.addElements(3,1,[5],[eTag],[eNode])
else: gmsh.model.mesh.addElements(3,1,[12],[eTag],[eNode])

# %% Extracts the scalars

eLen = elem.shape[0]
nLen = nXYZ.shape[0]

# Stores the displacement field
    
ux = [[u[i,0]] for i in range(nLen)]
uy = [[u[i,1]] for i in range(nLen)]
uz = [[u[i,2]] for i in range(nLen)]

# Stores the stress field

vm = [[s[i]] for i in range(eLen)]

# %% Plots the solution

gmsh.view.add("ux",1)
gmsh.view.add("uy",2)
gmsh.view.add("uz",3)

# Stress field

gmsh.view.add("vm",4)

# Writes the displacement data in the model
    
gmsh.view.addModelData(1,0,"Nascam","NodeData",nTag,ux)
gmsh.view.addModelData(2,0,"Nascam","NodeData",nTag,uy)
gmsh.view.addModelData(3,0,"Nascam","NodeData",nTag,uz)

gmsh.view.write(1,r"..\output\disp-X.msh")
gmsh.view.write(2,r"..\output\disp-Y.msh")
gmsh.view.write(3,r"..\output\disp-Z.msh")

# Writes the stress data in the model

gmsh.view.addModelData(4,0,"Nascam","ElementData",eTag,vm)
gmsh.view.write(4,r"..\output\stress.msh")