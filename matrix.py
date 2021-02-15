import numpy as np
import os

# %% Functions

def fixed(K,B,dim,node,val):
    
    for i in node:
        
        B -= K[:,i+dim*nLen]*val
        B[i+dim*nLen] = val
        
        K[i+dim*nLen,:] = 0
        K[:,i+dim*nLen] = 0
        K[i+dim*nLen,i+dim*nLen] = 1
    
    return K,B

# Pair = [xi,xj] => xi' = xi-xj

def delta(K,B,dim,pair):
    
    i = pair[0]
    j = pair[1]
    
    K[j+dim*nLen,:] += K[i+dim*nLen,:]
    K[:,j+dim*nLen] += K[:,i+dim*nLen]
    B[j+dim*nLen] += B[i+dim*nLen]
    
    return K,B

def coupled(K,B,dim,coup):
    
    i = coup[0]
    j = coup[-1]
    
    for i in coup:
        if(i!=j):
            
            K[j+dim*nLen,:] += K[i+dim*nLen,:]
            K[:,j+dim*nLen] += K[:,i+dim*nLen]
            B[j+dim*nLen] += B[i+dim*nLen]
            
            K[i+dim*nLen,:] = 0
            K[:,i+dim*nLen] = 0
            K[i+dim*nLen,i+dim*nLen] = 1
            B[i+dim*nLen] = 0
    
    return K,B

def makeFull(K):
    
    K2 = K.T.copy()
    for i in range(len(K)): K2[i,i] = 0
    K += K2
    return K

    
# %% Load files

U = np.loadtxt("output/displacement.txt",delimiter=",")
nXYZ = np.loadtxt("output/coordinates.txt",delimiter=",")
nLen = U.shape[0]

K = np.loadtxt("output/K.txt")
B = np.loadtxt("output/B.txt")

K_old = np.loadtxt("output/K_old.txt")
B_old = np.loadtxt("output/B_old.txt")

# Re-build the full matrix

K_old = makeFull(K_old)
K = makeFull(K)

sameK = np.allclose(K_old,K)
sameB = np.allclose(B_old,B)

# %% Old Solution

u = np.zeros(3*nLen)
for i in range(nLen):
    u[i+2*nLen] = U[i,2]
    u[i+nLen] = U[i,1]
    u[i] = U[i,0]

# %% Test Periodic

# pair0 = np.zeros((9,2)).astype(int)
# pair0[:,0] = [18,19,20,21,22,23,24,25,26]
# pair0[:,1] = [0,1,2,3,4,5,6,7,8]
# for i in range(pair0.shape[0]): Kf,B = delta(Kf,B,0,pair0[i])

# pair1 = np.zeros((9,2)).astype(int)
# pair1[:,0] = [2,5,8,11,14,17,20,23,26]
# pair1[:,1] = [0,3,6,9,12,15,18,21,24]
# for i in range(pair1.shape[0]): Kf,B = delta(Kf,B,1,pair1[i])

# coup0 = [18,19,20,21,22,23,24,25,26]
# Kf,B = coupled(Kf,B,0,coup0)


coup2 = [1,3,5,7,9,11,13,15,17]
K,B = coupled(K,B,2,coup2)

K,B = fixed(K,B,0,[0,1,2,3,4,5,6,7,8],0)
K,B = fixed(K,B,2,[0,3,6,9,12,15,18,21,24],0)


u = np.linalg.solve(K,B)

# for i in range(len(coup0)): u[coup0[i]] = u[coup0[-1]]
# for i in range(len(coup1)): u[coup1[i]+nLen] = u[coup1[-1]+nLen]

# for i in range(pair0.shape[0]): u[pair0[i,0]] += u[pair0[i,1]]
# for i in range(pair1.shape[0]): u[pair1[i,0]+nLen] += u[pair1[i,1]+nLen]

# %% Prints

print("\n")
for i in range(nLen):print("ux[",i,"] = ",round(u[i],4))
print("\n")
for i in range(nLen):print("uy[",i,"] = ",round(u[i+nLen],4))
print("\n")
for i in range(nLen):print("uz[",i,"] = ",round(u[i+2*nLen],4))

file = open("output/displacement.txt", "w")
for i in range(nLen):
    file.write(str(u[i]))
    file.write(",")
    file.write(str(u[i+nLen]))
    file.write(",")
    file.write(str(u[i+2*nLen]))
    file.write("\n")
               
file.close()

file = open("output/coordinates.txt", "w")
for i in range(nLen):
    file.write(str(nXYZ[i,0]+u[i]))
    file.write(",")
    file.write(str(nXYZ[i,1]+u[i+nLen]))
    file.write(",")
    file.write(str(nXYZ[i,2]+u[i+2*nLen]))
    file.write("\n")
               
file.close()

