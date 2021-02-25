import matplotlib.pyplot as plt
import numpy as np

# %% Functions

def fixed(K,B,dim,node,val):
    
    nLen = B.shape[0]//3
    
    for i in node:
        
        B -= K[:,i+dim*nLen]*val
        B[i+dim*nLen] = val
        
        K[i+dim*nLen,:] = 0
        K[:,i+dim*nLen] = 0
        K[i+dim*nLen,i+dim*nLen] = 1
    
    return K,B

# Pair = [xi,xj] => xi' = xi-xj

def delta(K,B,dim,pair):
    
    nLen = B.shape[0]//3
    
    i = pair[0]
    j = pair[1]
    
    K[j+dim*nLen,:] += K[i+dim*nLen,:]
    K[:,j+dim*nLen] += K[:,i+dim*nLen]
    B[j+dim*nLen] += B[i+dim*nLen]
    
    return K,B

def coupled(K,B,dim,coup):
    
    nLen = B.shape[0]//3
    j = coup[-1]
    i = coup[0]
    
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
M = np.loadtxt("output/M.txt")
B = np.loadtxt("output/B.txt")

# K_old = np.loadtxt("output/K_old.txt")
# M_old = np.loadtxt("output/M_old.txt")
# B_old = np.loadtxt("output/B_old.txt")

# Re-build the full matrix

# K_old = makeFull(K_old)
# K = makeFull(K)

# sameK = np.allclose(K_old,K)
# sameB = np.allclose(B_old,B)

# %% Solves

dt = 0.001
nLen = B.shape[0]
u1 = np.zeros(nLen)
u2 = np.zeros(nLen)
u3 = np.zeros(nLen)

M = M+M.T
K = K+K.T

for i in range(nLen): M[i,i] /= 2
for i in range(nLen): K[i,i] /= 2

dx = []

for i in range(100):

    M_new = M.copy()
    RHS = dt**2*B+M.dot(2*u2-u1)-dt**2*K.dot(u2)
    
    u12 = 2*u2-u1
    Mu12 = M.dot(2*u2-u1)
    
    fixed(M_new,RHS,0,[0,1,2,3],0)
    fixed(M_new,RHS,1,[0,1,4,5],0)
    fixed(M_new,RHS,2,[0,2,4,6],0)
    
    coupled(M_new,RHS,0,[4,5,6,7])
    
    u3 = np.linalg.solve(M_new,RHS)
    u3[[4,5,6]] = u3[7]
    u1 = u2.copy()
    u2 = u3.copy()
    
    dx.append(u3[4])





plt.plot(dx)



