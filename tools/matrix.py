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

U = np.loadtxt("../output/disp.txt",delimiter=",")
nXYZ = np.loadtxt("../output/node.txt",delimiter=",")

K = np.loadtxt("../output/K.txt")
B = np.loadtxt("../output/B.txt")
# M = np.loadtxt("../output/M.txt")
