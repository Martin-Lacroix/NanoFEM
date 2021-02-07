import numpy as np

# %% Load files

u = np.loadtxt("output/displacement.txt")
nXYZ = np.loadtxt("output/coordinates.txt",delimiter=",")

K = np.loadtxt("output/K.txt")
B = np.loadtxt("output/B.txt")


K_old = np.loadtxt("output/K_old.txt")
B_old = np.loadtxt("output/B_old.txt")

sameK = np.allclose(K,K_old)
sameB = np.allclose(B,B_old)

K2 = K.T.copy()
for i in range(len(K)): K2[i,i] = 0
K2 = K+K2

zero = K_old-K2
print(np.max(abs(zero)))