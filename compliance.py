import numpy as np

# %% Input Data

Lx = 200
Ly = 100
Lz = 270

sxx = 100
syy = 100
szz = 100

sxy = 100
szy = 100
szx = 100

# %% Simulation Data

dx_X = 0.276
dy_X = -0.0395
dz_X = -0.0905

dx_Y = -0.072
dy_Y = 0.14
dz_Y = -0.0905

dx_Z = -0.0714
dy_Z = -0.0375
dz_Z = 0.36

dx_ZX = 0.93
dy_XY = 0.751
dy_ZY = 0.936

# %% Strain components

exx_X = dx_X/Lx
exx_Y = dx_Y/Lx
exx_Z = dx_Z/Lx

eyy_X = dy_X/Ly
eyy_Y = dy_Y/Ly
eyy_Z = dy_Z/Ly

ezz_X = dz_X/Lz
ezz_Y = dz_Y/Lz
ezz_Z = dz_Z/Lz

gzy = dy_ZY/Lz
gzx = dx_ZX/Lz
gxy = dy_XY/Lx

# %% Compliance Tensor

S11 = exx_X/sxx
S22 = eyy_Y/syy
S33 = ezz_Z/szz

S12 = exx_Y/syy
S21 = eyy_X/sxx
S13 = exx_Z/szz
S31 = ezz_X/sxx
S23 = eyy_Z/szz
S32 = ezz_Y/syy

S44 = gzy/szy
S55 = gzx/szx
S66 = gxy/sxy

# Builds the matrix

S = np.zeros((6,6))
S[0,0] = S11
S[0,1] = S12
S[0,2] = S13
S[1,0] = S21
S[1,1] = S22
S[1,2] = S23
S[2,0] = S31
S[2,1] = S32
S[2,2] = S33
S[3,3] = S44
S[4,4] = S55
S[5,5] = S66

# Stiffness Tensor

C = np.linalg.inv(S)
A = np.tensordot(C,S)-6