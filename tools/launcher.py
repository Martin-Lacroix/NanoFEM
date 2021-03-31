import numpy as np
import subprocess
import shutil
import time
import os

# %% Functions

def replace(file,num,text):
    
    lines = open(file,'r').readlines()
    lines[num] = text
    out = open(file, 'w')
    out.writelines(lines)
    out.close()

# %% Main Code

step = 10
scale = 0.34
height = np.arange(10,1310,step)
path = r'C:/Users/ORBBE/Desktop/CGFEM-Nano/'
crop = (1439-height)*scale


for i in range(height.shape[0]):
    
    # Updates the parameters

    start = time.time()
    replace(path+'input.txt',1,'1;0.34;'+str(crop[i])+'!General Parameters\n')
    
    # Runs the algorithm
    
    FNULL = open(os.devnull,'w')
    subprocess.call(path+'Mechanics.exe',stdout=FNULL,stderr=subprocess.STDOUT,cwd=path)
    subprocess.call(['python','cleaner.py'],stdout=FNULL,stderr=subprocess.STDOUT,cwd=path+'viewer/')
    
    # Moves the results
    
    out = r'E:\Results\Big Double Cu Zr\Bulk\Height '+str(height[i])
    shutil.move(path+'output',out)
    
    # Timer and iterations
    
    current = height[i]
    laps = (time.time()-start)/60
    print('Height =',height[i],'\tTime =',round(laps,1),'min')