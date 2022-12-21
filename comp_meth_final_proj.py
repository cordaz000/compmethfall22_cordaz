#final project, attempt at solution of Barotropic Potential Vorticity equation

#handle imports 
import numpy as np
from random import seed,random
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import time as Time
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

from functions_final_proj import spatial_grid_and_time_steps
from functions_final_proj import GFD_constants
from functions_final_proj import initial_conditions
from functions_final_proj import main_time_loop

Rearth = 6.371*10**6 
SpatialArraysTimeStepsnSpatialScales = spatial_grid_and_time_steps(1, 100, 100, 1/12, 0, 5500,Rearth,1)  

xx = SpatialArraysTimeStepsnSpatialScales[0]
yy = SpatialArraysTimeStepsnSpatialScales[1]
TimeSteps = SpatialArraysTimeStepsnSpatialScales[2]
spatialExtentX = SpatialArraysTimeStepsnSpatialScales[3]
spatialExtentY = SpatialArraysTimeStepsnSpatialScales[4]

ComputationConstants = GFD_constants(R_Earth,45,5500)
CoriolisParam = ComputationConstants[0]
BetaParam = ComputationConstants[1]
DeformationWaveVector = ComputationConstants[2]

Psi_ini = initial_conditions(100,100, xx,yy,350,spatialExtentX,spatialExtentY,5500,CoriolisParam,BetaParam,30)

main_time_loop(Psi_ini, numberoftimes, 100, 100,xx,yy,times2plot,DeformationWaveVector**2)