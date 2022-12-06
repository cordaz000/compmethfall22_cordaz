
import numpy as np
import matplotlib.pyplot as plt


from functions_hw8 import Second2FirstOrderODE
from functions_hw8 import Runge_Kutta_4th_Order
from functions_hw8 import TrajectoryPlotter


#define constant mass
mass1 = 1 #units of kg

a = 0.0
b = 200.0
N = 1000.0
h = (b-a)/N

tpoints = np.arange(a,b,h)
xpoints = list()
ypoints = list()

v_x_ini = 100*np.cos(np.deg2rad(30))
v_y_ini = 100*np.sin(np.deg2rad(30)) 

r = np.array([0.0,v_x_ini,0.0,v_y_ini],float) #vector holding initial conditions
epsilon = 10e-5 #use this small number to compare floats

exes,whys = Runge_Kutta_4th_Order(tpoints,r,h,xpoints,ypoints,Second2FirstOrderODE,epsilon,mass1)

TrajectoryPlotter(exes,whys)
plt.show()

#see effect of different masses on trajectory
masses = np.linspace(0.2, 0.8,4) #this has 4 masses

print('\n')
legends = ['ZeroPointTwo_kg', 'ZeroPointFour_kg', 'ZeroPointSix_kg', 'ZeroPointEight_kg']

dic_holder = {mass:{'x': list(), 'y':list()} for mass in legends}

m = 0
for key, value in dic_holder.items():
    Runge_Kutta_4th_Order(tpoints,r,h,value['x'],value['y'],\
        Second2FirstOrderODE,epsilon,masses[m])
    m = m + 1

fig, ax = plt.subplots()
ax.scatter(dic_holder['ZeroPointTwo_kg']['x'],dic_holder['ZeroPointTwo_kg']['y'],label = '0.2 kg')
ax.scatter(dic_holder['ZeroPointFour_kg']['x'],dic_holder['ZeroPointFour_kg']['y'],label = '0.4 kg')
ax.scatter(dic_holder['ZeroPointSix_kg']['x'],dic_holder['ZeroPointSix_kg']['y'],label = '0.6 kg')
ax.scatter(dic_holder['ZeroPointEight_kg']['x'],dic_holder['ZeroPointEight_kg']['y'],label = '0.8 kg')
ax.set_xlabel('Distance [m]')
ax.set_ylabel('Height [m]')
ax.set_title("Different masses cannonball's trajectories")
ax.legend()
plt.show()
