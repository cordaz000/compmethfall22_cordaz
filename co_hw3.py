#
import numpy as np


from hw3_functions import potential_calculator
from hw3_functions import potential_plotter
from hw3_functions import partial_derivative
from hw3_functions import electric_field_plotter


#define grid
gx = np.linspace(-0.50, 0.50, 101)
gy = np.linspace(-0.50, 0.50, 101) 

#grid is defined in meters
gxx, gyy = np.meshgrid(gx, gy) 

#first charge is one Coulumb
q_1 = 1 
qloc1 = 0, 0.005 # place first charge at x = 0, y = 0.5 cm

#second charge is negative one Coulumb
q_2 = -1
qloc2 = 0, -0.005 #place second charge at x = 0, y = -0.5 cm

pot1 = potential_calculator(gxx,gyy, q_1,qloc1)
pot2 = potential_calculator(gxx, gyy, q_2,qloc2)

#superposition, add the potential fields from each charge
Tot_pot = pot1 + pot2

levels = [np.percentile(Tot_pot,percent) for percent in np.linspace(5,95, 19)]
potential_plotter(gxx,gyy, Tot_pot, levels)

#compute the gradient of the Potential field (Tot_pot)
#do it wrt to x, or columns
partial_x = -partial_derivative(Tot_pot,x=True, y=False)
partial_y = -partial_derivative(Tot_pot, x = False, y = True)

#reduce the size of the arrays for plotting the vector field
p_x_r = partial_x[::2, ::2] # this is a (51, 51)
p_y_r = partial_y[::2,::2] #same size as as above

#make a grid that is (51, 51)
grid_array = np.linspace(-0.50, 0.50, 51) #start stop # of steps
xxx, yyy = np.meshgrid(grid_array,grid_array)

electric_field_plotter(xxx,yyy,p_x_r,p_y_r)