import numpy as np
import matplotlib.pyplot as plt

from scipy.constants import epsilon_0 as ep

def potential_calculator(xx,yy,q_val,q_loc):

    """
    Provided a grid (two arrays), value of charge and its 
    location, return the potential field on the
    specified grid.
    """
    const = (4*np.pi*ep)**-1 #units of kg m**2 / s**4 A**2

    #q_loc is a tuple, with first entry x coord and 2nd y val
    distance = np.sqrt( (xx - q_loc[0])**2 + (yy - q_loc[1])**2 )
    inverted_distance = distance**-1

    Potential = const*inverted_distance*q_val #units of kg m / s**3 A = J/C = V
    return Potential

def potential_plotter(xx,yy,Pot_vals,levels):
    """
    Plot potential values.
    """
    fig, ax = plt.subplots()
    pot1 = ax.contourf(xx,yy, Pot_vals,levels = levels, cmap = 'bwr', extend = 'both')
    _= ax.contour(xx,yy,Pot_vals,levels = levels,cmap = 'bwr', linewidths = 3)
    ax.set_title('Potential for two point charges [V]')
    
    #remove axis labels
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.tick_params(bottom = False)
    ax.tick_params(left = False)
    fig.colorbar(pot1)
    plt.show()

def partial_derivative(A,x,y):
    """
    Input an array to compute its partial derivative.
    That array is A--what the input function takes. Specify whether you want
    to differentiate wrt to x or wrt y in the second and third arg respectively.
    """
    #make an array to hold values obtained via central diff method
    fill_me = np.empty((A.shape))

    if x:
        #handle first edge, where you can do central diff, and outer edge
        for i in range(A.shape[1]): #columns indexed
            for j in range(A.shape[0]): #rows indexed
                if i == 0:
                    fill_me[j,i] = (A[j,i + 1] - A[j,i])/2 #first column
                    
                elif (i != 0) and (i != A.shape[1] - 1): #central diff applicable
                    fill_me[j,i] = (A[j,i + 1] - A[j,i - 1])/2
                    
                elif i == A.shape[1] - 1: #the edge
                    fill_me[j,i] = (A[j,i] - A[j,i - 1])/2
    if y:
        #handle first row, where you can do central diff, and bottomost row
        for y in range(A.shape[0]): #rows indexed
            for x in range(A.shape[1]): #columns indexed
                if y == 0:
                    fill_me[y,x] = (A[y+1,x] - A[y,x])/2 #first row
                    
                elif (y != 0) and (y != A.shape[1] - 1): #central diff applicable
                    fill_me[y,x] = (A[y+1,x] - A[y-1,x])/2
                    
                elif y == A.shape[1] - 1: #bottom row
                    fill_me[y,x] = (A[y,x] - A[y-1,x])/2

    return fill_me

def electric_field_plotter(xxx,yyy,x_comp,y_comp):
    """
    Given x and y cordinate arrays, the x & y components,
    make a streamplot of the electric field.
    """
    fig,ax = plt.subplots()

    plt.streamplot(xxx,yyy, x_comp,y_comp, color = 'k', linewidth=1.5, arrowsize=1)
    
    ax.set_title('Electric field as function of position')
    #remove axis labels
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.tick_params(bottom = False)
    ax.tick_params(left = False)
    
    plt.show()