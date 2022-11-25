import numpy as np
import matplotlib.pyplot as plt

#definition of constants
g = 10.0 #units of m/s^2
rho = 1.22 # units of kg/ m^3
C = 0.47 #unitless coefficient of drag
radius_ball = 0.08 #in meters


def Second2FirstOrderODE(r,mass,t): # Coupled differential equations
    """
    Input a vector r with order of entries as follows:
    i)initial position x
    ii) initial velocity x
    iii)initial position y
    iv) initial velocity y
    t is the provided time step at which the function is evaluated.||
    """
    global g; global rho; global C;
    global radius_ball;
    drag_const = (np.pi*(radius_ball**2)*rho*C)/(2*mass)
    pos_x = r[0]
    vel_x = r[1]
    pos_y = r[2]
    vel_y = r[3]
    fvx = vel_x
    fvy = vel_y
    fx = -drag_const*vel_x*np.sqrt(vel_x**2 + vel_y**2)
    fy = -g - drag_const*vel_y*np.sqrt(vel_x**2 + vel_y**2)
    
    return np.array([fvx,fx,fvy,fy],float)

def Runge_Kutta_4th_Order(tPoints,inputVect, stepSize,\
    horizPoints,vertPoints,func,SmallParam,mass):

    """
    This RK4th needs:
    -list with time intevals to run RK 4th order recipe
    -Input vector to hold initial conditions
    -step size to do the RK 4th order algorithm
    -two lists to contain horizontal points to plot trajectory
    -function to turn second order ODE to 1st order
    -Small parameter to stop algorithm (comparison)
    It returns the two lists with x and y coordinates of trajectory
    """
    for t in tPoints:
        horizPoints.append(inputVect[0])
        vertPoints.append(inputVect[2])

        k1 = stepSize*func(inputVect,mass,t) #check

        k2 = stepSize*func(inputVect+ 0.5*k1 ,mass,t + 0.5*stepSize) #check

        k3 = stepSize*func(inputVect+ 0.5*k2 ,mass,t + 0.5*stepSize) #check 

        k4 = stepSize*func(inputVect + k3 ,mass, t + stepSize)

        inputVect = inputVect + (k1 + 2*k2 + 2*k3 + k4)/6

        if inputVect[2] <= SmallParam:
            print("The total distance travelled with " + str(np.round(mass,2)) + " kg ball is about " \
                + str(np.round(inputVect[0],1)) + " meters.")
            break

    return horizPoints,vertPoints


def TrajectoryPlotter(xcoords,ycoords):
    """
    Plot trajectory of a single Cannonball.
    """
    fig,ax = plt.subplots()
    ax.plot(xcoords,ycoords,'ko' ,linewidth=1.5)
    ax.set_ylim(-0.5,57)
    ax.set_ylabel('Height [m]')
    ax.set_xlabel('Distance [m]')
    ax.set_title('Cannonball Trajectory')
    