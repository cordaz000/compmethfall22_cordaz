#bring in packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#Part a) plots x = 2 cos(theta) + cos (2 theta), y = 2 sin (theta) -2 sin (2 theta)
#with theta being in range 0 <= theta <= 2 pi

#create a bunch of values in given range
theta1 = np.linspace(0,2*np.pi, 100) #radians, as np.sin and np.cos like

#function to remove tick and axis labels
def tick_label_remover(ax):
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.tick_params(bottom = False)
    ax.tick_params(left = False)


#initialize canvas
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3)

#each independent and dependent variable when plotting will be followed
#by the number of its corresponding axis above, e.g ax1 has x1, y2, etc.

#values for first plot
x1 = 2*np.cos(theta1) + np.cos(2*theta1)
y1 = 2*np.sin(theta1) - np.sin(2*theta1)
ax1.set_xlabel('')
ax1.set_ylabel('')
ax1.set_title('Deltoid Curve')
tick_label_remover(ax1)

ax1.plot(x1,y1)

#make the second plot, r = (theta)**2
theta2 = np.linspace(0,10*np.pi, 100) #radians
r1 = pd.Series([theta**2 for theta in theta2], name  = 'r')
dheta = pd.Series(data = theta2, name = 'theta')

#calculate x and y using values of r & theta
x2 = pd.Series([r*np.cos(theta) for r,theta in zip(r1,dheta)], name = 'x')
y2 = pd.Series([r*np.sin(theta) for r,theta in zip(r1,dheta)],name = 'y')
x_v_y = pd.DataFrame(data = [x2,y2]).transpose() #tabular data x vs y
ax2.set_title('Galilean Spiral')
tick_label_remover(ax2)


ax2.plot(x2,y2)


#plot the last equation r = exp( np.cos(theta) - 2*np.cos(4*theta))
# + (sin(theta/12))**5
theta3 = np.linspace(0, 24*np.pi, 1000)

r3 = np.array([np.exp(np.cos(theta)) - 2*np.cos(4*theta) + \
    (np.sin(theta/12))**5 for theta in theta3])

x3 = [r*np.cos(theta) for r,theta in zip(r3, theta3)]
y3 = [r*np.sin(theta) for r, theta in zip(r3, theta3)]
ax3.plot(x3,y3)
ax3.set_xlabel('')
ax3.set_ylabel('')
ax3.set_title("Fey's Function!")
tick_label_remover(ax3)
plt.tight_layout()
plt.show()