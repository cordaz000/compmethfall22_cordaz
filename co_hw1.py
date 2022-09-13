#bring in packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


#Part a) plots x = 2 cos(theta) + cos (2 theta), y = 2 sin (theta) -2 sin (2 theta)
#with theta being in range 0 <= theta <= 2 pi

#create a bunch of values in given range
theta1 = np.linspace(0,2*np.pi, 100) #radians, as np.sin and np.cos like

#initialize canvas
fig, (ax1,ax2, ax3) = plt.subplots(nrows = 1, ncols = 3)

#values for first plot
x = 2*np.cos(theta1) + np.cos(2*theta1)
y = 2*np.sin(theta1) - np.sin(2*theta1)
ax1.set_title('Deltoid Curve')
ax1.plot(x,y)

#make the second plot, r = (theta)**2
theta2 = np.linspace(0,10*np.pi, 100) #radians
r1 = pd.Series([theta**2 for theta in theta2], name  = 'r')
dheta = pd.Series(data = theta2, name = 'theta')

#calculate x and y using values of r & theta
x1 = pd.Series([r*np.cos(theta) for r,theta in zip(r1,dheta)], name = 'x')
y1 = pd.Series([r*np.sin(theta) for r,theta in zip(r1,dheta)],name = 'y')
x_v_y = pd.DataFrame(data = [x1,y1]).transpose() #tabular data x vs y
ax2.set_title('Galilean Spiral')
ax2.plot(x1,y1)
plt.tight_layout()
plt.show()