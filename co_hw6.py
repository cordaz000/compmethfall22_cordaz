#random walk HW

import matplotlib.pyplot as plt 
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import random

x = 0
y = 0

x_data=[]
y_data=[]

fig, ax = plt.subplots()
ax.set_xlim(-10,10)
ax.set_ylim(-10,10)
line, = ax.plot(0,0)
counter = list()
def animate_rw(i):
    global x # kept getting a "local variable referenced before assignment error
    global y
    global counter
# No need for the move_random function
    direction = random.randint(1, 4)
    if direction == 1:
        x += 1
    elif direction == 2:
        y += 1
    elif direction == 3:
        x += -1
    elif direction == 4:
        y += -1

    x_data.append(x)
    y_data.append(y)

    line.set_xdata(x_data)
    line.set_ydata(y_data)

    return line,

N = 100 # number of times animate function gets called--this was supposed to be a million, too high
# Did not write in anything for frames, since it defaults to passing itertools.count
anim = FuncAnimation(fig, animate_rw, frames = N, interval = 600) #frames controls how many random walks your program makes

plt.show()