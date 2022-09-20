import numpy as np
import matplotlib.pyplot as plt


from integrating_functions import trap_integral
from integrating_functions import simpson_integral
from integrating_functions import integral_area_function_x


#integrate exp(-t**2)dt from 0, 3 in steps of 0.1
#I'll do this with both Trapezoidal and Simpson way and 
#estimate the error both ways

def g(t):
    return np.exp(-t**2) # function to integrate

#find how many bins you need to have h be 0.1
b = 3; a = 0; h = 0.1;
N = int((b-a)/h) #thirty steps

#first answer with Trapezoidal method
I_trap_1 = trap_integral(a,b,N,g)

#second way--double N and halve step size
I_trap_2 = trap_integral(a,b,2*N,g)


#calculate the error
Error_Trap = (1/3)*(I_trap_2 - I_trap_1)
print("The value of the integral using 1st order Newton-Cotes is " + str(round((I_trap_1),3))+ \
        " with an error of " + str(Error_Trap) + ".")
print('\n')

#repeat using Homer Simpson's integration scheme 
I_simp_1 = simpson_integral(a,b,N,g)
I_simp_2 = simpson_integral(a,b,2*N,g)

Error_Simp = (1/15)*(I_simp_2 - I_simp_1)

print("Value of integral using 2nd order Newton-Cotes integration is " \
        + str(np.round(I_simp_1,3)) + " and the error is " + str(Error_Simp))

intervals = np.linspace(0,3, 31)

dxs = list()
for i in range(1,len(intervals)):
    dxs.append((intervals[i-1], intervals[i])) # [ (0, 0.1), (0.1, 0.2), etc]

y_vals = integral_area_function_x(dxs,trap_integral, g)

#plot the value of the integral as a value of x

fig,ax = plt.subplots()

ax.scatter(intervals[1:],y_vals, color = 'k')
ax.set_ylabel('Area as function of x')
ax.set_title('Value of the integral as a function of x')
plt.show()
