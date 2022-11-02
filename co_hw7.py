#importance sampling integration


import numpy as np
from integrating_functions import trap_integral as TI #inputs are [a,b], N, integrand

def troubled_func(x):
    return 1/(x**(1/2))

N = 10000
#calucate the denominator of p(x)
denominator = round(TI(10e-5,1,N,troubled_func)) #round to get 2...can't put 0 as first arg

#When we multiply the divergeing function by denominator we get (1/2)*(x**0.5) --> p(x)
#next, we find a transformation formulat for generating random numbers between zero and one
#from this distribution

#starting from p(x)dx = q(z)dz, we get x = z**2, using p(x) from above, z is a random number [0,1]

def nice_func(x):
    return 1/(np.exp(x) + 1)

num_rand_points = 1000000
counter = 0
for i in range(num_rand_points):
    z = np.random.uniform(low = 0, high = 1)
    x = z**2 #this is the transformation function
    counter = counter + nice_func(x) #summation happening here

val_of_integral = (1/num_rand_points)*counter*denominator # this is (1/N)*Sum of g(x_i) * integral w(x)
print('\n')
print("The value of the integral is ", np.round(val_of_integral,2))
print('\n')