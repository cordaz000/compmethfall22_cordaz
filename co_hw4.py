#Fourth HW assignment

from scipy.constants import G # this has units of (m**3)(kg**-1)(s**-2)


import numpy as np


M_earth = 5.974*10**24 #units of kg
mass_moon = 7.348*10**22 # units of kg
R = 3.844*10**8 #units of m
omega = 2.662*10**-6 # units of Hz

#rearranging equation, we get G*M_earth*(R - r)**2 - G*mass_moon*r**2 - (omega**2)(r**3)*(R - r)**2 = 0
# above is f(r) = 0, we want to find r

#we calcuate derivative by hand of f(r) and get
# f_prime = 2*G*M_earth*(R - r) - 2*G*mass_moon*r + (omega**2)*(r**2)*(R - r)*(2*r - 3*(R - r))

#implement Newtown's method


r = R*(0.75) #initial guess for radius, 75% of distance from Earth to Moon

delta = 1*10**3 #initialize second term of equation (6.96) from textbook

accuracy = 1e-5 #accurate to at least 4 significant digits !?

i = 0 # counter, to not get stuck in a while loop forever if guess if off
while abs(delta) > accuracy:
    f = G*M_earth*(R - r)**2 - G*mass_moon*r**2 - (omega**2)*(r**3)*(R - r)**2 # this is f(r) set to zero

    f_prime = 2*G*M_earth*(r - R) - 2*G*mass_moon*r + (omega**2)*(r**2)*(R - r)*(2*r - 3*(R - r)) # divisor in second term

    delta = f/f_prime
    r = r - delta # equation (6.96) in textbook, implementation of Newton's method
    i = i + 1
    if i >= 100:
        print("One-hundred guesses, off.")
        break

r_km = r/1000
print("The Lagrangian point (L1) is approximately " + str(np.round(r_km,4)) + " km away from Earth.")
print('\n')
