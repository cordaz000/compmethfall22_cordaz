#make functions to integrate
import numpy as np

#trapezoidal rule integrator
def trap_integral(a,b,N,f):
    """
    Provided initial point a, final point b,
    number of bins N and analytic function f;
    find the approximation to the area under the 
    curve as per the trapezoidal rule.
    """
    #find bin width
    h = (b - a) / N
    s = 0.5*f(a) + 0.5*f(b) #first and last heights get counted once

    for k in range(1,N):
        s = s + f(a + k*h) #intermediate heights get double counted
    
    approx_area = h*s # bin width (regular spaced)  multiplied by heights

    return approx_area

def simpson_integral(a,b,N,f):
    """
    Provided initial point a, final point b, 
    number of bins N and analytic function f,
    find the approximation to the area under the curve
    using Simpson's rule.
    """

    h = (b - a)/N #uniform bind width
    
    #handle odd terms
    s_odd = 0
    for k in range(1,N - 1,2):
        
        s_odd = (s_odd + f(a + k*h) )

    #handle even terms
    s_even = 0
    for k in range(2,N - 2,2):
        s_even = (s_even + f(a + k*h))

    return h*( (f(a)+f(b))/3 + (4/3)*s_odd + (2/3)*s_even)

def integral_area_function_x(delxs, integrator_func, func):
    """
    Provided a list with tuples of intervals (in this case spaced
    with del_x = 0.1, thirty of them), a function to calculate the 
    integral via Trapezoid rule and a function to integrate, 
    return an array as long as the # of tuples having the area of 
    all previous trapezoids summed together, i.e, the last element in
    the array will be the total area of the integral over interval 
    specified.
    """

    area_holder = list()
    
    for tup in delxs:
        area_holder.append(integrator_func(tup[0], tup[1], 1, func))
    
    area_holder = np.array(area_holder)


    return np.array( [ 0 if i== 0 else np.sum(area_holder[:i]) for i in range(len(area_holder)) ])