import numpy as np
from random import seed,random
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import time as Time
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D


def spatial_grid_and_time_steps(num_forecast_days, spat_res_x, spat_res_y,\
                      t_step_hrs, mean_zonal_wind, mean_height,EarthRadius,
                      fracR_forSpatialDomain):
    """
    Set up the initial flow. Provide the number of forecast days.
    Provide the number of grids desired in x and y.
    Provide the fraction of hour you'd like to time step.
    Provide an average zonal wind.
    Provide a mean height at which you'd like to solve BVE
    Input the Radius of the Earth and it's fraction over which
    you wish to set-up as the spatial domain.
    Return coordinate matrices xx and yy (from meshgrid), along
    with the total number of time-steps the loop will do.
    """
    
    #specificy initial conditions for the flow
    fcast_days = num_forecast_days #how long you want your forecast to be
    nx = spat_res_x ####----------------Used as input below-----------########
    ny = spat_res_y ####----------------Used as input below-----------########
    
    time_step_per_hour = t_step_hrs #how often are you time-stepping each hour
    u_avg = mean_zonal_wind #mean zonal flow 
    h_avg = mean_height #vertical level where we are to evaluate the flow
    
    #set up the constants and domain
    hours_2_seconds = 60*60 #(60 min/hour) * (60 seconds / min)
    seconds_in_one_day = 24*hours_2_seconds 
    seconds_in_forcast = fcast_days*seconds_in_one_day
    time_step = time_step_per_hour*hours_2_seconds
    nt = seconds_in_one_day/time_step #bound to set up time vector
    
    time_vector = np.linspace(0, nt, int(nt + 1))*time_step #second increments all the way up to 24hrs
    numberofTimeStepsperDay = seconds_in_one_day/time_step #time steps in a day
    tot_time_steps = int(seconds_in_one_day/time_step_per_hour) #total number of times time-loop will iterate---used below
    
    R_Earth = EarthRadius
    x_length = fracR_forSpatialDomain*R_Earth #set spatial scale at which to resolve problem
    y_length = fracR_forSpatialDomain*R_Earth
    delta_x = x_length/nx
    delta_y = y_length/ny
    
    #define grid
    x = np.linspace(0, nx, nx+1)*delta_x
    y = np.linspace(0, ny, ny+1)*delta_y
    print(x_length,nx, delta_x)
    x_mesh, y_mesh = np.meshgrid(x,y)
    
    xx = np.transpose(x_mesh)
    yy = np.transpose(y_mesh)
    
    return xx,yy, tot_time_steps,x_length,y_length

def GFD_constants(EarthRadius, latitudeOfInterest,mean_height):
    """
    Return Geophysical Fluid Dynamics constants of interest for this computation.
    Return a big tupple of these values
    
    """
    hours_2_seconds = 60*60 #(60 min/hour) * (60 seconds / min)
    Omega = (2*np.pi)/(24*hours_2_seconds) #Angular velocity of the Earth
    phi_not = np.deg2rad(latitudeOfInterest) #latitude where BVE will be solved
    
    #Coriolis and Beta parameter follow
    f_not = 2*Omega*np.sin(phi_not)
    beta_not = 2*Omega*np.cos(phi_not)/EarthRadius
    lil_g = 9.81 #units of m/s**2
    lambda_r = np.sqrt(lil_g*mean_height)/f_not #this is the Rossby radius of deformation
    kappa_d = 1/(lambda_r) # deformation wavevector, kappa_d
    return f_not, beta_not,kappa_d

def initial_conditions(spat_res_x, spat_res_y, coord_matrix_X, coord_matrix_Y,pertubation_depth,\
                      xlength, ylength,mean_height, f_not, beta_not,mean_zonal_wind):
    """
    Provided the spatial x and y resolution, a desired perturbation
    depth, the mean heigth at which BVE is being solved, along with 
    the Coriolis and Beta parameter, seet up the initial conditions 
    of the problem.
    This will return plots of your initial state and 
    """
    g = 9.81 # m/s**2
    z_not = np.zeros((spat_res_x+1,spat_res_y+1))
    depth = pertubation_depth
    xx = coord_matrix_X; yy = coord_matrix_Y
    
    perturbation = -depth*np.exp(  (-(xx/xlength - 0.2)**2 - (yy/ylength - 0.5)**2)/0.05 )
    z_not = z_not + perturbation
    
    #set first and last two columns of z_not to be the same
    z_not[:,spat_res_x] = z_not[:,spat_res_x - 1]
    z_not[:,0] = z_not[:,1]
    
    z_plus = z_not + mean_height
    
    height_hypsometric = (f_not/g)*mean_zonal_wind*yy #units of meters, should descend towards pole, increase 
    
    #towards the equator
    z_tot = z_plus - height_hypsometric
    
    #couple of plots to include initial configuration of system
    xp = xx/1000; yp = yy/1000 #so that you get km dimensions on axis
    
    plt.figure()
    CS = plt.contourf(xp, yp, z_not, cmap = 'plasma')
    plt.xlabel('x [km]',fontsize=15)
    plt.ylabel('y [km]',fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    cb=plt.colorbar()
    cb.set_label('Perturbation Height (m)',size=15)
    cb.ax.tick_params(labelsize=15)
    plt.title('500 hPa Geopotential Perturbation',fontsize=15) #initial Geopotential perturbation
    
    plt.figure()
    plt.contourf(xp, yp, z_tot, cmap = 'plasma')
    plt.xlabel('x [km]',fontsize=15)
    plt.ylabel('y [km]',fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    cb=plt.colorbar()
    cb.set_label('Total Geopotential Height (m)',size=15)
    cb.ax.tick_params(labelsize=15)
    plt.title('500 hPa Geopotential Height',fontsize=15) #initial Geopotential height
    
    psi_not = z_not #this has units of distance, this will turn into the stream function
    psi_tot = (z_not +  mean_height - height_hypsometric)* (g/f_not) #QG Stream function
    
    #plot stream function as surface
    fig = plt.figure()
    ax = plt.axes(projection = '3d')
    surface = ax.plot_surface(xp, yp, psi_not, cmap = 'plasma', linewidth = 1)
    fig.colorbar(surface, shrink = 0.5, aspect = 8, label = 'Geo. pot height w/respect to 500hPa')
    plt.title('Prescribed initial stream-function')
    plt.show(block = False)
    
    return psi_not

def main_time_loop(psi_not,time_steps, spatResX,spatResY,cordMatX, cordMatY,times2show,kappaDSqrd):
    """
    Return visualization of geopotential heights and zonal winds.
    Input initial streamfunction,desired time-steps, and desired
    spatial resolution --number of rows and columns wanted
    """
    #needed arrays
    nx = spatResX; ny = spatResY
    xx = cordMatX; yy = cordMatY
    xm = xx/1000; ym = yy/1000
    w = psi_not
    R = np.zeros((nx,ny))
    zonalwind = list()
    
    dwdx = np.zeros((nx+1,ny+1))
    dwdy = np.zeros((nx+1,ny+1))
    gradsq = np.zeros((nx+1,ny+1))
    d2wdx2 = np.zeros((nx+1,ny+1))
    d2wdy2 = np.zeros((nx+1,ny+1))
    laplac = np.zeros((nx+1,ny+1))
    dlapdx = np.zeros((nx+1,ny+1))
    dlapdy = np.zeros((nx+1,ny+1))
    wdldx = np.zeros((nx+1,ny+1))
    wdldy = np.zeros((nx+1,ny+1))
    dwdxl = np.zeros((nx+1,ny+1))
    dwdyl = np.zeros((nx+1,ny+1))
    dwdldydx = np.zeros((nx+1,ny+1))
    dwdldxdy = np.zeros((nx+1,ny+1))
    ddwdxldy = np.zeros((nx+1,ny+1))
    ddwdyldx = np.zeros((nx+1,ny+1))
    Jac1 = np.zeros((nx+1,ny+1))
    Jac2 = np.zeros((nx+1,ny+1))
    Jac3 = np.zeros((nx+1,ny+1))
    Jarakawa = np.zeros((nx+1,ny+1))
    
    plt.figure()
    for n in range(1,time_steps):

    #Take derivatives using finite-difference method, set periodic boundary
    #conditions in x and symmetric in y.


    #x derivative of w
        dwdx[1:nx,0:ny+2]= (w[2:nx+1,0:ny+2]-w[0:nx-1,0:ny+2])/(2*Delta_x) #rows that you can do central diff
        dwdx[0,0:ny+2] = (w[1,0:ny+2]-w[nx-1,0:ny+2])/(2*Delta_x) #treat first row differently
        dwdx[nx,0:ny+2] = dwdx[0,0:ny+2] #set last row equal to first
        dwdx[0:nx+1,ny]=dwdx[0:nx+1,ny-1] #set last two columns equal to each other
        dwdx[0:nx+1,0]=dwdx[0:nx+1,1] #set first two columns equal to each other
    #y-derivative of w
        dwdy[0:nx+2,2:ny-1] = (w[0:nx+2,3:ny]-w[0:nx+2,1:ny-2])/(2*Delta_y);
        dwdy[0:nx+2,1]= (w[0:nx+2,2]-w[0:nx+2,1])/Delta_y
        dwdy[0:nx+2,ny-1]=(w[0:nx+2,ny-1]-w[0:nx+2,ny-2])/Delta_y
        dwdy[0:nx+2,0] =  dwdy[0:nx+2,1]
        dwdy[0:nx+2,ny] = dwdy[0:nx+2,ny-1];
    #Square of the gradient of w
        gradsq = dwdx**2+dwdy**2;
    #Second x-derivative of w
        d2wdx2[1:nx,0:ny+2] = (w[2:nx+1,0:ny+2]+w[0:nx-1,0:ny+2]-2*w[1:nx,0:ny+2])/(Delta_x**2)
        d2wdx2[0,0:ny+2] = (w[1,0:ny+2]+w[nx-1,0:ny+2]-2*w[0,0:ny+2])/(Delta_x**2)
        d2wdx2[nx,0:ny+2] = d2wdx2[0,0:ny+2]
        d2wdx2[0:nx+1,ny]=d2wdx2[0:nx+1,ny-1]
        d2wdx2[0:nx+1,0]=d2wdx2[0:nx+1,1]
    #Second y-derivative of w
        d2wdy2[0:nx+2,2:ny-1] = (w[0:nx+2,3:ny]+w[0:nx+2,1:ny-2]-2*w[0:nx+2,2:ny-1])/(Delta_y**2)
        d2wdy2[0:nx+2,1]=(w[0:nx+2,3]+w[0:nx+2,1]-2*w[0:nx+2,2])/(Delta_y**2)
        d2wdy2[0:nx+2,0] = d2wdy2[0:nx+2,1]
        d2wdy2[0:nx+2,ny-1]=(w[0:nx+2,ny-3]+w[0:nx+2,ny-1]-2*w[0:nx+2,ny-2])/(Delta_y**2)
        d2wdy2[0:nx+2,ny] = d2wdy2[0:nx+2,ny-1]
        laplac = d2wdx2+d2wdy2;
    #x-derivative of laplacian
        dlapdx[1:nx,0:ny+2] = (laplac[2:nx+1,0:ny+2]-laplac[0:nx-1,0:ny+2])/(2*Delta_x)
        dlapdx[0,0:ny+2] = (laplac[1,0:ny+2]-laplac[nx-1,0:ny+2])/(2*Delta_x)
        dlapdx[nx,0:ny+2] = dlapdx[0,0:ny+2]
        dlapdx[0:nx+1,ny]=dlapdx[0:nx+1,ny-1]
        dlapdx[0:nx+1,0]=dlapdx[0:nx+1,1]
    #y-derivative of laplacian
        dlapdy[0:nx+2,2:ny-1] = (laplac[0:nx+2,3:ny]-laplac[0:nx+2,1:ny-2])/(2*Delta_y)
        dlapdy[0:nx+2,1]= (laplac[0:nx+2,2]-laplac[0:nx+2,1])/Delta_y
        dlapdy[0:nx+2,ny-1]=(laplac[0:nx+2,ny-1]-laplac[0:nx+2,ny-2])/Delta_y
        dlapdy[0:nx+2,0] =  dlapdy[0:nx+2,1]
        dlapdy[0:nx+2,ny] = dlapdy[0:nx+2,ny-1];
        Jacobi = dwdx*dlapdy - dwdy*dlapdx
    #Compute the Arakawa Jacobian.
        Jac1 = Jacobi;
        wdldx = w*dlapdx
        wdldy = w*dlapdy

        dwdldydx[1:nx,0:ny+2] = (wdldy[2:nx+1,0:ny+2]-wdldy[0:nx-1,0:ny+2])/(2*Delta_x);
        dwdldydx[0,0:ny+2] = (wdldy[1,0:ny+2]-wdldy[nx-1,0:ny+2])/(2*Delta_x);
        dwdldydx[nx,0:ny+2] = dwdldydx[0,0:ny+2];
        dwdldydx[0:nx+1,ny]=dwdldydx[0:nx+1,ny-1]
        dwdldydx[0:nx+1,0]=dwdldydx[0:nx+1,1]

        dwdldxdy[0:nx+2,2:ny-1] = (wdldx[0:nx+2,3:ny]-wdldx[0:nx+2,1:ny-2])/(2*Delta_y)
        dwdldxdy[0:nx+2,1]= (wdldx[0:nx+2,2]-wdldx[0:nx+2,1])/Delta_y
        dwdldxdy[0:nx+2,ny-1]=(wdldx[0:nx+2,ny-1]-wdldx[0:nx+2,ny-2])/Delta_y
        dwdldxdy[0:nx+2,0] =  dwdldxdy[0:nx+2,1]
        dwdldxdy[0:nx+2,ny] = dwdldxdy[0:nx+2,ny-1];

        Jac2 = dwdldydx - dwdldxdy
        dwdxl = dwdx*laplac
        dwdyl = dwdy*laplac

        ddwdxldy[0:nx+2,2:ny-1] = (dwdxl[0:nx+2,3:ny]-dwdxl[0:nx+2,1:ny-2])/(2*Delta_y)
        ddwdxldy[0:nx+2,1]= (dwdxl[0:nx+2,2]-dwdxl[0:nx+2,1])/Delta_y
        ddwdxldy[0:nx+2,ny-1]=(dwdxl[0:nx+2,ny-1]-dwdxl[0:nx+2,ny-2])/Delta_y
        ddwdxldy[0:nx+2,0] =  ddwdxldy[0:nx+2,1]
        ddwdxldy[0:nx+2,ny] = ddwdxldy[0:nx+2,ny-1];

        ddwdyldx[1:nx,0:ny+2] = (dwdyl[2:nx+1,0:ny+2]-dwdyl[0:nx-1,0:ny+2])/(2*Delta_x);
        ddwdyldx[0,0:ny+2] = (dwdyl[1,0:ny+2]-dwdyl[nx-1,0:ny+2])/(2*Delta_x)
        ddwdyldx[nx,0:ny+2] = ddwdyldx[0,0:ny+2]
        ddwdyldx[0:nx+1,ny]=ddwdyldx[0:nx+1,ny-1]
        ddwdyldx[0:nx+1,0]=ddwdyldx[0:nx+1,1]

        Jac3 = ddwdxldy - ddwdyldx
        Jarakawa = (1/3)*(Jac1+Jac2+Jac3)
    #Use the energy and enstrophy preserving Jacobian.
        Jacobi = Jarakawa;

    #Compute the function to be stepped forward. The -F*w turn is to dampen out
    #small but quickly divergent modes
        Q_n = laplac - defWaveVec*w

    #First time through the loop:
        if n==1:
            Dt = Delta_t/2;
            Q_nm1 = Q_n;

            rmeshvec=np.linspace(0,nx-1,nx)
            smeshvec=np.linspace(0,ny-1,ny)
            [rmesh,smesh] = np.meshgrid(rmeshvec,smeshvec)

            rr = np.transpose(rmesh)
            ss = np.transpose(smesh)
            C_rs = 2*(np.cos(2*np.pi*rr/nx)-1)/Delta_x**2+2*(np.cos(2*np.pi*ss/ny)-1)/Delta_y**2- kappaDSqrd

    #Calculate the kinetic energy and enstrophy integrals, this allows us to see
    #how well energy is being conserved.
        Rgsq=np.zeros((nx,ny))
        Rgsq[0:nx,0:ny] = gradsq[0:nx,0:ny]
        Rwsq=np.zeros((nx,ny))
        Rwsq[0:nx,0:ny] = w[0:nx,0:ny]**2
        
        Rgsq[0:nx,0:ny] = laplac[0:nx,0:ny]
        Rwsq[0:nx,0:ny] = w[0:nx,0:ny]
        
        umax = abs(dwdy).max()
        vmax = abs(dwdx).max()
        maxx=np.array((umax,vmax))
        VMAX = maxx.max()
        Q_np1 = Q_nm1 - (2*Dt)*((lil_g/f_knot)*Jacobi + b_not*dwdx + Ubar*dlapdx)

    #Section 3.3: Solve the Helmholtz Equation (Del^2-F)w = R.Compute the fft of the
    #right hand side (strip off additional row and column).
        R[0:nx,0:ny] = Q_np1[0:nx,0:ny]
        R_hat = np.fft.fft2(R)

    #Compute the transform of the solution
        W_hat = R_hat/C_rs
    #Compute the inverse transform to get the solution at (n+1)*Delta_t.
        w_new = np.real(np.fft.ifft2(W_hat)) # We assume w is real
        w[0:nx,0:ny] = w_new
        w[nx,0:ny] = w[0,0:ny]     # Fill in additional column at east.
        w[0:nx+1,ny]=w[0:nx+1,0]   # Fill in additional row at north.

        #Add the term for the zonal mean flow.
        wtotal = w + Hbar - (f_knot/lil_g)*Ubar*yy # here is that YY 
        zonalwind.append(-lil_g*dwdy/f_knot)

    #Shuffle the fields at the end of each time-step
        Dt = Delta_t
        Q_nm1 = Q_n
        timeestep=n
       
        #Plot the new 500 mb geopotential height contour
        plt.clf()
        plt.contourf(xm, ym, wtotal)
        plt.colorbar()
        plt.title('500 mb Geopotential Height and Zonal Wind Contours')

        #Plot east-west wind strength contours
        plt.draw()
        CS = plt.contour(xm, ym, -lil_g*dwdy/f_knot,colors='k')
        plt.clabel(CS, inline=1, fontsize=10)
        plt.pause(.00001)
        
    return zonalwind

