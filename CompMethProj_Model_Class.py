#barotropic_spectral.py rewrite from https://github.com/lmadaus/Barotropic-Python/blob/master/barotropic_spectral.py

import numpy as np
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyspharm as spharm
import os
from hyperdiffusion import del4_filter, apply_des_filter
import namelist as NL # <---- IMPORTANT! Namelist containing constants and other model parameters


class Model:
    """
    Class that stores and plots the flow fields for a constant density, 
    non-divergent (incompressible) fluid on a sphere and integrates them
    forward with the barotropic vorticity equation.
    """

    def __init__(self, ics, forcing = None):
        """
        This initializes the model. It needs:
        -ics: A dictionary of linearized fields and space and time dims,
              its keys are : u_bar, v_bar, u_prime, v_prime, lats, lons and a start time.

        -forcing: a 2D array with same shape as model fields that has a vorticity tendency [1/s^2]
                  that is imposed at each integration time step.
        
        """

        # 1) Store the space and time variables, i.e. dimensions
        self.lats = ics['lats']
        self.lons = ics['lons']
        self.start_time = ics['start_time']
        self.curtime = self.start_time

        # 2) Set up and save the non-divergent (incompressible) state
        # This begins to be mysterious: 'setting up the spherical harmonic transform object'
        self.s = spharm.Sharmt(self.nlons(), self.nlats(), rsphere=NL.Re,
                                gridtype = 'regular', legfunc = 'computed')

        #Truncation for the spherical transformation
        if NL.M is None:
            self.ntrunc = self.nlats()

        else:
            self.ntrunc = NL.M

        #Use the object--s-- to get initial conditions
        vortb_spec, div_spec = self.s.getvrtdivspec(ics['u_bar'], ics['vbar'])
        vortp_spec, div_spec = self.s.getvrtdivspec(ics['u_prime'], ics['v_prime'])
        div_spec = np.zeros(vortb_spec.shape) #Only want non-diverge

        #Reconvert to horizontal winds to get the non-divergent component of field
        self.ub, self.vb = self.s.getuv(vortb_spec, div_spec)
        self.up, self.vp = self.getuv(vortp_spec, div_spec)

        #Use these winds to get the streamfunction Psi, and the velocity potential Chi
        self.psib,chi = self.s.getpsichi(self.ub, self.vb) # mean streamfunction
        self.psip, chi = self.s.getpsichi(self.up, self.vp) # perturbation streamfunction

        #convert spectral vorticity to grid
        self.vort_bar = self.s.sepctrogrd(vortb_spec) #mean relative vorticity
        self.vortp = self.s.spectrogrd(vortp_spec) # perturbation relative vorticity

        # 3) Some more variables to store
        self.bmaps = create_basemaps(self.lons, self.lats)
        #get the vorticity tendency forcing (if any) for integration
        self.forcing = forcing

    #some simple dimensional functions
    def nlons(self):
        return len(self.lons)
    def nlats(self):
        return len(self.lats)

    def integrate(self):
        """
        Integrates the barotropic model using spherical harmonics.
        The configuration is set in the namelist.py file
        """

        lat_list_r = [x * np.pi/180 for x in self.lats]
        lon_list_r = [x * np.pi/180 for x in self.lons]

        #meshgrid
        lons, lats = np.meshgrid(self.lons, self.lats)
        lamb, theta = np.meshgrid(lon_list_r, lat_list_r)

        #Save for derivatives later
        dlamb = np.gradient(lamb)[1]
        dtheta = np.gradient(theta)[0]

        #plot initial conditions
        if NL.plot_freq != 0:
            self.plot_figures(0)

        #Now loop through the timesteps
        for n in range(NL.times):
            
            #Here's where it's at, the vorticity tendency is computed
            #compute tendency with beta as only forcing
            vort_tend = -2*NL.omega/ (NL.Re**2) * d_dlamb(self.psip + self.psib, dlamb) - \
                        Jacobian(self.psip + self.psib, self.vortp + self.vort_bar, theta, dtheta, dlamb)

            #hyperdiffusion can be applied for smoothing if required
            if NL.diff_opt == 1:
                vort_tend = vort_tend - del4_filter(self.vortp, self.lats, self.lons)
            elif NL.diff_opt == 2 :
                vort_tend = apply_des_filter(self.s, self.vortp, vort_tend, self.ntrunc, 
                            t = (n + 1) *  NL.dt / 3600).squeeze()
            
            #add any imposed vorticity tendency forcing
            if self.forcing is not None:
                vort_tend = vort_tend + self.forcing

            if n == 0: #first step do forward difference, vorticity at next step is vort + vort_tend * dt
                vortp_next = self.vortp + vort_tend * NL.dt
            
            else: # otherwise do the leapfrog
                vortp_next = vortp_prev + vort_tend * 2 * NL.dt

            
            #Invert this new vorticity to get the streamfunction Psi ( the u and v winds)
            #first, go back to spectral space
            vortp_spec = self.s.grdtospec(vortp_next)
            div_spec = np.zeros(np.shape(vortp_spec)) #Divergence is zero in the barotropic vorticity (incompressible)

            #Now use the spharm methods to update the u and v grid
            self.up, self.vp = self.s.getuv(vortp_spec, div_spec)
            self.psip, chi = self.s.getpsichi(self.up,self.vp)

            #Change vort_now to vort_prev
            #if it isn't the first step, add a filter ('Robert filter') to dampen out crazy modes
            if n == 0:
                vortp_prev = self.vortp
            else:
                vortp_prev = (1.-2.*NL.r)*self.vortp + NL.r*(vortp_next + vortp_prev) #what is the first part of this expression doing!?

            #update the vorticity
            self.vortp = self.s.spectogrd(vortp_spec)

            #update the current time
            cur_fhour = (n + 1)* NL.dt / 3600.
            self.curtime = self.start_time + timedelta(hours = cur_fhour)

             # make figure(s) every  "plot_freq" hours
            if NL.plot_freq != 0 and cur_fhour % NL.plot_freq == 0:
                #Go from Psi to geopotential
                print("Plotting hour", cur_fhour)
                self.plot_figures(int(cur_fhour))

    def plot_figures(self, n, winds = 'total', vorts = 'total', psis = 'pert', showforcing = True,
                    vortlevs = np.array([-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]))*1e-5,
                    windlevs = np.arange(20,61,4), hgtlevs = np.linspace(-500, 500, 26),
                    forcelevs = np.array([-15, -12, -9, -6, -3, 3, 6, 9, 12, 15])*1e-10):

        """
        Makes global and regional plots of flow.
        Inputs:
        n -> timestep number
        winds, vorts, psis -> which we plotting : mean, perturbed, or total field
        showforcing -> if this is True, contour the vorticity tendency forcing
        *levs -> contour/fill levels for the respective variables
        """

        #specify wind component(s) that are being plotted
        if winds == 'pert': u = self.up; v = self.vp
        elif winds == 'mean': u = self.ub; v = self.vb
        else:   u = self.up  + self.ub; v = self.vp + self.vb

        #specify vorticity plotting
        if vorts == 'pert': vort = self.vortp
        elif vorts == 'mean' : vort = self.vort_bar
        else:               vort = self.vortp + self.vort_bar
        
        #what streamfunction are we plotting
        if psis = 'pert': psi = self.psip
        elif psis == 'mean': psi = self.psib
        else:               psi = self.psip + self.psib

        #make global zeta and wind barb map
        fig, ax = plt.subplots(figsize = (10,8))
        fig.subplots_adjust(bottom = 0.2, left = 0.05, right = 0.95)

        xx, yy = self.bmaps['global_x'], self.bmaps['global_y']
        cs = ax.contourf(xx,yy, vort, vortlevs, cmap=plt.cm.RdBu_r, extent = 'both')
        self.bmaps['global'].drawcoastlines()
        ax.quiver(xx[::2,::2], yy[::2,::2], u[::2,::2], v[::2,::2])

        #plot the forcing
        if showforcing and self.forcing is not None:
            ax.contour(xx, yy, self.forcing, forcelevs, linewidths = 2, colors = 'darorchid')
            ax.set_title('relative vorticity [s$^{-1}$] and winds [m s${-1}$] at %03d hours' % n)

            #colorbar
            cax = fig.add_axes([0.05, 0.12, 0.9, 0.03])
            plt.colorbar(cs, cax = cax, orientation = 'horizontal')
            #save figure
            if not os.path.isdir(NL.figdir + '/global'): os.makedirs(NL.figdir + '/global')
            plt.savefig('{}/global/zeta_wnd_{:03d}.png'.format(NL.figdir,n), bbox_inches = 'tight')
            plt.close()




