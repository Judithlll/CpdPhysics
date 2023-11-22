import numpy as np
# from scipy.integrate import odeint
import sys
import cgs
import ode
import matplotlib.pyplot as plt
import disk_properties as dp 
import functions as f

class System(object):
    rout=27*cgs.RJ
    alpha=1e-4
    fraction=0.02 
    dgratio=0.01  #dust to gas ratio
    sigmamol=2e-15
    rhoint = 1.4
    deltaT=0.0
    mtot1 = 1e24 #total mass a single superparticle represents
    daction={}

    def __init__(self,Rdi,nini,time=0.0):

        #initialize parameter from txt file // disk.py


        self.Rdi=Rdi  #initial radius of particles
        self.nini=nini
        self.mini = 4*np.pi/3*self.rhoint *Rdi**3
        self.time=time  #initial time
        self.ntime = 0

        # define a disk class

        self.disk = Disk ()

        # define class of superparticles here
        self.particles = Superparticles(nini,self.mini,self.disk.rinn,self.disk.rout,self.mtot1)

        #the amount of solid mass that has crossed into the domain
        self.Minflux = 0
        self.Mcp=self.disk.Mcp_t(self.time)

    def Mcp(self,Mcp0=0.4*cgs.MJ):

        Mcp=self.disk.Mcp_t(self.time)
        return Mcp

    def update(self,tEnd):
        """
        How to evolving the gas in the disk is still not sure

        Evolving system to the self.time

        """

        #time derivative of particles
        dydtP = self.particles.dY2d_dt(self.particles.Y2d,self.time,self.disk)
            
        #timescale
        tscaleArr = np.abs(self.particles.Y2d/dydtP)
        deltaT = np.nanmin(0.2*tscaleArr)
        # import pdb; pdb.set_trace()

        if self.time+deltaT>tEnd:
            deltaT = tEnd - self.time

        #update particle properties
        Yt=self.particles.update(self.time,self.time+deltaT,self.disk)
        
        #TBD: find better way to integrate (try: scipy.integrate.quad)
        self.Minflux += self.disk.M_influx(self.time,self.time+deltaT)

        
        #post_process particles
        self.post_process()
        
        if 'remove' in self.daction.keys():
            #remove the particles from Y2d!
            
            self.particles.remove_particles(self.daction['remove'])
            self.nini-=len(self.daction['remove'])
            # import pdb; pdb.set_trace()

        if 'add' in self.daction.keys():
            self.particles.add_particles(self.daction['add'])
            self.nini+=self.daction['add']

        self.time += deltaT
        self.deltaT = deltaT
        self.ntime += 1
        
        return Yt

    
    def post_process (self):

        self.daction = {}

        loc = self.particles.Y2d[0]

        #particles that cross the inner disk edge
        idx, = (loc<self.disk.rinn).nonzero()
        if len(idx)>0:
            self.daction['remove'] = idx
        
        #particles that enter the domain
        Nadd = 0
        while self.Minflux>self.mtot1:
            Nadd += 1
            self.Minflux -= self.mtot1
        
        if Nadd>0:
            self.daction['add'] = Nadd

        #particles that are eaten by the planet
        # need to use pebble accretion rate 
        #....

        #particles that have become too big
        #....




        

#perhaps this class object is not necessary...
class Disk(object):

    def __init__(self):
        
        self.alpha = dp.alpha
        self.rout = dp.rout
        self.rinn = dp.rinn
        #self.t=t

    def Mcp_t(self,time):
        return dp.Mcp_t(time)
    
    def dotMg(self,time):
        return dp.dot_Mg(time)
    
    def dotMd(self,time):
        return dp.dot_Md(time)
    
    def M_influx(self,t0,tEnd):
        return dp.M_influx(t0,tEnd)

    def Sigmag (self,loc,time):
        
        return dp.Sigma_g(loc,time)
    
    def OmegaK(self,loc,time):
        return dp.Omega_K(loc,time)
    
    def Td(self,loc,time):
        return dp.T_d(loc,time)
    
    def cs(self,loc,time):
        return dp.c_s(loc,time)
    
    def vth(self,loc,time):
        return dp.v_th(loc,time)
    
    def Hg(self,loc,time):
        return dp.H_g(loc,time)
    
    def nu(self,loc,time):
        return dp.viscosity(loc,time)
    
    def rhog(self,loc,time):
        return dp.rho_g(loc,time)
    
    def lmfp(self,loc,time):
        return dp.l_mfp(loc,time)
    
    def vK(self,loc,time):
        return dp.v_K(loc,time)
    
    def eta(self,loc,time):
        return dp.eta(loc,time)
    
    def update(self,deltaT):
        pass


class Superparticles(object):
    rhoint=1.4
    pi=np.pi
    Sto=0.0001  
    error=1e-8


    def __init__(self,nini,mini,rinn,rout,mtot1):
        """
        systems initial properties

        nini: initial number of the particles
        mini: initial mass for every particles
        rinn: inner edge of the disk
        rout: outer edge of the disk
        """
        self.nini=nini
        self.mini=mini #physical mass
        self.rinn=rinn
        self.rout=rout
        self.mtot1 = mtot1

        ndim = 3# location // mass // total mass

        self.Y2d = np.empty((ndim,nini))
        self.Y2d[0] = 10**np.linspace(np.log10(rinn),np.log10(rout),nini)
        self.Y2d[1] = mini
        self.Y2d[2] = self.mtot1

    def dY2d_dt (self,Y2d,t,disk):
        """
        input:
            Y2d -- state vector
            time -- time
            disk -- disk object
        """

        #unpack the state vector
        r, mphy, mtot = self.Y2d   #maybe the total mass needn't to be put in Y2d

        Rd=(mphy/(self.rhoint*4/3*np.pi))**(1/3)

        eta=disk.eta(r,t)
        v_K=disk.vK(r,t)
        v_th=disk.vth(r,t)
        lmfp=disk.lmfp(r,t)
        rho_g=disk.rhog(r,t)
        Omega_K=disk.OmegaK(r,t)
        H_g=disk.Hg(r,t)
        dotMd=disk.dotMd(t)

        St,v_r = f.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rd)

        v_dd=np.abs(v_r)/2
        H_d=H_g*(1+St/disk.alpha*(1+2*St)/(1+St))**(-0.5)

            

        drdt = v_r
        #dR_ddt= v_dd*dot_M_d/4/pi**(3/2)/rho_int/H_d/r/v_r**2 *dr_dt

        sigD = dotMd /(-2*r*self.pi*v_r) #v_r<0
        
        dmdt=2*np.sqrt(np.pi)*Rd**2*v_dd/H_d*sigD  #eq. 5 of Shibaike et al. 2017

        Y2ddt = np.zeros_like(self.Y2d)
        Y2ddt[0] = drdt
        Y2ddt[1] = dmdt
        Y2ddt[2] = 0.0

        return Y2ddt 
    
    
    def update(self,t0,tFi,disk,nstep=10):
        """
        this integrate the particles until tFi
        -- d: disk object
        """

        tSpan=np.array([t0,tFi])
        tstep=(tFi-t0)/nstep #why 100?
    
        Y2copy = np.copy(self.Y2d)

        #integrates system to tFi
        Yt = ode.ode(self.dY2d_dt,Y2copy,tSpan,tstep,'RK5',disk)

        self.Y2d = Yt[-1,:,:]

        return Yt
    
    def remove_particles(self,remove_idx):
        self.Y2d = np.delete(self.Y2d, remove_idx, 1)

    def get_stokes_number(self,disk,t):
        r, mphy, mtot = self.Y2d
        Rd=(mphy/(self.rhoint*4/3*np.pi))**(1/3)

        eta=disk.eta(r,t)
        v_K=disk.vK(r,t)
        v_th=disk.vth(r,t)
        lmfp=disk.lmfp(r,t)
        rho_g=disk.rhog(r,t)
        Omega_K=disk.OmegaK(r,t)

        St,v_r = f.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rd) 
        return St,v_r      

    def add_particles(self,Nadd):

        ynew = np.array([[self.rout],[self.mini],[self.mtot1]])

        self.Y2d = np.append(self.Y2d, ynew,1)

        if Nadd!=1:
            print('[core]Error:can only add 1 particle // reduce timestep')
            sys.exit()

    
        #for moment put the implemented particles all at the outmost edge
        #lociL=np.linspace(self.rout,self.rout,add_number)
        #miL=np.linspace(self.mini,self.mini,add_number)
        #mtot=np.nanmean(self.Y2d[2])+add_number*self.mini
        #mtotL=np.linspace(mtot,mtot,add_number)
        #self.Y2d=np.append(self.Y2d,np.array([lociL,miL,mtotL]),1)
        #self.Y2d[-1]=mtot


        return