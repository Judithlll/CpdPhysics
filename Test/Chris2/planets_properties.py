import numpy as np
import cgs
import functions as f
## import core# disk should go
import physics
import disk_properties as dp

def planet_migration (gas,planetLoc,planetMass,time,rhopl):
    #aerodynamic drag migration
    out = gas.get_key_disk_properties (planetLoc, time)

    mcp = dp.Mcp_t(time)

    disk = physics.DISK (*out, planetLoc, time, mcp) #pro
    # disk.add_auxiliary ()
    disk.add_uservar (dp.user_add_var())    #variables
    disk.add_userfun (dp.user_add_fun())
    #disk = core.DISK (*out, planetLoc, time)
    #disk.add_auxiliary ()
    #disk.user_difined ()

    #CWO: only add what you need!
    #TBR
    if False:
        eta=disk.eta
        v_K=disk.vK
        v_th=disk.vth
        lmfp=disk.lmfp
        rho_g=disk.rhog
        Omega_K=disk.OmegaK
        dotMg=disk.dotMg
        Mcp=disk.Mcp
        Sigmag=disk.sigmaG
        cs=disk.cs

    mg = disk.mu*cgs.mp
    dotMg = dp.dot_Mg(planetLoc)
    
    ## CWO: Stokes number for planetesimals/planets a bit weird...
    #for the moment, I put it infinite
    if False:
        Rpl=(planetMass/(4/3*np.pi*rhopl))**(1/3)
        Stpl,vd=f.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rpl)

    vd = 0

    #Type I migration
    ## CWO: too messy and too long...
    qr = -0.14493*dotMg*cgs.gC*mcp*(-0.206349*planetLoc**(5.5)+planetLoc**4.5*disk.rout)/disk.alpha/planetLoc**(8.5)/disk.rout/np.sqrt(cgs.kB*(dotMg*cgs.gC*mcp/planetLoc**3/cgs.sigmaSB)**(1/4)/mg) #pressure gradient

    CI = 0.1
    bt1 = CI*(2.7+1.1*qr)   #a constant Ogihara 2014
    vt1 = bt1*(planetMass/mcp)*(disk.sigmaG*planetLoc**2/mcp)*(disk.v_K(planetLoc, time)/disk.c_s(disk.temp))**2 *disk.v_K(planetLoc, time)

    v_mig=vt1+vd
    return v_mig


class Planets(object):
    """
    get the planets
    """

    def __init__(self,disk,particles,loc_planets,m_initial):
        """
        initial parameters:

        m_initial: the initial mass of planets
        disk: the disk object
        particles: the superparticles object
        loc_planets: the initial location of planets
        time: the system time 
        """
        self.m=m_initial
        self.disk=disk
        self.particles=particles
        self.r=loc_planets
        self.t=0.0
        
        self.eta=disk.eta(self.r,self.t)
        self.v_K=disk.vK(self.r,self.t)
        self.v_th=disk.vth(self.r,self.t)
        self.lmfp=disk.lmfp(self.r,self.t)
        self.rho_g=disk.rhog(self.r,self.t)
        self.Omega_K=disk.OmegaK(self.r,self.t)    
       
        self.St=[]
        self.v_r=[]  

    def hill_radius(self):
        """
        get the hill radius of the seeds
        """
        hill_r=self.r*(self.m/3/self.disk.Mcp_t(self.t))
        return hill_r
    
    def get_effective_pebbles_index(self):
        """
        get the superparticle index inside of the planets accration scale
        """

        index_to_be_accreted=np.argwhere((self.particles.Y2d[0]<self.r+100*self.hill_radius())&(self.particles.Y2d[0]>self.r-100*self.hill_radius()))  
        return index_to_be_accreted 

    def get_pebbles_size(self):

        self.effective_particleP=self.particles.Y2d[:,index]
        
        self.RdL=(self.effective_particleP[1]/(self.particles.rhoint*4/3*np.pi))**(1/3)
        
        self.St,self.v_r=f.St_iterate(self.eta,self.v_K,self.v_th,self.lmfp,self.rho_g,self.Omega_K,self.RdL)


    def Sigma_p(self):
        """
        pebble surface density. 
        !!!Now it's Larger than Sigmag, abnormal!!
        """
        
        self.get_pebbles_size()
        sigmap=self.disk.dotMd(self.t)/2/np.pi/self.r/abs(self.v_r)
        return sigmap

    def pebble_scale_height(self):
        Hg=self.disk.Hg(self.r,self.t)
        Hp=Hg*(1+self.St/self.disk.alpha*(1+2*self.St)/(1+self.St))
        return Hp

    def critical_mass(self):
        self.M_critical=1/8*self.eta**3*self.St *self.disk.Mcp_t(self.t) #Shibaike 2019

    def Pebble_accretion_rate(self):
        """
        abnormal now
        
        different particles has different size and so different accretion rate, so maybe should change the total mass of particles???
        """
        index=self.get_effective_pebbles_index()
        if len(index)==0:
            print('WARNING: NO pebbles will be accreted, maybe some problems')
            P_eff=0
        else:
            self.get_pebbles_size()
            mus=self.m/self.disk.Mcp_t(self.t)
            hp=self.pebble_scale_height()/self.r

            delv_o_vK=0.52*(mus*self.St)**(1/3)+self.eta/(1+5.7*(mus/self.eta**3*self.St))
            
            P_eff=1/np.sqrt((0.32*np.sqrt(mus*delv_o_vK/self.St/self.eta**2))**(-2)+(0.39*mus/self.eta/hp)**(-2)) #Liu & Ormel 2018
        return P_eff
        
        
