import numpy as np
import cgs
import functions as f
import core

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
        
        
