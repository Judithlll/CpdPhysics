from os import wait
import numpy as np
from scipy.integrate import odeint
import sys
import cgs
import ode
import matplotlib.pyplot as plt
import disk_properties as dp 
import functions as ff
import copy
import parameters as pars
from gas import GAS
import scipy.optimize as sciop
import scipy.integrate as sciint
import time
import scipy.optimize as sciop
import physics
import userfun
import resample

class COPY (object):
    """
    This object used to back up the system state of last time point
    """

    def __init__ (self, state, attributeL):

        for attr in attributeL:
            setattr(self, attr, copy.deepcopy(getattr(state,attr)))

class Mintimes (object):
    """
    store the minimam times
    """
    def __init__(self, mintimeL, jumpfracD={}):
        tevol = []
        tminL = []
        nameL = []

        for i in range (len(mintimeL)):
            tmin = mintimeL[i]['tmin']
            name = mintimeL[i]['name']
            setattr(self,name , tmin)
            
            tminL.append(tmin) 
            nameL.append(name)
            #TBD:need to be more general 
            if i>0:#don't consider the particles timescale here
                if mintimeL[i]['name'] == 'PlanetsRes':
                    tevol.append( jumpfracD['PlanetsRes']*tmin )
                else:
                    tevol.append( jumpfracD['general']*tmin)
        
        self.tminarr = np.array(tminL)

        self.nameL = nameL

        if len(tevol)>0:
            self.max_tevol = np.min(tevol)

        #first element is particles
        self.dpart = mintimeL[0]



class Messages (object):
    """
    this class collects messages (experimental)
    """
    def __init__ (self):
        self.msgL = [] #the message list

    def add_message (self, ntime, mtype, msg):
        dmsg = {'ntime':ntime, 'type':mtype, 'msg': msg}
        self.msgL.append(dmsg)

    def flush (self):
        if len(self.msgL)>0:
            sfmt = '{:8d} {:12s} {:s}'
            with open ('log/messages.log', 'a') as file:
                for dmsg in self.msgL:
                    line = sfmt.format(dmsg['ntime'],dmsg['type'],dmsg['msg'])
                    file.write(line+'\n')

            self.msgL = []#reset


class System(object):

    """
    SYSTEM: integrate the every part of the PPD (or CPD) and update them with time.
    """

    ## TBD-later: think about how to calculate rhoPlanet

    daction={}

    def __init__(self, timeini=0.0, rhoPlanet=1.9):
        """
        System initiallization parameters:
            
            timeini: [float]
                the initial time with default value 0.0(maybe not necessary)
            rhoPlanet: [float]
                internal density of planets with default value 1.9g/cm^3
        """
        self.rhoPlanet = rhoPlanet

        self.time=timeini
        self.ntime = 0  
        self.timeL = []
        # define a disk class
        self.gas = self.init_gas ()
        self.init_centralbody()

        #the amount of solid mass that has crossed into the domain
        self.dotMg = dp.dot_Mg(self.time)
        self.Minflux = 0
        self.Minflux_step = 0

        self.Moutflux = 0.

        #initiallize the old state
        self.oldstate=None

        #self.planetMassData=[]
        self.njump = 0
        self.njumptime = 0
        self.jumptime = 0.0

        self.resam_time = np.array([])

        self.milestones = {}

        self.doJump = False

        #some resonance information
        self.dres = self.make_resL (jmax=10)
        self.res_setL = []

        self.nplanet = 0 #start with 0 planets (??)
        
        #TBD:-later: make this more general
        #   like: r_crit={'cavity':...}
        self.rinn = self.get_rinn()
        self.rout = self.get_rout()

        #[24.05.12]cwo: create a messaging class
        self.messages = Messages ()

        self.con_Jump =[]

        self.remidx = None


    def add_message (self, mtype, message):
        """
        just passes the message on...
        """
        self.messages.add_message (self.ntime, mtype, message)


    def get_auxiliary (self, time= None):
        '''
        Do the get_auxiliary of particals
        the disk need to be updated according to the time, and so the 
        mcp 
        '''

        disk = self.get_disk(time)
        self.specloc = np.append([i.loc for i in self.icelineL], [p.loc for p in self.planetL])
        if self.specloc is None: 
            import pdb;pdb.set_trace()
        self.particles.get_auxiliary(disk, time, specloc=self.specloc)


    def re_sample (self):

        mphy = self.particles.massL
        if np.all(np.diff(np.log10(mphy[10:1000]))<0)==False and False:
            print('physical mass not in order')
            import pdb; pdb.set_trace()

        if pars.resampleMode=='splitmerge':

            #the particles crossing is not physical
            #try:
            #    assert( np.all(np.diff(self.particles.locL)>0.) )
            #except:
            #    print('[core.re_sample]: Particles are crossing, please check')
            #    import pdb;pdb.set_trace()

            newarr = resample.re_sample_splitmerge(self,self.particles,nsampleX=2,full_output=True, **pars.dresample)
        elif pars.resampleMode == 'dropmerge':
            newarr = resample.re_sample_dropmerge(self,self.particles,nsampleX=2,**pars.dresample)
        
        elif pars.resampleMode == 'global_resample':
            #LZX [24.11.05]: this is not stable now
            newarr = resample.global_resample(self, self.particles, **pars.dresample)

        #[25.01.01]cwo: variation on the above
        elif pars.resampleMode == 'global_resample2':
            newarr = resample.global_resample2(self, self.particles, **pars.dresample)

        #[25.01.13]lzx: just for test 
        elif pars.resampleMode == 'global_resample3':
            newarr = resample.global_resample3(self, self.particles, **pars.dresample)

        #[25.01.01]cwo: another variation...
        elif pars.resampleMode == 'global_resample4':
            newarr = resample.global_resample4(self, self.particles, **pars.dresample)

        #[25.01.21]cwo: fixed_resampling (should be similar to global_resample)
        elif pars.resampleMode == 'fixed_resample':
            newarr = resample.fixed_resample(self, self.particles, self.specloc, **pars.dresample)

        elif pars.resampleMode == 'local_splitmerge':
            newarr = resample.local_splitmerge(self, self.particles, **pars.dresample)

        #[25.01.18]cwo: another variation...
        elif pars.resampleMode == 'new_splitmerge_chris':
            newarr = resample.new_splitmerge_chris(self, self.particles, **pars.dresample)
        elif pars.resampleMode == 'new_splitmerge_zxl':
            newarr = resample.new_splitmerge_zxl(self, self.particles, full_output=True, **pars.dresample)

        else:
            newarr = None

        if newarr is not None:
            #the global resample will sometimes have some mass loss so the Moutflux should be updated here 
            if pars.resampleMode == 'global_resample':
                self.Moutflux += np.sum(self.particles.mtotL) - np.sum(newarr[1])

            self.resam_time= np.append(self.resam_time, self.time)

            if np.all(np.diff(newarr[0])>0.):

                # if self.time>17*cgs.yr: 
                # #[24/11/11]let's compare the properties before and after the resampling here 
                #     mtotn = newarr[1]
                #     locn = newarr[0] 
                #     sfdn = ff.sfd_simple(mtotn, locn)
                #
                #     plt.figure()
                #     plt.xscale('log')
                #
                #     plt.xlabel('location')
                #     plt.ylabel('sfd')
                #     plt.plot(locn/cgs.RJ, sfdn, '.-', c='r', label='after')
                #     plt.plot(self.particles.locL/cgs.RJ, self.particles.sfd, '.-', c ='b', label='before')
                #     plt.legend()
                #     plt.savefig('sfd_ba_resample.png')
                #     plt.close()
                #
                #     plt.figure()
                #     plt.xlabel('location')
                #     plt.ylabel('cumsum')
                #     plt.loglog(locn/cgs.RJ, np.cumsum(mtotn), '.-', c='r', label='after')
                #     plt.loglog(self.particles.locL/cgs.RJ, np.cumsum(self.particles.mtotL), '.-', c ='b', label='before') 
                #     plt.legend()
                #     plt.savefig('cumsum_ba_resample.png')
                #     plt.close()
                #     import pdb;pdb.set_trace()


                #assign the key properties 
                self.particles.locL,self.particles.mtotL,self.particles.massL,self.particles.fcomp = newarr[0:4]
                self.particles.num = len(self.particles.locL)
                #also we need to get other properties of particles 
                #idea is make an function to get all these auxiliary
                self.get_auxiliary(self.time)

                if np.isnan(self.particles.locL).any():
                    import pdb; pdb.set_trace()
                

                isL, imL = newarr[4:6]
                # the index of new particles 
                self.remidx = isL[0]+1 
                pars.resampleMode = None

            else:
                print('some crossing caused by resampling')

        if self.remidx is not None: 
            if 'remove' in self.daction:
                self.remidx -= len(self.daction['remove']) 

            userfun.check_split(self, self.remidx)
            import pdb;pdb.set_trace()
            #plot the mphy, vr, sfd, around the new particle 
                
                

    #ZL-TBD: put this stuff in fileio.py
    #[24.04.21]CWO: I've again removed the logdir. It should ALWAYS be local!       
    def update_log(self, djump={}, logdir = './log/', init=False):
        import subprocess as sp

        loglist = ['system_evol.log', 'planets.log', 'jump.log']
        pathlist = [logdir+log for log in loglist]
        
        nameL = self.minTimes.nameL
        tminarr = self.minTimes.tminarr
        
        imin = tminarr.argmin() #minimum evolution
        nameM = nameL[imin]
        tmin = tminarr[imin]

        #give nameM more information
        if nameM == 'particles':
            #partpropL = ['drift','growth','collision']
            nameM = nameM + '-' + str(self.minTimes.dpart['imin'])

        # get the lines will be written into files
        efmt = '{:10d} {:10.2f} {:20d} {:10.2e} {:20s}'       
        jfmt = '{:10d} {:20.2f} {:20.2f} {:30s}'
    
        line_evol =efmt.format(self.ntime, self.time/cgs.yr, self.particles.num, self.deltaT, nameM.rjust(20))+''.join('{:7s}'.format(str(con).rjust(7)) for con in self.con_Jump)

        #line_evol =efmt.format(self.ntime, self.time/cgs.yr, self.particles.num, self.Moutflux/self.minitDisk, nameM.rjust(20)) 

        #[24.05.12]cwo:please, don't print to screen here..
        #print (line_evol)
        line_plan = '{:10.2f}'.format(self.time/cgs.yr)+''.join('{:10d} {:20.2f} {:20.2f} {:20.2e} {:10.5f}'.format(self.planetL[i].number, self.planetL[i].max_jumpT, self.planetL[i].loc/cgs.RJ,self.planetL[i].mass, self.planetL[i].acc_rate) for i in range(self.nplanet))
        #if self.nplanet >0:
        #    import pdb;pdb.set_trace()

        # get the line of jump
        if self.doJump:
            nameJ = djump['tjumpkeys'][np.argwhere(djump['tjumparr']== self.jumpT)[0][0]] 
            #TBD: also need to know which resonance and which system_evol
            if nameJ == 'milestones':
                nameJ = nameJ + '-'+ self.milestones[self.time]

            self.jump_limitation = nameJ

            line_jump = jfmt.format(self.njump, self.jumptime/cgs.yr, self.jumpT/cgs.yr, nameJ.rjust(30))

            print(f'[core.system.jump]:at {self.timeL[-2]/cgs.yr:5.2f} yr jumped by {self.jumpT/cgs.yr:5.2f} yr')
            print(f'[core.system_jump]:jump time limited by: '+nameJ)
            print(f'[core.system_jump]:min. evolution time ({nameM}) {tmin/cgs.yr:9.2e} yr')
            #import pdb;pdb.set_trace()
            lines = [line_evol, line_plan, line_jump,]
        else:
            lines = [line_evol, line_plan]


        # write to files
        if init:
            #[24.05.26]cwo: it is cleaner to put stuff like "clear_dir" in separate fileio module
            ff.clear_dir(logdir)

            #[24.05.26]I also initialize message file
            sp.run(['touch', logdir+'messages.log'])

            # generate the title foemation
            etfmt = '{:10s} {:10s} {:20s} {:10s} {:20s} {:7s}'
            jtfmt = '{:10s} {:20s} {:20s} {:30s}'
            ptfmt = '{:10s} {:10s} {:20s} {:20s} {:20s} {:10s}'
            line_evol_title = etfmt.format('ntime'.rjust(10), 'time'.rjust(10), 'particles_number'.rjust(20), 'deltaT'.rjust(10), 'restrict_factor'.rjust(20), 'conditions'.rjust(7))
            line_jump_title = jtfmt.format('njump'.rjust(10), 'time_begin_jump'.rjust(20), 'time_jumped_over'.rjust(20), 'jump_retrict_factor'.rjust(30))
            line_plan_title = ptfmt.format('time'.rjust(10), 'planet_num'.rjust(10), 'planet_max_jumpT'.rjust(20), 'location'.rjust(20), 'mass'.rjust(20), 'acc_rate'.rjust(10))
            titles = [line_evol_title, line_plan_title, line_jump_title]

            for i,p in enumerate(pathlist):
                with open (p, 'w') as file:
                    file.write(titles[i] + '\n')
        else:
            for i,line in enumerate(lines):
                with open (pathlist[i], 'a') as file:
                    file .write(line + '\n')

        #[24.05.26]cwo: added writing to messages
        self.messages.flush()
            

    def get_rinn(self):
        """
        get the rinn from the disk_properties or the central object's radius
        """
        try:
            rinn = dp.rinn
        except:
            rinn = self.centralbody.r

        return rinn

    def get_rout(self):
        """
        get the rout from the disk_properties.py if it exists 
        """
        try:
            rout = dp.r_out(self.time)
        except:
            rout = None 

        return rout

    def make_resL (self,jmax=4):
        dres = {}
        dres['res'] = []
        dres['prat'] = []
        for j in range(1,jmax):
            sj = str(j)
            sj1 = str(j+1)
            dres['res'].append(sj1+':'+sj)
            dres['prat'].append((j+1)/j)

        #final one (hack)
        dres['res'].append('inf')
        dres['prat'].append(1.0)

        dres['prat'] = np.array(dres['prat'])
        return dres

    ## CWO: let's add planet by the system like this..
    #   This is copied from another program from mine
    def add_planet (self, planet):
        self.planetL.append(planet)
        planet.number = self.nplanet
        self.nplanet += 1
        self.planetD.setdefault(planet.number, planet)
        locpl = planet.loc
        
        #sort the planetL (TBD)

        locL = [self.planetL[i].loc for i in range(self.nplanet)]
        iloc = locL.index(locpl)#provides the new sorted location of the added planet

        #LZX [24.08.04]: seems we defaultly think that the period ratio between planets will decrease, is that right?
        #there is an interior planet
        if iloc>0:
            prat = (locpl/self.planetL[iloc-1].loc)**1.5
            inxt = (self.dres['prat']<prat).argmax()
            planet.inxt = inxt  #next upcoming res.
            planet.resS = -1    #-1: not yet in resonance

        #there's an exterior planet
        if iloc<self.nplanet-1:
            prat = (planetL[iloc+1]/locpl)**1.5
            inxt = (self.dres['prat']<prat).argmax()
            planetE = self.planetL[iloc+1]
            planetE.inxt = inxt  #next upcoming res.
            planetE.resS = -1    #-1: not yet in resonance

    def init_centralbody(self):
        self.centralbody = CentralBody(self.time, dp.Mcp_t(0.0), cgs.rhoRJ)

    def init_particles (self, dparticleprops={}):
        """
        because we need to consider iceline, so separatly 
        initiallize the particles, for now just water iceline is considered  

        history:
        [23.12.30]CWO: instead of forwarding diskmass, I supply self.gas to the superparticle class
        [24.01.08]LZX: mtot1 is generated from Superparticles, but is used much 
                        in post_process, so get this from superparticles for now
        """

        #[24.07.21]CWO: not sure if we need dp.rout here...
        #               this self.rinn and dp.rout is weird...

        self.particles = Superparticles(self.rinn,self.rout,self.dcomposL,self.gas, **dparticleprops)
        self.get_auxiliary(self.time)

        self.minitDisk = sum(self.particles.mtotL)

        ##TBD: check the surface density
        #loc = self.particles.locL
        #sigG, *dum = self.gas.get_key_disk_properties(loc,0.0)

        ##LZX[24.08.06] uncommented for now, but maybe used later
        #sigP = sigG *0.01
        #mtot = self.particles.mtotL
        #boundaries = np.sqrt(loc[1:]*loc[:-1])
        #boundaries = np.append(dp.rinn,boundaries)
        #boundaries = np.append(boundaries,dp.rout)
        #warr = np.diff(boundaries)
        #sigma = mtot /(2*np.pi*loc*warr)
        #import pdb; pdb.set_trace()
        

    def init_gas (self, gasL=None, dcomposL=None, dgrid={}):
        """
        init_gas is called at initialization. It adds an instance of the GAS class

        history:
        [23.12.13]:this is copied from /NewLagrange code base...
                  :for the moment gasL and dcomposL are put to None

        """
        #also initializes a gas class (nested class)
        if pars.gasmodel=='gridstatic':
            dgas = pars.dgasgrid
        else:
            dgas = {}

        dum = GAS (gasL, dcomposL, mode=pars.gasmodel, time=self.time, **dgas)
        return dum


    def get_disk (self, time=None, loc=None):
        #TBD: generalize this..

        if time is None:
            time = self.time
            mcp = self.centralbody.m 
        else: 
            mcp = self.centralbody.get_mass(time)

        if loc is None:
            loc = self.particles.locL


        #we need these 2 things to initalize the class object
        out = self.gas.get_key_disk_properties (loc, time)

        disk = physics.DISK (*out, loc, time, mcp) #pro
        disk.add_auxiliary ()
        userparL = disk.add_uservar (dp.user_add_var())    #variables
        userfuncL = disk.add_userfun (dp.user_add_fun())    #functions only
        userevalL = disk.add_user_eval (dp.user_add_eval()) #evaluations

        disk.userprops = userparL+userfuncL+userevalL

        return disk


    def update_particles (self, timestepn=3, **kwargs):
        """
        Integrate the particles forward by amount deltaT

        Evolving system to self.time
        """
        # prevent the large value 'eats' the small value
        #if self.deltaT/self.time <1e-15:
        #    import pdb;pdb.set_trace()


        Y2d = self.particles.make_Y2d()
        Y0 = np.copy(Y2d)
        t0 = self.time
        tn = t0 +self.deltaT
        #NOTE: need to update the properties of particles according to the time of every 
        #      midium step

        if pars.dtimesteppars['itgmethod']=='Euler':
            Yn = Y0 +self.particles.dY2d_dt(Y0,t0) *self.deltaT

        
        elif pars.dtimesteppars['itgmethod']=='Heun':
            Y1 = Y0 +self.particles.dY2d_dt(Y0,t0) *self.deltaT

            self.particles.locL = Y1[0]
            self.particles.massL = Y1[1]
            self.get_auxiliary(tn)

            Y2 = Y0 +self.particles.dY2d_dt(Y1,t0) *self.deltaT
            Yn = 0.5*(Y1+Y2)

        elif pars.dtimesteppars['itgmethod']=='RK4':
            k1 = self.particles.dY2d_dt(Y0,t0)
            Y1 = Y0 +k1*self.deltaT/2 
            self.particles.locL = Y1[0]
            self.particles.massL = Y1[1]

            tmid = t0+self.deltaT/2 
            self.get_auxiliary(tmid)

            k2 = self.particles.dY2d_dt(Y1, tmid)
            Y2 = Y0 + k2*self.deltaT/2
            self.particles.locL = Y2[0]
            self.particles.massL = Y2[1]

            self.get_auxiliary(tmid)
            k3 = self.particles.dY2d_dt(Y2, tmid)
            Y3 = Y0 + k3*self.deltaT 
            self.particles.locL = Y3[0]
            self.particles.massL = Y3[1]
            
            self.get_auxiliary(tn)
            k4 = self.particles.dY2d_dt(Y3,   tn)
            Y4 = Y0 + k4*self.deltaT

            Yn = Y0 +1/6 *(k1 +2*k2 +2*k3 +k4) *self.deltaT
            
            #Yn = 1/6*(Y1+2*Y2+2*Y3+Y4)

        else:
            print('[core-update_particles]: the {} is not a valid integration method, please check'.format(pars.dtimesteppars))
            sys.exit(1)

        self.particles.locL = Yn[0]
        self.particles.massL = Yn[1]

        #Yt = self.particles.update(self.time,self.time+self.deltaT,self.disk,timestepn)

        if self.time==np.nan or self.deltaT==np.nan:
            print('hello')
            import pdb;pdb.set_trace()

        return Yn


    def back_up_last_data(self):
        """
        copies present state to "old" 
        """
        copy_list = ['time', 'particles', 'planetL', 'nplanet', 'icelineL', 'Minflux_step', 'dotMg', 'deltaT','centralbody']
        self.oldstate = COPY (self, copy_list)

    def remove_planet(self):
        """
        remove the planet

        seems very stupid now, let's think about this later TBD:
        """

        for num, pl in self.planetD.items():
            if pl.loc <= self.centralbody.r:
                #let's test what if we change the resonance state first
                remove_setL = []
                for ss in self.res_setL:
                    if {num}.issubset(ss):
                        plnum = list(ss-{num})[0]
                        self.planetD[plnum].resS = -1
                        remove_setL.append(ss)
                
                #then remove the planet 
                self.planetL.remove(pl)
                self.nplanet -= 1
                print('[system.remove_planet]: planet '+str(pl.number)+ ' is removed.')

                for rs in remove_setL:
                    self.res_setL.remove(rs)
        self.planetD = {p.number:p for p in self.planetL}


        if False:
            for planet in self.planetL:
                if planet.loc <= 4*cgs.RJ:
                    #remove the planet firstly
                    self.planetL.remove(planet)
                    self.nplanet -= 1

                    #remove this planet from the resonance sets and delete conresponding sets
                    uname = planet.number
                    for ss in self.res_setL:
                        if {uname}.issubset(ss):
                            # change the resS of the planet within this set to -1 
                            for p in self.planetL:
                                if {p.number}.issubset(ss):
                                    p.resS = -1

                            self.res_setL.remove(ss)


    def post_process (self):
        """
        system elements (particles, planets, icelines...) have been evolved to the next state

        Look for following changes:

        - adding and removing of particles that crossed border
            * first indicate which (daction dict)
            * remove/add them
        - resample particle distribution (Nplevel or splitmerge algorithms)
        - make new particle state vector
        - remove planets
        - update central mass

        LZX[24.08.06]: put add_planet here
        """

        self.daction = {}

        #loc = self.particles.Y2d[0] #TBR
        loc = self.particles.locL

        #particles that cross the inner disk edge
        idx, = (loc<self.rinn).nonzero()
        if len(idx)>0:
            self.daction['remove'] = idx
        
        #delmgasIn = dp.ratio*sciint.quad(dp.dot_Mg,self.time,self.time+self.deltaT)[0]
        #self.dotMd = dp.ratio*dp.dot_Mg(self.time)

        #update the gas mass flow
        self.dotMg = dp.dot_Mg(self.time)


        #get the dust mass flow
        delmdustIn = sciint.quad(dp.dot_Md, self.time, self.time+self.deltaT)[0]

        self.Minflux_step = 0# delmdustIn

        #LZX: [24.04.27] make the mtot1 the largest among these three values, so that the Nadd will 
        #      never be lager than 1
        #however: this results in jumping of mtot1, which destabelizes the scheme
        #TBD: restrict timestep such that Nadd<=1 is guaranteed

        self.particles.mtot1 = max(self.particles.mtot1, self.Minflux, self.Minflux_step)
        #but the self.mtot1 cannot be lower than (self.Minflux+self.Minflux_step)/2 
        #because this can make the adding two particles once
        self.Minflux += self.Minflux_step

        Nadd = 0#particles that enter the domain
        if pars.resampleMode=='Nplevel':
            #if the total mass has exceeded some threshold "mtot1"
            #create a new particle
            mtot1 = self.particles.mtot1
            while self.Minflux> mtot1:
                Nadd += 1
                self.Minflux -= mtot1
        elif pars.resampleMode=='splitmerge' or pars.resampleMode == 'dropmerge' or\
             pars.resampleMode in ['new_splitmerge_chris','fixed_resample','local_splitmerge'] or\
             pars.resampleMode in ['new_splitmerge_zxl'] or\
             pars.resampleMode in ['global_resample','global_resample2', 'global_resample3', 'global_resample4'] and self.rout is not None:
            mtot1 = self.particles.mtot1
            while self.Minflux> mtot1:
                Nadd += 1
                self.Minflux -= mtot1
            #dis = np.log(self.rout/self.particles.locL[-1])
            #if dis >= np.sqrt(pars.dresample['fdelS']*pars.dresample['fdelM']):
            #    Nadd += 1 
            #    #if add particles here, the total mass are always changed
            #    self.particles.mtot1 = self.Minflux 
            #    self.Minflux = 0.
        elif pars.resampleMode==None:
                self.Minflux = 0.
        else:
            print('[core.py]:No valid resampleMode, choose from: [Nplevel,splitmerge,dropmerge,global_resample,None]')
            sys.exit()

        
        if Nadd>0:
            self.daction['add'] = Nadd
        
        if Nadd>=2:
            print('there are two super particles will be add once, please check')
            import pdb;pdb.set_trace()
        
        #post_process particles
        if 'remove' in self.daction.keys():
            #remove the particles from Y2d!
            
            nrem,mrem = self.particles.remove_particles(self.daction['remove'])
            self.Moutflux += mrem

            self.add_message('remove', str(nrem)+' particles lost to inner edge crossing')


        if 'add' in self.daction.keys():
            self.particles.add_particles(self.daction['add'])


        #[24.01.04]
        #it is really difficult to stabalize particle numbers, b/c 
        #of the huge lag... I think the below algorith accomplishes smth
        #but is a bit ugly

        #[24.05.25]:splitmerge should already stabilize particle number. Do one or the other
        if pars.resampleMode=='Nplevel':
            mtot1 = self.particles.mtot1
            #to be brief..
            part = self.particles

            Np = len(part.massL)

            nch = 4
            if Np==part.ninit: part.Nplevel = part.ninit #reset
            
            #particle level is decreasing...
            if Np<part.Nplevel -nch and Np<part.ninit:
                part.Nplevel -= nch
                eps = abs(part.ninit -part.Nplevel) /part.ninit
                mtot1 *= 1 -eps

            #particle level is decreasing, but above ninit: modest decrease mtot1
            elif Np<part.Nplevel -nch and Np>part.ninit:
                part.Nplevel -= nch
                eps = nch /part.ninit
                mtot1 *= 1 -eps

            #particle level is increasing, but below ninit: modest increase mtot1
            elif Np>part.Nplevel +nch and Np<part.ninit:
                part.Nplevel += nch
                eps = nch /part.ninit
                mtot1 *= 1 +eps

            #particle level is increasing...
            elif Np>part.Nplevel +nch and Np>part.ninit:
                part.Nplevel += nch
                eps = abs(part.ninit -part.Nplevel) /part.ninit
                mtot1 *= 1 +eps

            #update the total mass 
            self.particles.mtot1 =mtot1

        #TBD: add particles when outmost particles is far from the rout

        #get the Y2d needed to be used next step
        #TBR
        #self.particles.make_Y2d()

        # get the property list to remove/add particles and select_single
        self.particles.propL = [attr for attr in dir(self.particles) if not attr.startswith('__') and isinstance(getattr(self.particles, attr), list) or isinstance(getattr(self.particles, attr), np.ndarray)]   
        self.particles.propSol = [attr for attr in dir(self.particles) if not attr.startswith('__') and isinstance(getattr(self.particles, attr), float)]
        
        # remove the planet if the planet location is too small, not sure if it should be here
        self.remove_planet()
        self.centralbody.update(self.time)

        #update the rinn and rout in case:
        ## 1.rinn may change with radius of central mass. 
        ## 2.rout may change b/c of viscous spreading.
        self.get_rinn()
        self.get_rout()
        


    def new_timestep (self, tEnd, deltaTfraction=0.2, afterjump = False, jumpfracD={},**kwargs):
        """
        - determine w/r planets end up in resonance
        - chooses a timestep

        As such, we look at the following processes
        - the change in the particles state vector Y/Y_t (tpart)
        - the relative motions among the particles (tcol -- this involves deltaTfraction)
        - the change in the planet's mass/location/composition (PxxxTscale)
        - the growth rate of the central mass object (McTscale]
        - The changing timescale of mass inflow rate (mdotgTscale)
        - the timescale on which resonances are approach (pratTscale)

        """
        #[24.07.26]TBD: this function is too long. It may become incomprehensible what's going on. 
        #               it needs to better commented and broken up where possible

        #organize the procedure a bit (for quasi-steady evolution... later!)
        mintimeL = []
        
        Y2d = self.particles.make_Y2d()
        Y2dp = self.particles.dY2d_dt(Y2d,self.time)

        #timescale for the particles
        tpart = np.abs(Y2d/Y2dp)

        #get the collision timescale of particles 
        #[2024.08.15]LZX: If the 0.5 here is absent, the particles crossing will happen
        tcol = np.append(np.inf, abs(np.diff(self.particles.locL)/np.diff(Y2dp[0])) *pars.dtimesteppars['coltimefrac'])
        #tcol = np.append(np.inf, abs(np.diff(self.particles.locL)/np.diff(Y2dp[0])))

        tpart = np.concatenate((tpart, tcol[np.newaxis,:]))

        #tmin = min([tpart[0].min(),tpart[1,2:].min(),tpart[2].min()])
        mintimeL.append({'name':'particles', 'tmin': np.nanmin(tpart), 
                                'imin':np.unravel_index(tpart.argmin(),tpart.shape)})


        #if self.ntime ==760:
        #    import pdb;pdb.set_trace()
        #after the jump, we only have particles drift timescale
        if afterjump:
            self.minTimes = Mintimes(mintimeL,jumpfracD)
            #self.mintimeL = mintimeL
            deltaT = deltaTfraction*tpart.min() #min([ob['tmin'] for ob in mintimeL])

            #make sure timestep doesn't jump as well...
            self.deltaT = min(deltaT, self.oldstate.deltaT*1.01)

            return


        McTscale = np.inf
        if self.time > 0:
            with np.errstate(divide='ignore', invalid='ignore'):
                McTscale = self.centralbody.m/ np.abs(self.centralbody.m - self.oldstate.centralbody.m) *self.deltaT
            mintimeL.append({'name': 'CentralMassGrowth', 'tmin': McTscale})
            # import pdb; pdb.set_trace()

            #rdum = (1-0.1*self.particles.delta)*self.rinn
            #tinn = (self.particles.locL[0] -rdum) /(-Y2dp[0,0])
            #mintimeL.append({'name': 'innercrossTime', 'tmin': tinn})


        #calculate mass flow change Timescale
        if self.oldstate is not None:   
            #Mass influx timescale
            #[24.05.12]cwo: I removed the "1e-100" from the denominator b/c when dotMg==0 this timescale should become infinite
            with np.errstate(divide='ignore', invalid='ignore'):
                mdotgTscale = (1e-100 + np.float64(self.dotMg)) / (abs(self.oldstate.dotMg - self.dotMg)) *self.deltaT
            mintimeL.append({'name': 'Mass_Influx', 'tmin': mdotgTscale})
            #timescale for the planets 
            # (including migration and mass growth)
        

        if pars.doPlanets and self.nplanet > 0:
            PmassTscale = np.inf*np.ones(self.nplanet) #mass growth T.
            PlocaTscale = np.inf*np.ones(self.nplanet) #migration T. 
            PcompTscale = np.inf*np.ones(self.nplanet) #composition T. (later)
            pratTscale = np.inf*np.ones(self.nplanet)


       
        #loop over all planets to monitor their changes in ...
        for i in range(self.nplanet):
            planet = self.planetL[i]

            #match the planet with its oldstate by their name (unique number)
            uname = planet.number #its unique name

            uoldL = [planet.number for planet in self.oldstate.planetL]
            try:
                iold = uoldL.index(uname)
            except:#it did not exist yet
                iold = -1


            #planet migration T.
            #oldstate needs to have the same shape (otherwise not the same planet)
            PlocaTscale[i] = np.float64(planet.loc)/abs(planet.dlocdt)                

            #TBD (maybe): merge the res.trapping stuff to post_process (?)
            if pars.doResonance:

                if i>0 and iold>0:
                    #print(self.planetL[1].inxt,(self.planetL[1].loc/self.planetL[0].loc)**(3/2))
                    jres = planet.inxt +1#for 3:2 reson, inxt=1, and j=2 (j+1:j)
                    prat = (planet.loc/self.planetL[i-1].loc)**(3/2)
                    pratold = (self.oldstate.planetL[iold].loc/self.oldstate.planetL[iold-1].loc)**(3/2) #please check...

                    #resonace offset (positive quantity)
                    pdel = prat -(jres+1)/jres
                    pdelold = pratold -(jres+1)/jres
                    if pdel == 0.:
                        import pdb;pdb.set_trace()
                    
                    #planets are in resonance
                    if planet.resS=='R':
                        pratTscale[i] = np.inf

                    # change the resonance state here, id the pdel is small enough, then teh resS = 'R'
                    # for now we don't consider the condition that the planet jump over the resonance.
                    elif pdel <= 0.0: #a little arbitrary now
                        
                        ta = PlocaTscale[i]
                        qinn = self.planetD[uname- 1].mass/ self.centralbody.m


                        disk = self.get_disk(loc = planet.loc)
                        Hg = disk.Hg
                        #Hg = f_Hg (planet.loc)
                        #plt.plot(self.particles.locL, self.particles.Hg)
                        #plt.scatter(planet.loc, f_Hg(planet.loc), )
                        haspect = Hg/planet.loc
                        tPer = 2*np.pi/physics.Omega_K(planet.loc, self.centralbody.m)
                        getTrapped = physics.crossedResonance (ta, jres, qinn, haspect, tPer)
                        
                        if getTrapped:
                            planet.resS = 'R'
                            # make sure every planet in the rasonance chain will be taken into consider when getting the migration rate
                            self.planetL[uname-1].resS = 'R'

                            print('[core.new_timestep]: planet {} and planet {} has been trapped in {}:{} resonance'.format(uname, uname-1, planet.inxt+2, planet.inxt+1))
                            
                            # make the resonance pairs into set list
                            res_set = {self.planetL[i-1].number, uname}
                            if res_set not in self.res_setL:
                                self.res_setL.append(res_set)
                        #trapping fails
                        else:
                            print ('[core.new_timestep]: failed trapping planet {} and planet {} into {}:{} resonance'.format(uname, uname-1, planet.inxt+2, planet.inxt+1))
                            planet.inxt += 1
                            # and need the pratTscale immediately to avoid to be jumped over
                            jres = planet.inxt +1
                            pdel = prat - (jres+1)/jres
                            pdelold = pratold - (jres+1)/jres
                            pratTscale[i] = np.float64(pdel) /(1e-100+pdelold-pdel) *self.deltaT 
                    

                    #calculate how fast the planets approach resonance 
                    #the timescale on which planets approach resonance
                    else:#pdel>0
                        # +1e-2 to prevent a too small value comes out
                        pratTscale[i] = np.float64(pdel + 1e-2) /(1e-100 +pdelold-pdel) *self.deltaT
                        #get the proper resonance inxt
                        planet.inxt = (self.dres['prat']<prat).argmax()



            #fit the planet growth by pebble accretion

            #store mass data first
            if iold>=0 and self.oldstate.planetL[i].mass != planet.mass:
                # self.masstime=
                planet.planetMassData.append([self.time , planet.mass])

            #[24.08.04] It seems more like smth for post_process 
            # if the planet cross the inner edge, then the accretion is False
            if planet.loc< self.rinn*(1-1e-4):
                planet.accretion =False

            #then try to fit the mass to a curve
            Npts = len(planet.planetMassData)
            Nfit = 10
            #consider the data is not large enough to make the fit

            #TBD:if particles evaporate before reaching the planet
            #if planet.loc < max(location_most_inner_iceline, cavity_radius):
            #    planet.dmdt = 0.0
            #    #max_jumpT = np.inf
            # if the planet mmigrates into cavity, ths max_jumpT should be infinity
            #if planet.loc < self.rinn:
            #    planet.dmdt = 0
            #    planet.relp_mass = np.inf
            #    planet.max_jumpT = np.inf
            planet.get_maxJumpT(Npts, Nfit)
            PmassTscale[i] = planet.massTscale

            #get the planet composition change time scale
            if False:
                if self.oldstate.planetL[i].fcomp[0] != planet.fcomp[0]:
                    planet.compData.append([self.time, planet.fcomp[0]])
                
                if len(planet.compData) >10:
                    timedots, compdots = np.array(planet.compData).T
                    timedots = np.log(timedots)
                    def comp_fit(t,a,b):
                        c = a*t+b
                        return c
                    
                    relp = np.inf
                    while Nfit<=Npts:
                        popt, pcov = curve_fit(mass_fit, timedots[-Nfit:], compdots[-Nfit:])
                        #line = '{:7.4f} {:10.3e} {:10.3e} {:6d}'.format(popt[0], np.sqrt(pcov[0,0]), np.sqrt(pcov[0,0])/popt[0], Nfit)
                        #print(line)
                        if np.sqrt(pcov[0,0])/np.abs(popt[0])<relp:
                            pidx = popt[0]
                            psig = np.sqrt(pcov[0,0])
                            relp = psig/np.abs(pidx)
                        else:
                            break
                        Nfit = int(Nfit*1.2) +1
                    planet.relp_comp = relp
                    
                    PcompTscale[i] = (np.exp(timedots[-1])-planet.starttime)* compdots[-1]/abs(popt[0]) 
                    planet.comp_tc = PcompTscale[i]

                else:
                    planet.relp_comp = np.nan
        
        if pars.doPlanets and self.nplanet >0:                        
            mintimeL.append({'name': 'planetsMigration', 'tmin': min(PlocaTscale)})
            mintimeL.append({'name': 'planetsGrowth', 'tmin': min(PmassTscale)})
            mintimeL.append({'name': 'planetsComp', 'tmin': min(PcompTscale)})
            
            if pars.doResonance:
                for key,value in self.milestones.items():
                    if value == 'resonance':
                        del self.milestones[key]
                        break
                 
                self.milestones[self.time+ 1e-3 +min(pratTscale)] = 'resonance'
                #TBD: find out WHY do we see 0 // perhaps put pratTscale->infinity when in resonance (delta=0)?
                mintimeL.append({'name': 'PlanetsRes', 'tmin': np.min(pratTscale[pratTscale >= 0.0])})
                
        #timescale for the icelines
        if pars.doIcelines and self.oldstate is not None:
            IlocaTscale=np.inf*np.ones_like(self.icelineL)
            for i,iceline in enumerate(self.icelineL):
                if not np.isnan(iceline.loc):
                    #[24.02.20]cwo:added a small number to the denominator
                    tscale = np.float64(iceline.loc)/(abs(self.oldstate.icelineL[i].loc-iceline.loc))*self.deltaT
                    iceline.loc_tc = tscale
                    IlocaTscale[i] = tscale
        
            mintimeL.append({'name': 'icelineloca', 'tmin': min(IlocaTscale)})


        ## We are (finally) ready to determine the new timestep (deltaT)
    
        # put mintimeL into system object for now to check
        self.minTimes = Mintimes(mintimeL, jumpfracD)

        #determine next timestep
        #[25.01.20]: maybe deltaTfraction not applied to ALL, but only to
        deltaT = deltaTfraction *min(self.minTimes.tminarr)

        #limit increase of deltaT to (some number)
        if self.ntime>0:
            deltaT = min(deltaT, self.oldstate.deltaT*1.05)

            
        if self.time+deltaT>tEnd:
            deltaT = tEnd - self.time

        self.deltaT = deltaT

        if self.time != tEnd:
            if self.deltaT/self.time <1e-15:
                print('[core.new_timestep]: warning: may lose precision b/c of small time step')
                import pdb;pdb.set_trace()
            if self.deltaT<=0:
                print('[core.new_timestep]warning:deltaT<=0')



    def query_system_jump (self):
        """
        investigates w/r we can jump and by how much
        ...

        returns True/False, {jump properties}
        """

        #[24.05.12]cwo:jump should be an option (I thought it was so already...)
        if pars.doJump==False:
            return {'jumpT':False}

        Tscale_ratio = []
        for t in self.minTimes.tminarr[1:]:
            Tscale_ratio.append( t/self.minTimes.particles)
        #print(self.minTimes.tminarr) 
        if len(self.minTimes.tminarr[1:])==0:
            self.doJump = False
            return {}

        #jump time is given by the min evol. timescales, excluding those of the particles
        max_tevol = self.minTimes.max_tevol
        #max_tevol = jumpfracD*min(self.minTimes.tminarr[1:])
        
        #(cannot jump over "important events" -> Milestones) 
        #adjust jumpT according to milestones
        timepoints = np.sort( np.array(list(self.milestones.keys())) )

        ii = timepoints >self.time
        if sum(ii)==0:
            max_tms = np.inf
        else:
            #if timepoints[ii][0] == pars.tmax and self.time>11.8e6*cgs.yr:
            #    import pdb;pdb.set_trace()
            max_tms = min(timepoints[ii]-self.time)
        
         
        #reached_ms = np.where( (timepoints - self.time - jumpT<0) & (timepoints - self.time>0) )[0] 
        #we cannot jump over milestones
        #if len(reached_ms) != 0:
        #   jumpT = timepoints[reached_ms[0]] - self.time
        
        #jumpT cannot exceed max_jumpT (see above)
        #NOTE: planet may not exist... hence the for/if construct
        max_tpl = np.inf
        for planet in self.planetL:
            if self.time>planet.starttime:
                max_tpl = min(max_tpl, planet.max_jumpT)


        #the jump time is the minimum
        #[24.01.21]I wrote this a bit more clearly...
        tjumpkeys = ['system-evol', 'milestones', 'planet-growth-fit']
        tjumparr = np.array([max_tevol, max_tms, max_tpl])
        jumpT = min(tjumparr)

        #evol.timescales >> drift timescales && ensure significant jump
        con0 = self.ntime >self.njumptime +100 #at least wait 100 steps
        #con1 = (min(Tscale_ratio) > 1e2) and (jumpT/self.deltaT>100) #should have more judgements
        con1 = (min(Tscale_ratio) > 1e2) #should have more judgements

        #perhaps
        if len(self.timeL)>100:
            con2 = (jumpT>self.time-self.timeL[-100])
        else:
            con2 = (jumpT>(self.time-self.timeL[0])*100/len(self.timeL))

        
        #initial relaxation time
        #[24.06.10]LZX: consider the iceline effect, the outflux should be larger 
        #than the total silicate mass
        con3 = self.Moutflux>self.particles.fcompini[0]*self.minitDisk
        print(self.Moutflux/self.minitDisk/0.5)


        #if self.ntime%1000==0: 
        #    print(con0, con1, min(Tscale_ratio), jumpT, np.argmin(tjumparr))

        if False:
            # the fit goodness of planets grow >0.95, if the r2 is np.inf, the planets isn't added, if the planet is already added but r2 is still np.nan, then the mass data points are not enough 
            relp_massL = np.array([p.relp_mass for p in self.planetL if p.starttime < self.time])
            if len(relp_massL) != 0:
                con2 = relp_massL.min() < 2e-3
            else:
                con2 = False
            #if self.time > 50*cgs.yr:
            #    import pdb; pdb.set_trace()
            # consider the fit goodness of the planet's composition change, similar with the mass
            relp_compL = np.array([abs(p.relp_comp) for p in self.planetL if p.starttime < self.time])
            if len(relp_compL) != 0:
                con3 = relp_compL.min() < 2e-3
            else:
                con3 = False

            self.doJump = con1 & con2 &con3
        else:
            self.doJump = con0 & con1 & con2 & con3

        
        #print([con0,con1,self.mintimeL[1:], self.time/cgs.yr]) 
        djump = {'jumpT':jumpT, 'tjumpkeys':tjumpkeys, 'tjumparr':tjumparr}

        self.con_Jump = [con0, con1, con2, con3]

        return djump


    def reset_after_jump(self):
        #reset the planet mass data
        for planet in self.planetL:
            planet.planetMassData = []
            planet.relp_mass = np.inf #Chris: never assign "nan" please
            planet.max_jumpT = 0.0
            #planet.compData = []
            #planet.relp_comp = np.nan
        
    def system_jump(self, djump):
        """
        execute the system jump
        """
        self.jumpT = djump['jumpT']
        
        # parameters needs to be updated:
        # planets: location and mass and composition(this maybe very complex)
        # icelines :location
        
        if self.nplanet>0:
            for planet in self.planetL:
                    
                planet.loc += planet.dlocdt *self.jumpT
                jumpmass = planet.dmdt* self.jumpT

               #TBD: generalize this. Perhaps best way is to make planet.dmdt a vector
                #       planet.dmdt = [dmdt comp 1, dmdt comp 2, ...]
                #paridx = np.argmin(abs(self.particles.locL - planet.loc))
                #if planet.loc < self.icelineL[0].loc*(0.5) and self.particles.fcomp[paridx][0] ==0.5:
                #    import pdb;pdb.set_trace()

                #NOTE:change the composition jump now more resonable but not general
                if pars.doIcelines:
                    if planet.loc <self.icelineL[0].loc:
                        fcomp=np.array([1,0])
                    else:
                        fcomp = np.array([0.5,0.5])
                else: 
                    fcomp = np.array([1,0])

                planet.fcomp = (fcomp*jumpmass +planet.mass*planet.fcomp)/(planet.mass+jumpmass)
                planet.mass += jumpmass
            #if self.planetL[0].loc <self.rinn:
            #    print('[core.system.jump]: the first planet migrates across the inner edge')

        if pars.doIcelines:
            for i, iceline in enumerate(self.icelineL):
                #to dicide whether the iceline moves inner or outer
                sign = (iceline.loc-self.oldstate.icelineL[i].loc)/abs(iceline.loc-self.oldstate.icelineL[i].loc)
                iceline.loc += sign*iceline.loc/self.minTimes.icelineloca *self.jumpT
        
        self.njump +=1
        self.njumptime = self.ntime
        self.jumptime = self.time

        #im = djump['tjumparr'].argmin()
        #maybe interesting to store and plot which factor limits the jumpT
        #self.jump_limitation = djump["tjumpkeys"][im]

        #nameL = [d['name'] for d in self.mintimeL]
        #tminarr = np.array([d['tmin'] for d in self.mintimeL])
        #imin = 1 +tminarr[1:].argmin() #minimum evolution

        #print(f'[core.system.jump]:at {self.time/cgs.yr:5.2f} yr jumped by {self.jumpT/cgs.yr:5.2f} yr')
        #print(f'[core.system_jump]:jump time limited by: {self.jump_limitation}')
        #print(f'[core.system_jump]:min. evolution time ({nameL[imin]}) {tminarr[imin]/cgs.yr:9.2e} yr')
        
           #import pdb;pdb.set_trace()

        #"erase" previous planet.crossL OR record the jump time to planet.
        #such that new fit for dm/dt starts w/ N=0 particles


def get_cross_idx (loc, locL, locLo, daction, locnew = None):
    """
    this is needed for the 'advance' things to avoid some error caused by 
    removing and adding particles
    """
    lrm = 0 #[25.01.23]cwo:??
    lad = 0
    #daction={}
    #make up the locL:
    #if 'remove' in daction.keys():
    #    for pos in sorted(daction['remove']):
    #        locL=np.insert(locL, pos, 0.)
    #if 'add' in daction.keys():
    #    lad = daction['add']
    #    locLo=np.append(locLo, [np.inf*lad])


    if locnew == None:
        idx,=np.nonzero((loc< locLo) & (loc>locL))
    else:
        idx,=np.nonzero((loc< locLo) & (locnew>locL))

    #[25.01.23]cwo:I dont understand this...
    idxD = {'idx_for_new': idx-lrm, 'idx_for_old': idx}

    #if loc<5.89*cgs.RJ:
    #    print(loc/cgs.RJ)
    #    print(np.append([0.]*lrm , locL)[0:4]/cgs.RJ)
    #    print(locLo[0:4]/cgs.RJ)
    #    import pdb;pdb.set_trace()

    return idxD
    

def advance_iceline (system):
    """
    for now particles directly lose the mass of water without any other effect
    """

    for k,iceline in enumerate(system.icelineL):
        ic = pars.composL.index(iceline.species) #refers to species index

        #[25.01.23]cwo: why this?
        if False:
            for i in range(system.particles.num):
                fice = system.particles.fcomp[i,ic]  #mass fraction in ice
                fremain = (1-fice)          #remain fraction

                if fremain < 1e-15:
                    fremain=0 #loss of numbers (!!)
                if fice!=0 and system.particles.locL[i]<iceline.loc:
                    system.particles.mtotL[i] *= fremain    #reduce masses accordingly
                    system.particles.massL[i] *= fremain
                    system.particles.fcomp[i,ic] = 0.      #gone is the ice!
                    #renormalize
                    system.particles.fcomp[i,:] = (system.particles.fcomp[i,:].T /(system.particles.fcomp[i,:].sum()+1e-100)).T
                    
                    #import pdb;pdb.set_trace()

            #renew the iceline location
            loc_pv = system.oldstate.icelineL[k].loc
            iceline.get_icelines_location(system.gas,system.time,bounds= (system.rinn, system.rout), guess=loc_pv)


        sploc = system.particles.locL
        sploc_old = system.oldstate.particles.locL
        idxD = get_cross_idx(iceline.loc,sploc,sploc_old, system.daction)
        idx = idxD['idx_for_new']
        if len(idx)!=0:     
           
            for ix in idx:
                fice = system.particles.fcomp[ix,ic]  #mass fraction in ice
                fremain = (1-fice)          #remain fraction

                if fremain < 1e-15:
                    fremain=0 #loss of numbers (!!)
                system.particles.mtotL[ix] *= fremain    #reduce masses accordingly
                system.particles.massL[ix] *= fremain
                system.particles.fcomp[ix,ic] = 0.      #gone is the ice!

               #renormalize
                system.particles.fcomp[ix,:] = (system.particles.fcomp[ix,:].T /(system.particles.fcomp[ix,:].sum()+1e-100)).T
        
        loc_pv = system.oldstate.icelineL[k].loc
        iceline.get_icelines_location(system.gas,system.time,guess=loc_pv)



def advance_planets (system):
    """
    [23.12.06]copied/edited from NewLagrange
    """
    res_chainL = ff.get_res_chain(system.res_setL)

    sploc = system.particles.locL
    sploc_old = system.oldstate.particles.locL


    for planet in system.planetL:

        #planet exists only after planet.time
        if planet.starttime<system.time:

            #particles that cross are those that
            idxD = get_cross_idx(planet.loc,sploc,sploc_old,system.daction)
            #idx, = np.nonzero( (planet.loc<sploc_old) & (planet.loc>sploc) )


            iterate = True
            niter = 0
            while iterate:


                crossL = []
                for ip in idxD['idx_for_old']:

                    spi = system.oldstate.particles.select_single(ip)
                    crossL.append(spi)

                #crossL=np.array(crossL)

                
                # if len(idx)!=0:
                #                     #this user-defined function should detail how the
                #planet properties (location,mass,composition) change
                #with time and how the crossed particles are affected... 


                #this is about planet migration...
                #TBD later...
                if False:
                    ## CWO: switch to the user-defined approach
                    loc_t, mass_t, fcomp_t = \
                    userfun.XY_planet (sim.time, planet.loc, planet.mass, planet.fcomp, 
                            crossL)
                else: 
                    
                    try:
                        resS = planet.resS
                    except:
                        resS = None

                    if resS == 'R':
                        chain = ff.locate_chain(res_chainL, planet.number) 
                        invtmigL = [] #inverse migration time
                        weightL = []
                        for num in chain:
                            p = system.planetD[num]
                            dum_t = userfun.planet_migration(system.gas,p.loc,p.mass, system.time, system.rhoPlanet)
                            invtmigL.append(dum_t/p.loc)
                            weightL.append(p.mass*p.loc**0.5)

                        #joint migration timescale
                        #TBD: this should not be like this, should get from the initial 
                        #magnetic field, make this an option
                        if 0.0 in invtmigL:
                            invmigtime =0.
                        else:
                            invmigtime = np.sum(np.array(weightL)*np.array(invtmigL)) /np.sum(np.array(weightL))
                        loc_t = planet.loc *invmigtime
                    else:
                        loc_t = userfun.planet_migration(system.gas,planet.loc,planet.mass, system.time, system.rhoPlanet)
                        
                    mass_t = 0.0    #gas accretion of planet, TBD:-later
                    fcomp_t = 0.0   #how its composition changes

                
                # set a milestone: when planet will reach to the inner edge
                msg = 'planet-reach-rinn-'+str(planet.number)
                #update the time corr. to the message
                #find the key corresponding to the value
                for key,val in system.milestones.items():
                    if val==msg:
                        del system.milestones[key] ## doesnt work?
                        break

                #add/update milestone
                #LZX: [24.07.30] add a pre-fractor here to prevent the overshooting 
                key = 0.9*np.float64((planet.loc-system.rinn))/abs(loc_t)
                system.milestones[key + system.time] = msg

                #update planet properties from rates supplied by user
                planet_loc_nw = planet.loc + loc_t *system.deltaT

                #particles that cross are those that
                idxND = get_cross_idx(planet.loc, sploc, sploc_old, system.daction, planet_loc_nw)
                #idxN, = np.nonzero( (planet.loc<sploc_old) & (sploc<planet_loc_nw) )
                if set(idxND['idx_for_new'])!=set(idxD['idx_for_new']):
                    idxD['idx_for_new'] = idxND['idx_for_new']
                    niter += 1
                else:
                    iterate = False


            #update planet properties from rates supplied by user
            planet.loc += loc_t *system.deltaT
            planet.dlocdt = loc_t
            
            ## TBD 
            # update some planet properties (like pdel) here?
            # seems most natural point (after planets have advanced)
            
            planet.time = system.time  ## add this time to planet

            #planet.mass += mass_t *system.deltaT
            #planet.fcomp += fcomp_t *system.deltaT

            #update s-particle properties from sp-crossings
            #assumes all particles are accreted (TBC!!)
            spN = system.particles

            for k, ip in enumerate(idxND['idx_for_new']):

                #calculate critical mass to verify if the pebble accretion can occur
                 
                #[24.01.05]
                ## CWO: I dont think it's a good idea to forward the entire system
                #       class to such simple functions..
                #
                #       (for the moment I ignore it... we can discuss)
                spk = crossL[k]
                Mc = userfun.M_critical(spk.eta, spk.St, spk.mcp)

                #mass*composition that is acccreted
                delmcomp = np.zeros((spN.fcomp.shape[1]))
                #TBD-later: I dont like the need for an M_critical userfun..
                #   instead, we can incorporate this into userfun.epsilon_PA
                #
                #On the other hand, particles may get stuck in pressure bump   
                #and fail to accrete and drift. This behavior would be nice to capture
                #but this requires more thinking..
                #
                #TBD-later: in reality particles may be "stuck" in pressure bump
                #           incorporate in planet object?
                #NOTE: this PIM is arbitrary now
                if Mc<planet.mass and planet.mass<userfun.PIM():                    
                    epsilon = min(0.99,userfun.epsilon_PA(planet.loc,planet.mass,spk))

                    planet.acc_rate =  epsilon
                    #accreted mass by composition
                    delmcomp += epsilon*crossL[k].fcomp *crossL[k].mtotL

                    #spN -> system.particles.Y2d...
                    spN.mtotL[ip] -= delmcomp.sum() #decrease mass sp
                    masscomp = planet.fcomp*planet.mass +delmcomp
                    planet.mass += delmcomp.sum()
                    planet.fcomp = masscomp /planet.mass
                    if spN.mtotL[ip]<=0. or epsilon>1:
                        import pdb;pdb.set_trace()


                elif planet.mass>userfun.PIM() and planet.accretion:
                    print('planet {:} has reach the PIM at {:10.1f}'.format(planet.number,system.time/cgs.yr ))
                    planet.accretion = False

                    #delm = epsilon*crossL[k].mtotL
                    #for i in range(len(crossL[k].fcomp)):
                    #    delmass[i] = epsilon* crossL[k].fcomp[i]*crossL[k].mtotL
                    
                    #planet.mass*planet.fcomp
                    
                    #planet.fcomp = [ (delmass[i]+planet.mass*planet.fcomp[i]) / (planet.mass+delmass.sum()) for i in range(len(delmass))]
                    #planet.mass += delmass.sum() #increase mass (pebble accretion)
                #else:
                    #print('[core]: pebble accretion can not happen')
                    #import pdb; pdb.set_trace()


                # planet.fcomp += 0.  #TBD !!
                

class SingleSP(object):
    """
    Used to get single Superparticle
    """
    def __init__ (self,**kwargs):
        for key,val in kwargs.items():
            setattr(self,key,val)


class Superparticles (object):

    def __init__(self, rinn, rout, dcomposL, gas, nini=40, Rdi=0.1, 
            user_init_radius = None,
            initrule='equalmass'):
        """
        systems initial properties

        initrule:[equalmass,equallogspace]
                :how initial parameters are distributed

        gas     :gas object (needed to get gas surface density)   

        nini    : initial number of the particles
        rinn    : inner edge of the disk
        rout    : outer edge of the disk
        dcomposL: the composition list

        """
        self.Rdi = Rdi
        
        self.rhocompos=[]
        for compos in dcomposL:
            if compos['name']!= 'gas':
                self.rhocompos.append(compos['rhoint'])
        
        self.ninit =nini
        self.Nplevel = nini

        #[23.12.30]this was commented out; now uncommented

        self.rinn=rinn
        self.rout=rout
        self.stokesOld = None
        
        #[23.12.30]:copied from /NewLagrance
        def construct_farr (dcomposL):
            fnL = []; miL = []
            for dcompos in dcomposL: 
                fnL.append(dcompos['Z_init'])
                miL.append(dcompos['mask_icl'])

            def f_arr (rad):
                farr = [fn(rad) for fn in fnL]
                mirr = [mi(rad) for mi in miL]
                return np.array(farr) * np.maximum(0, mirr)
            
            return f_arr

        #function gives the initial abundance of species in disk
        f_arr = construct_farr (dcomposL)

        def f_sample (r):
            #samples the initial solid mass distribution
            Zcom = f_arr(r)

            #get the initial surface density
            sigini, *dum = gas.get_key_disk_properties(r,0.0)

            #TBD: think about one way to remove the dp thing here
            #[24.04.21] CWO: This should be done arleady dp.ratio = Zcom.sum() !!
            #     (so I removed)       
            #sigini_d = dp.ratio*sigini
            return 2*np.pi*r*sigini *Zcom.sum()

        #this is the desired spacing among the particles in log-space
        self.delta = np.log(rout/rinn) /nini

        #grid used in calculating particle surface density
        #it should (?) be courser than aimed particle number
        self.pgrid = 10**np.linspace(np.log10(rinn), np.log10(rout), nini//5)
        
        #divide domain into pieces, as determined by iceline
        locspecL = []
        for comp in dcomposL:
            if comp['iceline'] == True:
                #locspec = comp['iceline_init'] ## perhaps iceline_init should be removed?
                locspec = comp['rice_init']
                if rinn<locspec<rout:
                    locspecL.append(locspec)
        locspecL.sort()

        print('[core.Superparticles.init]:initialization superparticles under rule:', initrule)
        radL = []
        rmidL = []
        msupL = []
        ncomp = len(pars.composL) #number of refractory +volatile species
        fcompL = []

        if pars.resampleMode=='fixed_resample' and initrule=='equallogspace':
            
            radL, locmid = self.loc_init (specL=locspecL)

            r0 = locmid[0]
            for k,r1 in enumerate(locmid[1:]): #the midpoints
                msup, err = sciint.quad(f_sample, r0, r1, limit =100)
                msupL.append(msup) #total mass
                fcdum = f_arr(np.sqrt(r0*r1))
                fc = fcdum[:ncomp] /fcdum.sum()
                fcompL.append(fc) #composition
                r0 = r1

            msup = np.array(msupL)

        elif initrule=='equallogspace':
            #put sp at equal distances in log space

            locspecL.insert(0, rinn)
            locspecL.append(rout)
            nspecial = len(locspecL) #includes the boundaries

            #piecewise..
            for iloc in range(nspecial-1):
                r0 = locspecL[iloc]
                xdum = np.log(locspecL[iloc+1]/r0)
                Nadd = round(xdum /self.delta)
                #[25.01.01]cwo: put particles at half-distance near boundaries
                rmid = r0 *np.exp(np.linspace(0,1,Nadd+1)*xdum)
                radL.append(np.sqrt(rmid[1:]*rmid[:-1]))

                for k,r1 in enumerate(rmid[1:]): #the midpoints
                    msup, err = sciint.quad(f_sample, r0, r1, limit =100)
                   #print('{:11.3e} {:11.3e} {:11.3e} {:11.3e}'.format(f_sample(r0), f_sample(r1), r1, msup))
                    msupL.append(msup) #total mass
                    fcdum = f_arr(np.sqrt(r0*r1))
                    fc = fcdum[:ncomp] /fcdum.sum()
                    fcompL.append(fc) #composition
                    r0 = r1

            radL = np.concatenate(radL) #single array
            msup = np.array(msupL)

        elif initrule=='equalmass':
            print('TBD: equalmass spacing')
            sys.exit()
            #puts superparticles at location such that they have
            #equald mass

            Mtot, err = sciint.quad(f_sample, rinn, rout)
            t0 = time.time()
            print('[core.Superparticles.init]:calling rout.sample_equal... this may take a while')
            radL = ff.sample_equal (f_sample, rinn, rout, nini)
            radL = np.array(radL)
            radL = np.sqrt(radL[1:]*radL[:-1])
            t1 = time.time()
            print('[core.Superparticles.init]:sampling done ({:4.1f} sec)'.format(t1-t0))

            #the mass of the super-particles are equally spread through
            msup = np.ones_like(radL) *Mtot/nini

        self.locL = np.array(radL)
        self.mtotL = np.array(msup)
        self.fcomp = np.array(fcompL)
        self.num = len(self.locL) #moved this below...

        #TBR
        if False:
            #[23.12.30]NEW:add composition data (fcomp)
            #[23.12.30]this looks a bit ugly...
            self.fcomp = np.empty((nini,len(pars.composL)))
            for k,rad in enumerate(radL):
                Zcomp = []
                for ic,scomp in enumerate(pars.composL):
                    Zcomp.append(dcomposL[ic]['Z_init'](rad)*max(0,dcomposL[ic]['mask_icl'](rad)))
                Zcomp = np.array(Zcomp)

                #Zcomp = np.array([dcompos['Z_init'](rad)*max(0,dcompos['mask_icl'](rad)) for dcompos in dcomposL])
                self.fcomp[k,:] = Zcomp/sum(Zcomp)


        #this creates self.rhoint
        self.get_rhoint()

        #if Stokes number is fixed, calculate initial radii
        if pars.dragmodel == 'fixed_St': 
            out = gas.get_key_disk_properties(self.locL, 0.0)
            Rdi = pars.fixed_St*2*out[0]/np.pi/self.rhoint

        if pars.fraginit:
            out = gas.get_key_disk_properties(self.locL, 0.0)

            #CWO? what is compmask?? Do we need it?
            compmask=np.array([])
            for i in range(nini):
                compmask = np.append(compmask,1-sum(self.fcompini[np.argwhere(self.fcomp[i]==0)]))
        

            #if we want the particles to initially reach the fragmentation velocity, then we solve for the initial Rdi with fsolve 
            from scipy.optimize import fsolve 
            def func(Rd,disk,rhoint):
                St,vr = ff.Stokes_number(disk, Rd, rhoint, Sto =0.03)
                return vr 

            def func2(Rdi, vfrag,disk,rhoint):
                rere = func(Rdi,disk,rhoint)+vfrag 
                return rere

            def get_temporary_disk(loc):
                out = gas.get_key_disk_properties(loc, 0.0)
                disk = physics.DISK(*out, loc, 0.0, dp.Mcp_t(0.0)) 
                disk.add_auxiliary()
                userparL = disk.add_uservar (dp.user_add_var())    #variables
                userfuncL = disk.add_userfun (dp.user_add_fun())    #functions only
                userevalL = disk.add_user_eval (dp.user_add_eval()) #evaluations
                return disk 

            idx = np.argwhere(compmask>0.5).min()

            Rdi = np.zeros(nini)
            for i in range(idx):
                disk = get_temporary_disk(self.locL[i])

                initguess = 0.01
                Rdi[i] = fsolve(func2, initguess, args=(pars.vc['silicates']*2, disk, self.rhoint[i]))[0]

            for i in range(idx,nini):
                disk = get_temporary_disk(self.locL[i])
                initguess = 0.1 
                Rdi[i] = fsolve(func2, initguess, args=((pars.vc['icy']*0.5+pars.vc['silicates']*0.5)*2,disk, self.rhoint[i]))[0]


        #[25.01.15]
        #Finally, determine the physical mass (massL)
        #we'll assume the mass corredsponding to the initial radius Rdi 
        #sets the contribution from the first species ("0") in the 
        #composition list 
        self.massL = self.rhocompos[0] * 4/3*Rdi**3*np.pi /self.fcomp[:,0]

        #[24.01.01]this is a bit ugly... but necessary for adding particles
        #TBd? could we give these more consistent names?
        self.fcompini = self.fcomp[-1]
        self.mtot1 = self.mtotL[-1] #for adding new particles
        self.mini = self.massL[-1]   #for adding particles

        def dm_dr(m,r):
            Rd = physics.mass_to_radius(m,self.rhoint[-1])
            out = gas.get_key_disk_properties(r, 0.0)

            disk = physics.DISK(*out,r, time, dp.Mcp_t(0.0))
            disk.add_auxiliary()

            userparL = disk.add_uservar (dp.user_add_var())    #variables
            userfuncL = disk.add_userfun (dp.user_add_fun())    #functions only
            userevalL = disk.add_user_eval (dp.user_add_eval()) #evaluations
            
            disk.userprops = userparL+userfuncL+userevalL

            self.get_auxiliary(disk, 0.0, Rd, self.rhoint[-1], r, mode='individual')

            Hd = userfun.H_d(self.St,disk)

            dmdr = 2*np.sqrt(np.pi)*Rd**2*self.sfd/2/Hd  
            return dmdr[-1]

        #massL = np.flip(odeint(dm_dr, self.mini, radL).T)*compmask
        #self.massL = massL[0]

        #TBR
        #self.make_Y2d()   #get a Y2d used to integrate
        for i in range(len(dcomposL)):
            del dcomposL[i]['Z_init']
            del dcomposL[i]['mask_icl']


        #<-- initialization (Superparticles)


    def loc_init (self, specL=[]):
        """
        provides the (initial) locations of particles and midpoints
        - accounting for refinement near special locations

        [25.01.21]: NOTE it is important to specify the locations first
                    and then the midpoints!!
        [25.01.23]: In this scheme it is important that grid points are fixed
                    so self.rinn and self.rout should remain fixed
        """

        ## first get the fixed positions
        xdum = np.log(self.rout/self.rinn)
        npar = int( xdum/self.delta)

        loc = self.rinn *np.exp((0.5+np.arange(npar))*self.delta)

        #insert the iceline point
        val = np.sort(specL)[::-1]
        ixL = list(np.searchsorted(loc, val))

        #this needs to be done in reverse order...
        #we add the iceline as a midpoint. It's possible that a
        #very small particle is created interior to it however... 

        if 'Xspecial' not in pars.dresample:
            nres = 1
        else:
            nres = pars.dresample['Xspecial'] #resolution enhancement at specials

        for i,ix in enumerate(ixL):
            #give finer resolution 
            locadd = loc[ix-1] *np.exp(self.delta*np.arange(1,2*nres)/nres)
            loc = np.concatenate((loc[:ix], locadd, loc[ix+1:]))

        locmid = np.concatenate(([self.rinn], 
                                 np.sqrt(loc[1:]*loc[:-1]), 
                                 [self.rout]))
        return loc, locmid


    def select_single (self, ix):

        kwargs = {}
        # select the properties that are list or numpy.ndarray
        #propL = [attr for attr in dir(self) if not attr.startswith('__') and isinstance(getattr(self, attr), list) or isinstance(getattr(self, attr), np.ndarray)]   
        # propL = ['locL','massL','mtotL','fcomp','St','eta'] maybe just select properties artificially is better
        # select the properties that are float

        for prop in self.propL:
            if len(getattr(self,prop)) == self.num:
                try:
                    kwargs[prop] = getattr(self,prop)[ix]
                except:
                    import pdb;pdb.set_trace()
        
        for prop in self.propSol:
            kwargs[prop] = getattr(self,prop)

        spi = SingleSP (**kwargs)
        return spi


    def make_Y2d (self):
        #let's say the second part in composition is always icy fraction 
        return np.array([self.locL, self.massL])

    def get_rhoint(self):
        "get the true fcomp according to fcomp"
        Volume = np.zeros_like(self.locL)
        for i in range(len(self.rhocompos)):
            Volume[:] += self.fcomp[:,i]/self.rhocompos[i]
        
        self.rhoint=1/Volume


    def get_radius(self):
        "update the rhoint according to new fcomp"
        self.get_rhoint()

        return (self.massL/(self.rhoint*4/3*np.pi))**(1/3)


    def get_auxiliary (self, disk, time, specloc, Rd = None, rhoint = None, loc = None, mode='group'):
        """
        Get the auxiliary properties of particles in disk, which 
        need disk properties in the particles's locations
        """

        #LZX[24.08.29]:this add disk properties things are moved to 
        ##             system.get_disk()

        if Rd is None:
            Rd = self.get_radius()

        if rhoint is None:
            rhoint = self.rhoint

        if loc is None:
            loc = self.locL 

        ##[25.01.20]let's to the if/else in functions

        St, v_r = ff.Stokes_number (disk, Rd, rhoint, Sto=self.stokesOld)

        if False:
            if pars.dragmodel=='Epstein':
                St = physics.Stokes_Epstein (Rd, self.rhoint, disk.vth, disk.rhog, disk.OmegaK)
                St *= np.sqrt(8/np.pi) #difference b/w sound speed and thermal velocity
            else:#default
                #obtain Stokes number by iterating on drag law
                #LZX [24.08.04]: insert the rhoint calculated from particles here
                St, v_r = ff.St_iterate (disk.eta,
                                         disk.vK,
                                         disk.vth,
                                         disk.lmfp,
                                         disk.rhog,
                                         disk.OmegaK,
                                         Rd,
                                         rhoint,
                                         Sto=self.stokesOld)

        if mode=='individual':
            sfd = disk.dot_Md(time)/(-2*loc*np.pi*v_r)
        else:
            if pars.sfdmode=='simple':
                #adds the surface to the particles
                #LZX[24.11.01] this fcomp thing should only be used in the steady mode 
                sfd = ff.sfd_simple (self.mtotL, loc, specloc)#/len(self.fcompini)*np.count_nonzero(self.fcomp, axis=1)
            elif pars.sfdmode=='special':
                sfd = ff.sfd_special (self.mtotL, loc, specloc)
            elif pars.sfdmode=='sfd_chris':
                sfd = ff.sfd_chris (self.mtotL, loc)
            elif pars.sfdmode=='sfd_spline':
                sfd = ff.sfd_spline (self.mtotL, loc)
            elif pars.sfdmode=='fixed_bin':
                sfd = ff.sfd_fixedbin (self.mtotL, loc, self.pgrid, specloc)
            elif pars.sfdmode=='steady':
                sfd = disk.dot_Md(time) /(-2*loc*np.pi*v_r)/len(self.fcompini)*np.count_nonzero(self.fcomp, axis=1) #v_r<0
                #sfd1= disk.dot_Md(time) /(-2*self.locL*np.pi*v_r)
                #import pdb;pdb.set_trace()
            else:
                sfd = None
                

        dauxi = {'Rd':Rd, 'St':St, 'v_r':v_r, 'mcp':disk.mcp, 'Hg':disk.Hg, 'sfd':sfd, 'temp':disk.temp} 
        for key in disk.userprops:
            dauxi[key] = getattr(disk, key)

        for key,val in dauxi.items():
            setattr(self, key, val)


    def dY2d_dt (self,Y2d,time):
        """
        input:
            Y2d -- state vector
            time -- time
            disk -- disk object
        """

        #unpack the state vector
        loc, mphy = Y2d
        self.stokesOld = self.St 


        #assume the relative velocity to be the half of radial velocity
        #v_dd = np.abs(v_r)/2    



        drdt = self.v_r
        
        #provide the composition as an argument (in general way)
        #[24.08.05]LZX: the Hd and delv have been integrated into dm_dt
        dmdt = userfun.dm_dt (self)

        Y2ddt = np.zeros_like(Y2d)
        Y2ddt[0] = drdt
        Y2ddt[1] = dmdt

        #[24.01.05]:also return additional particle properties
        #if returnMore:
            #[24.01.07]CWO: alpha cannot be returned here, b/c some disk don't have it!
            #
        #    dMore = {'Rd':Rd, 'St':St, 'v_r':v_r, 'mcp':disk.mcp, 'Hg':disk.Hg} 
        #    for key in userparL+userevalL:
        #        dMore[key] = getattr(disk,key)

        #    return Y2ddt, dMore

        #else:
        #    return Y2ddt  
        return Y2ddt
    
    def remove_particles (self,remove_idx):
        
        mrem = np.sum(self.mtotL[remove_idx])
        self.locL = np.delete(self.locL, remove_idx)
        self.massL = np.delete(self.massL, remove_idx)
        self.mtotL = np.delete(self.mtotL, remove_idx)
        self.fcomp = np.delete(self.fcomp, remove_idx, axis = 0)
        #for prop in self.propL:
        #    pL = getattr(self, prop)
        #    if len(pL)==self.num:
        #        pL = np.delete (pL , remove_idx, axis=0)
        #        setattr(self, prop, pL)

        nrem = len(remove_idx)
        self.num -= nrem

        

        return nrem,mrem


    def add_particles (self,Nadd):

        self.locL = np.append(self.locL, self.rout)
        self.massL = np.append(self.massL, self.mini)
        self.mtotL = np.append(self.mtotL, self.mtot1)  #mtot1 is needed here
        self.fcomp = np.append(self.fcomp, [self.fcompini], axis=0) #[24.01.01] added
        # because we need to use these properties somewhere, so these properties also needed to be postprocess, but for now this is very crude
        #self.Rd = np.append(self.Rd, self.Rdi)
        #self.St = np.append(self.St, self.St[-1])
        #self.v_r = np.append(self.v_r, self.v_r[-1])
        #self.Hg = np.append(self.Hg, self.Hg[-1])
        #self.eta = np.append(self.eta, self.eta[-1])
        #self.mg = np.append(self.mg, self.mg[-1])
        self.num += 1 

        if Nadd!=1:
            #[24.01.07]CWO: I dont understand why we need to add 2 particles sometimes?
            print('[Super.add_particles]WARNING:can only add 1 particle // reduce timestep')


        # For the situation where not only one particles will be added, maybe useful in future
             
        #for moment put the implemented particles all at the outmost edge
        #lociL=np.linspace(self.rout,self.rout,add_number)
        #miL=np.linspace(self.mini,self.mini,add_number)
        #mtot=np.nanmean(self.Y2d[2])+add_number*self.mini
        #mtotL=np.linspace(mtot,mtot,add_number)
        #self.Y2d=np.append(self.Y2d,np.array([lociL,miL,mtotL]),1)
        #self.Y2d[-1]=mtot

        return


class PLANET ():
    """
    Copied over from /NewLagrange
    """
    def __init__(self, time, loc, mplanet, fcomp):
        self.loc = loc          #location
        self.starttime = time        #time when it appears
        self.mass = mplanet     #its mass
        self.fcomp = fcomp      #its composition
        self.spCrossTime = [0.0]   #list when particles cross
        self.spCrossMass = [0.0]   #corresponding mass
        self.spCrossTau = [-1]   #corresponding stokes number
        self.ncross = 0
        self.planetMassData = [[time, mplanet]]
        self.relp_mass = np.inf
        #self.compData = [[time,fcomp[0]]]
        self.relp_comp = np.nan
        self.max_jumpT = 0.0
        self.dlocdt = 0.0
        self.acc_rate = 0.0
        #the fit paremeter in the mass growth 
        self.oldp0 = np.array([1,1])
        
        self.accretion = True

    def record_cross_sp (self, time, sp, idxcross):
        for k in idxcross:
            self.spCrossTime.append(time)
            self.spCrossMass.append(sp.msup[k])
            self.spCrossTau.append(sp.tau[k])
            self.ncross += 1

    def calc_mdot (self, time, Nast=15):
        tlba = time -np.array(self.spCrossTime)[::-1]   #lookback time
        tast = tlba[:Nast].max()                            #characteristic timescale
        mdotarr = np.array(self.spCrossMass)[::-1] /tast
        wi = np.exp(-tlba/tast)     #weights
        mdot = np.sum(mdotarr *wi)  #mass flux through iceline
        return mdot

    def get_maxJumpT(self, Npts, Nfit):
        if Npts >= Nfit and self.accretion:
            #better way to do
            timedots, massdots = np.log(np.array(self.planetMassData).T)
            #timedots = np.log10([self.selfMassData[j][0] for j in range(len(self.selfMassData))])
            def mass_fit(logt,p,b):
                logm = p*logt+b
                return logm 
            def jac_mass_fit (logt,p,b):
                jac = np.array([logt, np.ones_like(logt)]).T
                return jac

            #massdots = np.log10([self.selfMassData[j][1] for j in range(len(self.selfMassData))])

            relp = np.inf #relative error/uncertainty in pwl index
            #popt = self.oldp0#np.array([1,1]) #initial guess (p0)
            #mdots = massdots[-Nfit:]
            #tdots = timedots[-Nfit:]
            #k = (mdots[-1]-mdots[0])/(tdots[-1]-tdots[0])
            #b = mdots[-1]-k*tdots[-1]
            popt = np.array([0.,0.]) #initial guess (p0)
            
            while Nfit<=Npts:
                ##[24.01.20]CWO: this may help a bit, but can still be much faster b/c of linear square root

                mdots = massdots[-Nfit:]
                tdots = timedots[-Nfit:]
                mmean = np.mean(mdots)
                tmean = np.mean(tdots)

                #popt[0] = (np.sum((mdots*tdots))-Nfit*mmean*tmean)/(np.sum(tdots**2)-Nfit*tmean**2)
                den = np.sum((tdots-tmean)**2) #ensures it's positive
                nom = np.sum((mdots-mmean)*(tdots-tmean))
                popt[0] = nom/den
                popt[1] = mmean - popt[0]*tmean
                

                #Then we get the error of the 1st fit parameter
                sse = np.sum((mdots-(popt[0]*tdots+popt[1]))**2)
                kcov = np.sqrt(sse/((Nfit-2)*np.sum((tdots-tmean)**2)))
                #import pdb;pdb.set_trace()

                #try:
                #    popt, pcov = sciop.curve_fit(mass_fit, timedots[-Nfit:], massdots[-Nfit:], p0=popt)
                #except:
                #    import pdb;pdb.set_trace()
                #popt, pcov = sciop.curve_fit(mass_fit, timedots[-Nfit:], massdots[-Nfit:], p0=popt, jac=jac_mass_fit)

                #line = '{:7.4f} {:10.3e} {:10.3e} {:6d}'.format(popt[0], np.sqrt(pcov[0,0]), np.sqrt(pcov[0,0])/popt[0], Nfit)
                #print(line)
                #[15.05.2024]cwo: sometimes pcov is seen to evaluate to infinity due to bad initial popt.
                #                 so retry
                if kcov==np.inf or popt[0] ==np.inf:
                    import pdb;pdb.set_trace()
                    #the assumption/hope is that popt is OK
                elif np.sqrt(kcov)/np.abs(popt[0])<relp:
                    pidx = popt[0]
                    psig = np.sqrt(kcov)
                    relp = psig/np.abs(pidx)
                    Nfit = int(Nfit*1.5) +1 #[24.02.01]cwo increased to 1.5 to prevent noise
                else:
                    break


            self.relp_mass = relp
            self.oldp0 = popt
            #if Npts>=10 and self.ntime>1000:
            #    import pdb; pdb.set_trace()
            #plt.scatter(timedots, massdots)
            #t_list=np.linspace(timedots[0], timedots[-1], 30)
            #plt.plot(t_list, mass_fit(t_list, *popt))
            #plt.savefig('/home/lzx/CpdPhysics/Test/Zhixuan/test.jpg')
            #print(self.dmdt/3e23) 
            self.dmdt = pidx *self.mass/(self.time) 

            self.dmdt_err = abs(psig *self.mass/self.time) ##LZX: please check expressions
            #PmassTscale[i] = 1/abs(pidx)*(np.exp(timedots[-1]) - self.starttime)
            self.massTscale = self.mass/self.dmdt
            #self.tmass_err = abs(psig/popt[0] *PmassTscale[i])
            
            #jump time is limited by uncertainty in the fit
            denom = (self.dmdt_err - self.dmdt*pars.jumpfracD['thre_jump_max'])
            if denom<0:
                self.max_jumpT = np.inf
            else:
                self.max_jumpT = pars.jumpfracD['thre_jump_max']*self.mass /denom
                #print(f'{self.ntime} {i} {Nfit} {Npts} {self.max_jumpT/cgs.yr:10.3e}')
            #TBD: other conditions (migrate into icelines or no gas region...)
        elif self.accretion == False:
            self.dmdt = 0
            self.relp_mass = np.inf
            self.max_jumpT = np.inf
            self.massTscale = np.inf
        else:
            self.dmdt = 0
            self.relp_mass = np.inf 
            self.max_jumpT = 0.0 #np.nan
            self.massTscale = np.inf


        return self.max_jumpT
            


#class for the properties of the central object
class CentralBody (object):
    """
    Get properties of central body, star for PPD or planet for CPD
    """
    def __init__(self, time, mcp0, rho=cgs.rhoRJ):
        self.m = mcp0 
        #TBD: Getting radius should be user-defined, in the userfun.py.
        #If it's missing in the userfun, the radius should be None
        self.rho = rho 
        self.r = physics.mass_to_radius(mcp0,rho)
        self.time = time 

        #make here a function to be more general
        try:
            self.Mt = dp.Mcp_t
        except:
            def constant_m(time):
                return mcp0 
            self.Mt = constant_m

    def get_mass (self, time =None):
        if time is None:
            time = self.time 

        return self.Mt(time)

    def update (self, time):
        self.m = self.Mt(time)
        self.r = physics.mass_to_radius(self.m, self.rho)
        self.time = time

    def magneto_radius (self, dotMg, gas, B_mag=None):
        """
        Get the radius of magneto sphere, 
        only the fully ionized disk will have cavity. 
        
        Input parameters:
        B_mag: the magnetic field of central body, can be default 
        dotMg: the gas mass inflow rate. 
        gas: gas object
        """
        r_cav = userfun.magneto_radius(B_mag, self.rho, dotMg, gas, self.r, self.m, self.time)
        return r_cav

#make here a function to be more general

class ICELINE(object):
    def __init__(self,species,temp):
        self.species=species
        self.temp=temp
        #self.frac=frac
    
    def find_iceline (self,rad, time, gas):
        Tice = self.temp
        Tdisk = gas.get_key_disk_properties(rad,time)[1]
        return Tdisk -Tice

    def get_icelines_location (self,gas,time,bounds=None,guess=None):
        """
        get location of iceline, whose temperature is assumped as 160K
        """

        if guess is not None:
            dsol = sciop.root_scalar(self.find_iceline, x0=guess, args=(time,gas), 
                            method='secant', rtol=1e-6)
            self.loc = dsol.root
            return

        #change the bounds to make it general
        #[23.01.08]LZX: if we change alpha larger,there exist possibility that can't find the iceline location
        #               so for now make a try-exception here, if can't find, then set it to np.nan
        try:
            self.loc = sciop.brentq(self.find_iceline, *bounds, args=(time,gas))
        except:
            self.loc = np.nan
