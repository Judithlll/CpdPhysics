import cgs
import numpy as np
import subprocess as sp
import physics
import pylab as pl
import parameters as pars
import matplotlib.pyplot as plt

def init_compos (material):
    dcompos = {}
    if material=='silicates':
        dcompos['Zinit'] = 0.01
    elif material=='H2O':
        dcompos['name'] = 'H2O'
        dcompos['Zinit'] = 0.01
        dcompos['iceline'] = True
        dcompos['rhoint'] = 1.0
        dcompos['iceline_temp'] = 200
        #dcompos['iceline_init'] = 10*cgs.RJ ## do we need this?

    return dcompos

plotnum = 1
class PlotObj (object):

    def __init__(self):
        self.fg, self.axL = pl.subplots(2,2, figsize=(8,6))
        self.fg.subplots_adjust(bottom=0.05,left=0.07,right=0.93)

        axa, axb, axc, axd = self.axL.ravel()
        axa.set_ylim(0.03, 500)
        axb.set_ylim(1e-4, 1e1)
        axc.set_ylim(2e25, 1e28)
        axd.set_ylim(1e-8, 2e5)
        axa.set_ylabel('surface density')
        axb.set_ylabel('Stokes number')
        axc.set_ylabel('total mass')
        axd.set_ylabel('particle mass')

        axb.tick_params(axis='y', labelleft=False, labelright=True)
        axd.tick_params(axis='y', labelleft=False, labelright=True)

        for ax in [axa,axb,axc,axd]:
            ax.set_xlim(0.95*pars.dgasgrid['rinn']/cgs.au,pars.dgasgrid['rout']/cgs.au)
            #ax.set_xlim(1,5)

        for ax in [axb,axd]:
            ax.yaxis.set_label_position("right")

    def add_lines (self, system, iplot):
        axa, axb, axc, axd = self.axL.ravel()

        locL = system.particles.locL
        massL = system.particles.massL
        mtotL = system.particles.mtotL

        if iplot>1:
            for line in [self.linea, self.lineb, self.lined, self.linec]:
                line.pop(0).remove()

            #or more easily...
            for line in self.lineL:
                line.remove()

        self.lineL = []
        self.linea = axa.loglog(locL/cgs.au, system.particles.sfd, 'b.', ms=1, color='b')
        self.lineb = axb.loglog(locL/cgs.au, system.particles.St, color='r')
        self.lined = axc.loglog(locL/cgs.au, mtotL, '.', ms=1, color='k')
        self.linec = axd.loglog(locL/cgs.au, massL, '.', ms=1, color='k')

        for ax in [axa,axb,axc,axd]:
            line = ax.axvline(pars.dgasgrid['rinn']/cgs.au, color='k', lw=0.3)
            self.lineL.append(line)

        spatch = axa.text(0.1,0.1,'{:11.2e} yr'.format(system.time/cgs.yr), transform=axa.transAxes)

        for iceline in system.icelineL:
            for ax in [axa,axb,axc,axd]:
                line = ax.axvline(iceline.loc/cgs.au, color='k', lw=0.3)
                self.lineL.append(line)

        self.fg.savefig(f'data/plot{iplot:05d}.png')
        spatch.remove()


def do_stuff (system, init=False, final=False):
    global plotnum, plotobj

    plottimeL = np.array([0, 1e1, 2e1, 5e1, 1e2, 2e2, 3e2, 1e3, 1.4e3, 2e3, 
                          2.1e3, 2.2e3, 3e3, 4e3, 5e3, 5.3e3, 6e3, 7e3, 1e4, 2e4, 5e4, 1e5, 1.5e5, 2e5, 5e5, 1e6, 2e6, np.inf]) *cgs.yr

    if init:
        plotobj = PlotObj ()
    else:
        tkeyL = system.minTimes.nameL
        tminarr = system.minTimes.tminarr

        #partices drift/growth/rel.motion
        imin = system.minTimes.dpart['imin']

        sfmt = '{:8d} {:5d} {:10.2e} {:3d} {:2d} {:2d} {:10.2e}'
        line = sfmt.format(system.ntime, len(system.particles.massL), system.deltaT, 
                                            tminarr.argmin(), imin[0],imin[1], system.time/cgs.yr)


        output = sp.run('tail -n1 log/system_evol.log', shell=True)
        if system.time > plottimeL[plotnum]: #plot every 1 yr
            #plot_sfd(system.particles.locL, system.particles.sfd, system.time, system.minTimes.dpart['imin'], system.deltaT, system.timeL, system.resam_time)
            #my_plot(system, plotnum)
            plotobj.add_lines(system, plotnum)
            plotnum += 1




def add_planet(system):

    #First modify the candidate 

    #Then add the planet to the planet list
    if len(system.planet_candidate)>0:
        for planet in system.planet_candidate:
            if planet.starttime <= system.time:
                system.add_planet(planet)
                system.planet_candidate.remove(planet)

    #Finally we should sort the planets according to the location.

    return system


def del_v (St, particles):

    #turbulent relative veolocity
    cs = np.sqrt(cgs.kB*particles.temp/particles.mu/cgs.mp)
    vt = np.sqrt(3*particles.alpha*St)*cs

    return (vt**2 + (particles.v_r/2)**2)**0.5

def H_d (St, particles):
    return physics.H_d(particles.Hg, St, particles.alpha) 


def dm_dt (particles):

    Rd = particles.Rd 
    sigD = particles.sfd
    St = particles.St 
    fcomp = particles.fcomp
    delv = del_v(St, particles)
    Hd = H_d(St, particles)


    dmdt = Rd**2 *delv *sigD/Hd   #eq. 5 of Shibaike et al. 2017

    return 1e-1*dmdt


#def Stokes_number(v_r, Rd, v_th, lmfp, Omega_K, rho_g, rhoint):
#    return Rd*rhoint /(rho_g*v_th) *Omega_K

def Stokes_number (**kwargs):
    return physics.Stokes_Epstein(**kwargs)

#def Stokes_number (**kwargs):
#    return physics.Stokes_general(**kwargs)


#Let's make a function that show the comparison before and after the resample: 
def ba_resample(loc,locn, mphy, mphyn, mtot, mtotn, isL, imL, time):
    plt.close()
    fig, axL = plt.subplots(2, max(len(isL),len(imL)), figsize=(8,6)) 

    if len(axL.shape)==1: 
        axL = axL.reshape(2,1)

    #make the size of ticks smaller for all the subplots 
    for ax in axL.ravel():
        ax.tick_params(axis='both', which='major', labelsize=6)

    baloc =  loc.copy()/cgs.au 
    balocn = locn.copy()/cgs.au 
    
    #plot the particles around the split and merge locations 
    for id, iss in enumerate(isL):
        sidx = np.arange(max(iss-3,0),min(iss+3,len(loc)-1))
        axL[0, id].scatter(baloc[sidx], mphy[sidx], color='b', s=1)
        axL[0, id].set_xlim(baloc[sidx[0]], baloc[sidx[-1]])
        axL[0, id].set_ylim(0.8*min(mphy[sidx]), 1.1*max(mphy[sidx])) 

        sidx = np.append(sidx, sidx[-1]+1)
        axL[0, id].scatter(balocn[sidx], mphyn[sidx]-mphy[iss]/10, color='r', s=1) 

    for id, imm in enumerate(imL): 
        midx = np.arange(max(imm-3,0),min(imm+3,len(loc)-1)) 
        axL[1, id].scatter(baloc[midx], mphy[midx], color='b', s=1)
        axL[1, id].set_xlim(baloc[midx[0]], baloc[midx[-1]])
        axL[1, id].set_ylim(0.8*min(mphy[midx]), 1.1*max(mphy[midx]))

        midx = midx[:-1]
        axL[1, id].scatter(balocn[midx], mphyn[midx]-mphy[imm]/10, color='r', s=1)

    plt.savefig('./bc_rasample/'+'{:.2f}'.format(time/cgs.yr)+'.jpg')

    return
