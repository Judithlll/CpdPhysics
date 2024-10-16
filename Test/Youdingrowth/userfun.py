import cgs
import matplotlib.pyplot as plt
import parameters as pars 
import numpy as np
import subprocess as sp
import physics
import pylab as pl
import os 
import cv2
import disk_properties as dp

def init_compos_Z (material):
    """
    set the initial composition of the disk Z

    this is a 2D array, such that Z(i,j) gives
    the composition of particle i for species j
    """

    #repackage
    rc = 200*cgs.au
    if material=='silicates':
        def f_compos (rad):
            Zsil = 0.007  *np.exp(-(rad/rc)**2)
            return Zsil

    return f_compos

def add_planet(system):
    return system

def init_compos (compos):
    dcompos = {}
    if compos=='silicates':
        dcompos['Z_init'] = init_compos_Z (compos)

    return dcompos


def g_r (rad, gas, d=1.5):
    Z_0 = init_compos_Z ('silicates')
    Z0 = Z_0(rad)

    sig0, temp, mmw = gas.get_key_disk_properties(rad,0.0)
    #sig0, mmw = disk.sigma_gas_ini (rad)
    #sig0 = Z0*disk.sigma_gas_ini (rad)
    return rad**(d+1) *sig0*Z0


def sigma_rt (rad, vdr, time, gas, d=1.5):
    """
    Youdin & Shu solution
    """
    ri = rad *(1 -(d-1)*vdr*time/rad)**(-1/(d-1))
    sig = rad**(-d-1) *g_r(ri,gas)
    return sig



def do_stuff (system, init=False, final=False):
    global iplot, tplot, fg, ax, figname

    if init:
        pass
    else:
        tkeyL = system.minTimes.nameL
        tminarr = system.minTimes.tminarr

        #partices drift/growth/rel.motion
        imin = system.minTimes.dpart['imin']

        sfmt = '{:8d} {:5d} {:10.2e} {:3d} {:2d} {:2d} {:10.2e}'
        line = sfmt.format(system.ntime, len(system.particles.massL), system.deltaT, 
                                            tminarr.argmin(), imin[0],imin[1], system.time/cgs.yr)

        #output = sp.run('tail -n1 log/system_evol.log', shell=True)
        print(line)

        #sfmt = '{:10.3e} {:10.3e} {:10.3e}'        
        #line = sfmt.format(system.time, system.particles.v_r[-1], system.particles.locL[-1]/cgs.au)
        #print(line)

        #if len(system.messages.msgL)>0: import pdb; pdb.set_trace()

        if system.time>=tplot[iplot]:
            colL = ['k', '0.8', '0.6', '0.4', 'r', 'm', 'b', 'g', (0.8,0.4,0.),'y']
            print('should plot stuff')
            iplot += 1

            ax.loglog(system.particles.locL/cgs.au, system.particles.sfd, 
                        ms=2, lw=0.5, c=colL[iplot])

            axst.loglog(system.particles.locL/cgs.au, system.particles.St, c = colL[iplot], lw=0.5)

            #sigana = sigma_rt (system.particles.locL, -system.particles.v_r, system.time, system.gas)

            #ax.loglog(system.particles.locL/cgs.au, sigana, c='k', lw=0.5)
            #ax1.semilogx(system.particles.locL/cgs.au, system.particles.sfd/sigana-1, c=colL[iplot], lw=0.5)
            #ax2.semilogx(system.particles.locL/cgs.au, system.particles.sfd/sigana-1, c=colL[iplot], lw=0.5)

            fg.savefig(figname, dpi=180)



def plot_sfd(locL,sfd,time,imin,deltaT,timeL):
    plt.figure(figsize=(12,8))
    plt.subplot(211)
    plt.xlim(0,800)
    plt.ylim(1e-11, 1e3)
    plt.yscale('log')
    plt.title('Surface density profile at {:.2f}yr'.format(time/cgs.yr))
    plt.plot(locL/cgs.au, sfd, 'x-', label=str(imin)+'\n'+'{:.2f}'.format(deltaT))
    plt.scatter(locL[imin[1]]/cgs.au, sfd[imin[1]], c= 'red')
    plt.axvline(5.89, linestyle='dashed', color='black', linewidth = 1)
    plt.legend(loc='lower right')

    plt.subplot(212)
    plt.xlim(time/cgs.yr-500,time/cgs.yr+500)
    plt.xticks([time/cgs.yr], ['{:.2f}'.format(time/cgs.yr)])
    plt.ylim(1e5,1e10)
    plt.yscale('log')
    plt.plot(np.array(timeL)/cgs.yr, np.append(np.diff(timeL), deltaT) )
    plt.scatter(time/cgs.yr, deltaT, c='red')

    plt.savefig('./sfdevol/{:.2f}.png'.format(time))
    plt.close()

def del_v (St, disk):
    return np.abs(disk.v_r)/2


def H_d (St, disk):
    return disk.Hg


def dm_dt(particles):
    
    Rd = particles.Rd 
    sigD = particles.sfd
    St = particles.St 
    delv = del_v(St, particles)
    Hd = H_d(St, particles)

    return 2*np.sqrt(np.pi)*Rd**2*delv/Hd*sigD   #eq. 5 of Shibaike et al. 2017

#def Stokes_number(v_r, Rd, v_th, lmfp, Omega_K, rho_g, rhoint):
#    return Rd*rhoint /(rho_g*v_th) *Omega_K

#def Stokes_number (**kwargs):
#    return physics.Stokes_Epstein(**kwargs)

#def Stokes_number (**kwargs):
#    return physics.Stokes_general(**kwargs)
def make_animation(mp4name, path='./plot/satepart_splitmerge'):
    pic_list = []
    pics=os.listdir(path)
    pics_sorted=sorted(pics, key=lambda x: float(x[:-4]))
    frame = cv2.imread(os.path.join(path,pics_sorted[0]))
    height, width, layers = frame.shape
    video_name =  mp4name
    fps=10
    video_codec = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(video_name, video_codec, fps, (width, height))
    for pic in pics_sorted:
        video.write(cv2.imread(os.path.join(path,pic)))

    cv2.destroyAllWindows()

    #cv2.VideoWriter(video_name, video_codec, fps, (width,height))
    # import pdb;pdb.set_trace()
    #for pic in pics_sorted:
    #    im = imageio.imread(path+"/"+pic)
    #    pic_list.append(im)
    #imageio.mimsave(save_name_gif, pic_list, 'GIF', loop=0)

tplot = np.array([0,1e4,2e4,2.5e4,5e4,1e5,5e5,1e6, 1e99]) *cgs.yr
iplot = 0

fg, (ax,axst) = pl.subplots(2,1, figsize=(4,6), sharex=True, gridspec_kw={'height_ratios':[1,1]})

for aa in [ax, axst]:
    aa.set_xlim(0.4, 250)


ax.set_ylim(1e-3, 2e3)
#ax1.set_ylim(-0.1, 0.1)
#ax2.set_ylim(-0.01, 0.01)

axst.set_xlabel('distance')
ax.set_ylabel('surface density')
axst.set_ylabel('Stokes number')
if dp.dot_Mg(0.0) == 0.0:
    name1 = 'noimp'
else:
    name1 = 'imp'+'+'+str(dp.fraction)
figname = name1+'+'+str(pars.dresample['fdelS'])+'+'+str(pars.dresample['fdelM'])+'+'+str(pars.dresample['fdelDM'])+'+'+str(pars.dparticleprops['Rdi'])+'.jpg'

