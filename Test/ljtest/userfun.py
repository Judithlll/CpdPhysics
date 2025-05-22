import cgs
import matplotlib.pyplot as plt
import parameters as pars 
import numpy as np
import subprocess as sp
import physics
import pylab as pl
import os 
import disk_properties as dp

Z0 = 0.01
eps_d = 0.5/10 #Lamb&Jhon 2014 eq.9
def init_compos_Z (material):
    """
    set the initial composition of the disk Z

    this is a 2D array, such that Z(i,j) gives
    the composition of particle i for species j
    """

    #repackage
    rout = 1e3*cgs.au
    if material=='silicates':
        def f_compos (rad):
            Zsil = Z0*np.exp(-(rad/rout)**4)
            return Zsil

    return f_compos

def r_g(t):
    """
    pebble producing line 
    """
    return (3/16)**(1/3)*(cgs.gC*cgs.Msun)**(1/3)*(eps_d*Z0)**(2/3)*t**(2/3)

def add_planet(system):
    return system

def init_compos (compos):
    dcompos = {}
    if compos=='silicates':
        dcompos['Z_init'] = init_compos_Z (compos)

    return dcompos

def do_stuff (system, init=False, final=False):
    global iplot, tplot, fg, ax, figname, colL

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
        print('{:.2e} '.format(np.max(system.particles.sfd)), np.argmax(system.particles.sfd))

        #sfmt = '{:10.3e} {:10.3e} {:10.3e}'        
        #line = sfmt.format(system.time, system.particles.v_r[-1], system.particles.locL[-1]/cgs.au)
        #print(line)

        #if len(system.messages.msgL)>0: import pdb; pdb.set_trace()

        if system.time>=tplot[iplot]:
            print('should plot stuff')
            iplot += 1

            ax.loglog(system.particles.locL/cgs.au, system.particles.sfd, 
                        ms=2, lw=0.5, c=colL[iplot])

            sigmaG = system.gas.sigmaG
            ax.loglog(system.particles.locL/cgs.au, sigmaG, c='k', lw=0.5) 
            axst.loglog(system.particles.locL/cgs.au, system.particles.St, c = colL[iplot], lw=0.5)

            #sigana = sigma_rt (system.particles.locL, -system.particles.v_r, system.time, system.gas)

            #ax.loglog(system.particles.locL/cgs.au, sigana, c='k', lw=0.5)
            #ax1.semilogx(system.particles.locL/cgs.au, system.particles.sfd/sigana-1, c=colL[iplot], lw=0.5)
            #ax2.semilogx(system.particles.locL/cgs.au, system.particles.sfd/sigana-1, c=colL[iplot], lw=0.5)

            fg.savefig(figname, dpi=180)
            #if iplot>=3: import pdb; pdb.set_trace()



def plot_sfd(locL,sfd,time,imin,deltaT,timeL,restime):
    plt.figure(figsize=(6,4))
    plt.subplot(211)
    plt.xlim(0,800)
    plt.ylim(1e-11, 1e3)
    plt.title('Surface density profile at {:.2f}yr'.format(time/cgs.yr))
    plt.loglog(locL/cgs.au, sfd, '.-', label=str(imin)+'\n'+'{:.2f}'.format(deltaT))
    plt.scatter(locL[imin[1]]/cgs.au, sfd[imin[1]], c= 'red')
    plt.axvline(5.89, linestyle='dashed', color='black', linewidth = 1)
    plt.legend(loc='lower right')

    plt.subplot(212)
    plt.xlim(time/cgs.yr-500,time/cgs.yr+500)
    plt.xticks([time/cgs.yr], ['{:.2f}'.format(time/cgs.yr)])
    plt.ylim(1e2,1e9)
    plt.yscale('log')
    plt.plot(np.array(timeL)/cgs.yr, np.append(np.diff(timeL), deltaT) )
    for t in restime: 
        plt.axvline(t/cgs.yr, linestyle='dashed', color='black', linewidth = 1)
    plt.scatter(time/cgs.yr, deltaT, c='red')

    plt.savefig('./sfdevol/{:.2f}.png'.format(time))
    plt.close()

def del_v (St, disk):
    #from J&L 2014
    cs = np.sqrt(cgs.kB*disk.temp /(disk.m_gas*cgs.mp))
    return np.sqrt(3*disk.alpha*St)*cs

def H_d (St, disk):

    #from J&L 2014
    Hd = disk.Hg*np.sqrt(disk.alpha/St)
    
    return Hd


def dm_dt(particles):
    
    # from J&L 2014
    Rd = particles.Rd 
    sigD = particles.sfd
    St = particles.St 
    delv = del_v(St, particles)
    Hd = H_d(St, particles)

    sigmacol = 4*np.pi*Rd**2

    return sigD/(np.sqrt(2*np.pi)*Hd) *delv*sigmacol/4


def make_animation(mp4name, path='./plot/satepart_splitmerge'):
    import cv2
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

tplot = np.array([0,1e3,1e4,1e5,1e6,1e7, 1e99]) *cgs.yr
colL = ['k', '0.8', '0.6', '0.4', 'r', 'm', 'b', 'g', (0.8,0.4,0.),'y']
iplot = 0

fg, (ax,axst) = pl.subplots(2,1, figsize=(6,9), sharex=True, gridspec_kw={'height_ratios':[1,1]})

for aa in [ax, axst]:
    aa.set_xlim(0.4, 250)

rg = r_g(tplot[:-1])

ax.set_ylim(1e-4, 1e3)
ax.set_xlim(1, 1e3)
for r in rg: 
    ax.axvline(r/cgs.au, linestyle='dashed', color='gray', linewidth = 1)

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
        axL[0, id].set_ylim(0.3*min(mphy[sidx]), 1.1*max(mphy[sidx])) 

        sidx = np.append(sidx, sidx[-1]+1)
        axL[0, id].plot(balocn[sidx], mphyn[sidx],'x', color='r', label='insert particles at {:.4f}'.format(balocn[iss+1])) 
        axL[0, id].legend(fontsize=6)

    

    for id, imm in enumerate(imL): 
        midx = np.arange(max(imm-3,0),min(imm+3,len(loc)-1)) 
        axL[1, id].scatter(baloc[midx], mphy[midx], color='b', s=1)
        axL[1, id].set_xlim(baloc[midx[0]], baloc[midx[-1]])
        axL[1, id].set_ylim(0.8*min(mphy[midx]), 1.1*max(mphy[midx]))

        midx = midx[:-1]
        axL[1, id].plot(balocn[midx], mphyn[midx], 'x', color='r')

    plt.savefig('./ba_resample/'+'{:.2f}'.format(time/cgs.yr)+'.jpg')

    return

def check_split(sim, idx=None):
    #here we will make a space-time diagram with iso-mphy lines,to show what's the matter with the split 
    
    if idx is None:
        idxL = np.arange(0, 20)
    else:
        idxL = np.arange(max(idx-5,0),min(idx+10,len(sim.particles.locL)-1))
    loc = sim.particles.locL/cgs.au 
    mphy = sim.particles.massL 
    sfd = sim.particles.sfd 
    vr = sim.particles.v_r
    mtot = sim.particles.mtotL

    plt.close()
    fig, axL = plt.subplots(1,4, figsize=(16,4)) 

    #make the size of ticks smaller for all the subplots and set the limits 
    titles = ['mphy', 'sfd', 'vr', 'mtot']
    for ax in axL.ravel():
        ax.tick_params(axis='both', which='major', labelsize=6) 
        ax.set_xlim(loc[idxL[0]]*0.8, loc[idxL[-1]]*1.2) 
        ax.set_title(titles.pop(0))
        ax.axvline(sim.rinn/cgs.au, linestyle='dashed', color='black', linewidth = 1)

        for id in idxL:
            ax.axvline(loc[id], linestyle='dotted', color='gray', linewidth = 0.5)
        

    #plot the mphy with scatters 
    axL[0].scatter(loc[idxL], mphy[idxL], color='b', s=3) 
    axL[0].set_ylim(0.8*min(mphy[idxL]), 1.1*max(mphy[idxL])) 
    if idx is not None:
        axL[0].scatter(loc[idx], mphy[idx], color='r', s=5) 

    #plot the sfd with scatters 
    axL[1].scatter(loc[idxL], sfd[idxL], color='b', s=3) 
    axL[1].set_ylim(0.8*min(sfd[idxL]), 1.1*max(sfd[idxL])) 
    if idx is not None:
        axL[1].scatter(loc[idx], sfd[idx], color='r', s=5)


    #plot the vr with scatters 
    axL[2].scatter(loc[idxL], vr[idxL], color='b', s=3)
    axL[2].set_ylim(1.1*min(vr), 0.8*max(vr)) 
    if idx is not None:
        axL[2].scatter(loc[idx], vr[idx], color='r', s=5)

    #plot the mtot with scatters 
    axL[3].scatter(loc[idxL], mtot[idxL], color='b', s=3)
    axL[3].set_ylim(0.8*min(mtot[idxL]), 1.1*max(mtot[idxL]))
    if idx is not None:
        axL[3].scatter(loc[idx], mtot[idx], color='r', s=5)

    #remove the blank space 
    plt.tight_layout()

    plt.savefig('./split_check/'+'{:.2f}'.format(sim.time/cgs.yr)+'.jpg')
    plt.close()
    #plot the particles before and after the split around the split location 

    return

def loc_face(loc, face, msup,time ,mode='zi'): 
    plt.close() 
    fig, ax = plt.subplots(1,1, figsize=(5,2)) 
    ax.scatter(loc/cgs.au, msup, color='b', s=1) 
    ax.set_xscale('log')

    #set the fontsize of all ticks 
    ax.tick_params(axis='both', which='both', labelsize=6)

    if mode=='zi':
        ax.set_xlim(0.095, 0.125)
        ax.set_ylim(1e24,2e25)
    elif mode =='z':
        ax.set_xlim(0.095, 200)
        ax.set_ylim(1e23, 1e26)
    elif mode =='sfd':
        ax.set_xlim(0.095, 10)
        ax.set_ylim(1e-1, 100)
        ax.set_yscale('log')
    elif mode =='sfdzi':
        ax.set_xlim(0.095, 0.125) 
        ax.set_ylim(1e-1, 100)
        ax.set_yscale('log')
    elif mode == 'mphy':
        ax.set_xlim(0.095, 0.125)
        ax.set_ylim(1e-1, 1e8)
        ax.set_yscale('log')
    elif mode == 'normal':
        pass

    #for f in face: 
    #    ax.axvline(f/cgs.au,linestyle='dashed', color='black', linewidth = 0.1)

    ax.vlines(face/cgs.au, 0, 1e28, colors='gray', linestyles='dashed', linewidth=0.1)
    
    plt.tight_layout()
    plt.savefig('./face_loc/'+'{:.2f}'.format(time/cgs.yr)+'.jpg')
    plt.close()

    return

def make_animation(mp4name, path='./plot/satepart_splitmerge'):
    import cv2
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
    return
