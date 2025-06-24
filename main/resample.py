import numpy as np
import parameters as pars
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import functions as ff
import physics
from scipy.optimize import fsolve 
import userfun

def face_splitmerge (sim, spN, fchange=0.5, fdelX=1, nsampleX=0, fdelDM=0.0001, full_output = False, **args):
    """
    got the inspiration from Arepo and Chengzhe 
    """

    #the criteria now should be related the width of particles
    loc = spN.locL 
    face = spN.get_face()
    
    #get the width of particles' space 
    fwdel = np.diff(face)/face[:-1]

    fdelS = sim.particles.fwdeli /fchange
    fdelM = sim.particles.fwdeli *fchange

    fdelXarr = np.ones_like(fwdel)
    # judge whether the particles need to be resampled 
    imL, = np.nonzero(fwdel<fdelM*fdelXarr) 

    #dont split where we merge 
    fdelXarr[imL] = np.inf 
    fdelXarr[imL+1] = np.inf 
    fdelXarr[imL-1] = np.inf 
    fdelXarr[0] = np.inf; fdelXarr[-1] = np.inf #don't split 1st/last particles 
    
    isL, = np.nonzero(fwdel>fdelS*fdelXarr) 

    if len(isL)>0: 
        print('Resample happens') 

        mtot = spN.mtotL 
        mphy = spN.massL 

        locn = loc.copy() 

        #only consider split now 
        facen = face.copy()
        #facen = np.insert(facen, isL+1, loc[isL])
        #add face at the center 
        faceaddS = np.sqrt(face[isL]*face[isL+1])
        facen = np.insert(facen, isL+1, faceaddS)


        cglocS = np.sqrt(face[isL]*faceaddS)
        locn[isL] = cglocS 
        addlocS = np.sqrt(faceaddS*face[isL+1])
        locn = np.insert(locn, isL+1, addlocS)

        # sevaral ideas for mtot: 
        mtotn = mtot.copy() 

        # 1. split the mass equally 
        #mtotn[isL] = mtot[isL]/2 
        #mtotn = np.insert(mtotn, isL+1, mtot[isL]/2)

        # 2. split the mass according to the width 
        #frac = (loc[isL]-face[isL])/(face[isL+1]-face[isL])
        #mtotn[isL] = mtot[isL]*frac 
        #mtotn = np.insert(mtotn, isL+1, mtot[isL]*(1-frac))

        # 3. conserve the center of mass 
        #mtotn[isL] = mtot[isL]*(addlocS-face[isL])/(addlocS-cglocS)
        #mtotn = np.insert(mtotn, isL+1, mtot[isL]-mtotn[isL])

        # 4. keep the eom simply 
        # the location should be euqally around the original particle
        # locn = loc.copy()
        # sinterv = np.array([])
        # for i in isL:
        #     if loc[isL]-face[isL]>face[isL+1]-loc[isL]:
        #         sinterv =np.append(sinterv, np.sqrt(face[isL+1]*loc[isL])-loc[isL])
        #     else:
        #         sinterv =np.append(sinterv, loc[isL]-np.sqrt(face[isL]*loc[isL]))
        # locn[isL] = loc[isL] - sinterv 
        # locn = np.insert(locn, isL+1, loc[isL] + sinterv)
        #
        # mtotn[isL] = mtot[isL]/2 
        # mtotn = np.insert(mtotn, isL+1, mtot[isL]/2)

        # 4. according to the sfd 
        newlocs = np.concatenate((cglocS, addlocS))
        sfdn = np.interp(newlocs, loc, spN.sfd)
        
        #mtotn[isL] = sfdn[:len(isL)]*2*np.pi*cglocS*(loc[isL] - face[isL])
        #mtotn = np.insert(mtotn, isL+1, sfdn[len(isL):]*2*np.pi*addlocS*(face[isL+1]-loc[isL]))

        # cgmtotS = sfdn[:len(isL)]*2*np.pi*cglocS*(loc[isL] - face[isL])
        # addmtotS = sfdn[len(isL):]*2*np.pi*addlocS*(face[isL+1]-loc[isL])
        # frac = cgmtotS/(cgmtotS+addmtotS)
        # mtotn[isL] = mtot[isL]*frac 
        # mtotn = np.insert(mtotn, isL+1, mtot[isL]*(1-frac))

        # new idea
        cgmtotS = sfdn[:len(isL)]*2*np.pi*cglocS*(faceaddS - face[isL])
        addmtotS = sfdn[len(isL):]*2*np.pi*addlocS*(face[isL+1]-faceaddS)
        frac = cgmtotS/(cgmtotS+addmtotS)
        mtotn[isL] = mtot[isL]*frac 
        mtotn = np.insert(mtotn, isL+1, mtot[isL]*(1-frac))
        # mphyn = mphy.copy() 
        # mphyn[isL] = mphy[isL] 
        # mphyn = np.insert(mphyn, isL+1, mphy[isL])

        #interpolate the mphyn 
        mphyn = np.interp(locn, loc, mphy)
        # import cgs
        # plt.close()
        # plt.plot(loc[:5]/cgs.au, mtot[:5], 'x-', label='old')
        #
        # sfdnm = ff.sfd_face(mtotn[:6], locn[:6], facen[:7])
        # plt.plot(locn[:6]/cgs.au, mtotn[:6], 'o-', label='new')
        # plt.xlabel('loc')
        # plt.ylabel('sfd')
        #
        # plt.legend()
        # plt.show()
        # import pdb;pdb.set_trace()

        #userfun.ba_resample(loc, locn, mtot, mtotn, mtot, mtotn, isL, imL, sim.time)
    else:
        facen = face.copy() 
        locn = loc.copy()
        mtotn = spN.mtotL.copy() 
        mphyn = spN.massL.copy() 

    if len(imL)>0:
        print('Merge happens')

        #get a new imL from ...n things 
        fwdel = np.diff(facen)/facen[:-1]
        imL = np.argwhere(fwdel<fdelM)[0]

        facenm = facen.copy()
        locnm = locn.copy()
        mtotnm = mtotn.copy()
        mphynm = mphyn.copy()

        #get the cloest particle index 
        imcL = np.array([]).astype(int)
        for im in imL:
            if abs(locn[im]-locn[im+1])<abs(locn[im-1]-locn[im]):
                imcL = np.append(imcL, im+1)
                facenm = np.delete(facenm, im+1)
            else:
                imcL = np.append(imcL, im-1)
                facenm = np.delete(facenm, im)

        locnm[imL] = np.sqrt(locn[imL]*locn[imcL])
        locnm = np.delete(locnm, imcL)

        #deal with the mtot 
        mtotnm[imL] = mtotn[imL] + mtotn[imcL]
        mtotnm = np.delete(mtotnm, imcL) 

        #deal with the mphynm: get the weighted mean mass
        mphynm[imL] = (mphyn[imL]*mtotn[imL] + mphyn[imcL]*mtotn[imcL])/mtotnm[imL] 
        mphynm = np.delete(mphynm, imcL)



    else: 
        facenm = facen.copy()
        locnm = locn.copy()
        mtotnm = mtotn.copy()
        mphynm = mphyn.copy()

    #then do the direct merge 
    #direct merge is decided by the distance b/w particles
    fdelDM = spN.fdeli*fchange/3 
    fdel = np.diff(locnm)/locnm[:-1]
    idmL = np.argwhere(fdel<fdelDM)

    if len(idmL)>0: 
        print('Direct merge happens')

        locnmd = locnm.copy() 

        locnmd[idmL] = np.sqrt(locnm[idmL]*locnm[idmL+1]) 
        locnmd = np.delete(locnmd, idmL+1) 

        facenmd = facenm.copy()
        facenmd = np.delete(facenmd, idmL+1) 

        mtotnmd = mtotnm.copy() 
        mtotnmd[idmL] = mtotnm[idmL] + mtotnm[idmL+1]
        mtotnmd = np.delete(mtotnmd, idmL+1) 

        mphynmd = mphynm.copy()
        mphynmd[idmL] = (mphynm[idmL]*mtotnm[idmL] + mphynm[idmL+1]*mtotnm[idmL+1])/mtotnmd[idmL]
        mphynmd = np.delete(mphynmd, idmL+1)
        
    else:
        locnmd = locnm.copy() 
        facenmd = facenm.copy()
        mtotnmd = mtotnm.copy() 
        mphynmd = mphynm.copy()

    if len(isL)>0 or len(imL)>0 or len(idmL)>0: 
        doResample = True 
    else:
        doResample = False       

    if doResample:
        npar = len(locnmd)

        fcomp = spN.fcomp #composition fraction 
        ncomp = len(fcomp[0]) 

        #now treat the fcomp 
        fcompnmd = np.empty((npar,ncomp))  
        if len(sim.specloc) ==0: 
            fcompnmd[:] = fcomp[0]


        #check the mass conservation 
        delm = (np.sum(mtotnmd) - np.sum(spN.mtotL))/np.sum(spN.mtotL)
        if np.abs(delm)>1e-15:
            print(f'mass conservation error: {delm}')
            import pdb;pdb.set_trace()


        if full_output: 
            return facenmd, locnmd, mtotnmd, mphynmd, fcompnmd, isL, imL 
        else:
            return facenmd, locnmd, mtotnmd, mphynmd, fcompnmd
    else: 
        return None




def v_rad (marr,disk,rhoint,Stg=None,vaim=0):
    """
    obtain radial velocity of particles
    - radarr    :input locations
    - disk      :disk object at locations
    - rhoint    :internal densities
    """

    sarr = physics.mass_to_radius(marr,rhoint)
    St, vr = ff.Stokes_number(disk, sarr, rhoint, Sto=Stg)

    return vr -vaim


def local_splitmerge (sim, spN, **kwargs):
    """
    similar to new_splitmerge_chris
    """
    loc = spN.locL 
    mtot = spN.mtotL
    mphy = spN.massL
    fcomp = spN.fcomp #composition fraction

    ncomp = len(fcomp[0])
    xdel = np.diff(np.log(loc))

    fdelXarr = np.ones_like(xdel)
    fdelS = 2*sim.particles.delta
    isL, = np.nonzero(xdel>fdelS*fdelXarr)
    imL, = np.nonzero(xdel<fdelXarr*sim.particles.delta*2/3)

    #splitting: add the locations
    addlocS = np.sqrt(loc[isL]*loc[isL+1])
    addlocM = np.sqrt(loc[imL]*loc[imL+1])

    if len(isL)>0 or len(imL)>0:
        doResample = True
    else:
        doResample = False

    if doResample:#or len(imL)>0:
        locmidext = locmid_ext (loc)
        cummtot = np.concatenate(([0],np.cumsum(mtot)))

        #merging: remove the locations from loc (TBD)
        locn = loc.copy()
        locn[imL] = addlocM
        locn = np.delete(locn,imL+1)

        #a bit weird
        locn = np.concatenate((locn,addlocS))
        locn.sort()
        npar = len(locn) #new number of particles

        locmidnext = locmid_ext (locn)
        locmidnext[0] = locmidext[0] #hack
        cummtotn = np.interp(locmidnext, locmidext, cummtot)
        mtotn = np.diff(cummtotn)

        #composition... TBD
        fcompn = np.empty((npar,ncomp))
        for k in range(ncomp):
            cummass = np.concatenate(([0], np.cumsum(mtot*fcomp[:,k])))
            cummassn = np.interp(locmidnext, locmidext, cummass)
            fcompn[:,k] = np.diff(cummassn) /mtotn

        #the physical mass
        mphyn = interp_mtot_weighted (locmidnext, locmidext, mphy, mtot, mtotn)

        #print(mtotn)
        #sfdnew = ff.sfd_fixedbin (mtotn, locn, sim.particles.pgrid, sim.specloc)
        return locn, mtotn, mphyn, fcompn
    else:
        return None


def new_splitmerge_chris (sim, spN, fdelS, fdelM=0., fdelX=1, nsampleX=0, fdelDM=0.0001):
    """
    [25.01.18]: new splitmerge based on cumulative mass function
                for the moment w/o special locations functionality
    """
    loc = spN.locL 
    mtot = spN.mtotL
    mphy = spN.massL
    fcomp = spN.fcomp #composition fraction

    ncomp = len(fcomp[0])
    xdel = np.diff(np.log(loc))


    #now the masses
    #ydel = np.diff(np.log(mphy))
    #isL2, = np.nonzero(np.abs(ydel)>0.5*fdelXarr)
    #isL = np.union1d(isL1, isL2)

    #merging
    #fdelXarr[0] = 0; fdelXarr[-1] = 0
    fdelXarr = np.ones_like(xdel)
    imL, = np.nonzero(xdel<fdelXarr*sim.particles.delta*2/3)

    #switch off merging
    imL = np.array([],dtype=np.int64)


    fdelXarr[imL] = np.inf #dont split where we merge
    fdelXarr[imL+1] = np.inf
    fdelXarr[imL-1] = np.inf

    #midpoint locations elligible for splitting (isL)
    #don't split 1st/last particles (i.e., around special locations; TBD)
    #fdelXarr[:2] = np.inf; fdelXarr[-2:] = np.inf
    fdelS = 2*sim.particles.delta
    isL, = np.nonzero(xdel>fdelS*fdelXarr)

    #splitting: add the locations
    addlocS = np.sqrt(loc[isL]*loc[isL+1])
    addlocM = np.sqrt(loc[imL]*loc[imL+1])

    if len(isL)>0:
        doResample = True
    else:
        doResample = False


    if doResample:#or len(imL)>0:
        locmidext = locmid_ext (loc)
        cummtot = np.concatenate(([0],np.cumsum(mtot)))

        #merging: remove the locations from loc (TBD)
        locn = loc.copy()
        locn[imL] = addlocM
        locn = np.delete(locn,imL+1)

        #a bit weird
        locn = np.concatenate((locn,addlocS))
        locn.sort()
        npar = len(locn) #new number of particles

        locmidnext = locmid_ext (locn)
        locmidnext[0] = locmidext[0] #hack
        cummtotn = np.interp(locmidnext, locmidext, cummtot)
        mtotn = np.diff(cummtotn)

        #composition... TBD
        fcompn = np.empty((npar,ncomp))
        for k in range(ncomp):
            cummass = np.concatenate(([0], np.cumsum(mtot*fcomp[:,k])))
            cummassn = np.interp(locmidnext, locmidext, cummass)
            fcompn[:,k] = np.diff(cummassn) /mtotn

        mflux = sim.particles.v_r *sim.particles.sfd *sim.particles.locL
        if np.any(np.diff(mflux[1:20])<0) and False:
            print('mflux not in order')


        if True:
            # this mixing scheme for the physical mass works surprisingly well
            # if the power-law equals 0.4. I have no clue why
            pm = 0.5
            dum = interp_mtot_weighted (locmidnext, locmidext, mphy**pm, mtot, mtotn)
            mphyn = dum**(1/pm)

            #the log-based mixing scheme isn't so good though
            #
            #logmphyn = interp_mtot_weighted (locmidnext, locmidext, np.log(mphy), mtot, mtotn)
            #mphyn = np.exp(logmphyn)

        #the following two mixing schemes do not yet account for merging
        elif False:
            # this is a more physical mixing scheme
            pm = 0.5
            mphyadd = mphy[isL]**pm *mphy[isL+1]**(1-pm)
            mphyn = np.insert(mphy, isL+1, mphyadd)

        else:
            #in this mixing scheme, we aim to insert the "right" velocity
            #unfortunately, we get the same issues as before

            #the desired velocities
            vraim = addlocS/2 *(sim.particles.v_r[isL]/loc[isL] +sim.particles.v_r[isL+1]/loc[isL+1])
            #vraim = (sim.particles.v_r[isL] +sim.particles.v_r[isL+1])/2
            #vraim = -np.sqrt(sim.particles.v_r[isL] *sim.particles.v_r[isL+1])

            ## search for proper particles mass
            disk = sim.get_disk (loc=addlocS)
            rhoint = sim.particles.rhoint[isL]  #this should be changed
            stgarr = sim.particles.St[isL]      #initial guess

            mphyadd = fsolve(v_rad, mphy[isL], args=(disk,rhoint,stgarr,vraim))
            mphyn = np.insert(mphy, isL+1, mphyadd)


        if loc[0]<sim.rinn:
            print('1st particle location too small')
            import pdb; pdb.set_trace()

        #print a warning message when the physical mass is not in order
        #I think that this shouldn't happen in reality...
        if np.all(np.diff(np.log10(mphyn[:100]))<0)==False and False:
            print('physical mass not in order')

        #sfdnew = ff.sfd_special (mtotn, locn, sim.specloc)

        return locn, mtotn, mphyn, fcompn
    else:
        return None

def v2mass (vr, loc, sim):
    """
    [25.01.20]: get the mass from the velocity
    vr: float
    loc: float
    """
    from scipy.optimize import fsolve 

    def func(Rd,disk,rhoint):
        St,v = ff.Stokes_number(disk, Rd, rhoint, Sto =0.03)
        return v 

    def func2(Rdi, vgive, disk, rhoint):
        rere = func(Rdi,disk,rhoint)-vgive
        return rere

    #get the disk object 
    disk = sim.get_disk(loc=loc, time = sim.time)
    #get the mass from the velocity 
    idx = np.searchsorted(sim.particles.locL,loc)
    rhoint = sim.particles.rhoint[idx]

    #use the mass of nearest particle as a initial guess
    Rd = physics.mass_to_radius(sim.particles.massL[idx],rhoint)
    rd = fsolve(func2, Rd, args=(vr,disk,rhoint))[0] 
    mass = physics.radius_to_mass(rd,rhoint)

    return mass

def hd2mass (hd, loc, sim):
    """
    [25.01.20]: get the mass from the dust scale height
    hd: float
    loc: float
    """
    from scipy.optimize import fsolve 

    def func(Rd,disk,rhoint, hg):
        St = ff.Stokes_number(disk, Rd, rhoint, Sto =0.03)[0]
        hdd = physics.H_d(hg, St, sim.particles.alpha)
        return hdd 

    def func2(Rdi, hdgive, disk, rhoint, hg):
        rere = func(Rdi, disk, rhoint, hg) - hdgive
        return rere

    #get the disk object 
    disk = sim.get_disk(loc=loc, time = sim.time)
    #get the mass from the velocity 
    idx = np.searchsorted(sim.particles.locL,loc)
    rhoint = sim.particles.rhoint[idx]

    #use the mass of nearest particle as a initial guess
    Rd = physics.mass_to_radius(sim.particles.massL[idx],rhoint)
    #if the iteration cannot be converged, stop here to check 
    rd, info, ier, msg = fsolve(func2, Rd, args=(hd,disk,rhoint,disk.Hg), full_output=True)
    rd = rd[0] 
    if ier != 1:
        print(f'fsolve error: {msg}')
        import pdb; pdb.set_trace()

    mass = physics.radius_to_mass(rd,rhoint)

    return mass 


def get_addmphy(iL, addloc, sim, loc, prop='rv'):
    """
    prop: the properties you want to get the mass from 
        'rv': radial velocity. 
        'hd': dust scale height.
    """
    addmphy = []
    if prop == 'rv':
        for i, idx in enumerate(iL):
            v1 = sim.particles.v_r[idx] 
            v2 = sim.particles.v_r[idx+1]
            r1 = loc[idx]
            r2 = loc[idx+1]
            vnew = (r1*v2+r2*v1)/2/addloc[i]

            addmphy.append(v2mass(vnew, addloc[i], sim))
    elif prop == 'hd':
        for i, idx in enumerate(iL):
            #we use interpolation to get the hd at the new location 
            hd = physics.H_d(sim.particles.Hg, sim.particles.St, sim.particles.alpha)
            hdnew = np.interp(addloc[i], loc, hd) 

            addmphy.append(hd2mass(hdnew, addloc[i], sim))
    return addmphy

def new_splitmerge_zxl (sim, spN, fchange=0.9, fdelX=1, nsampleX=0, fdelDM=0.0001, full_output = False,specloc=None, **args):
    """
    [25.01.18]: new splitmerge based on cumulative mass function
                for the moment w/o special locations functionality
    """
    loc = spN.locL 

    xdel = np.diff(np.log(loc))

    fdelS = sim.particles.delta*2 
    fdelM = sim.particles.delta/2


    #now the masses
    #ydel = np.diff(np.log(mphy))
    #isL2, = np.nonzero(np.abs(ydel)>0.5*fdelXarr)
    #isL = np.union1d(isL1, isL2)

    fdelXarr = np.ones_like(xdel)
    #don't resample the particles around the special locations
    if specloc is not None: 
        specidx = np.searchsorted(loc, specloc) 
        fdelXarr[specidx] = np.nan
        fdelXarr[specidx-1] = np.nan 

    #merging
    #don't merge 1st/last particles 
    #fdelXarr[0] = 0; 
    fdelXarr[-1] = 0
    #don't merge particles around special locations (TBD)

    imL, = np.nonzero(xdel<fdelXarr*fdelM)


    fdelXarr[imL] = np.inf #dont split where we merge
    fdelXarr[imL+1] = np.inf
    fdelXarr[imL-1] = np.inf
    fdelXarr[0] = np.inf; fdelXarr[-1] = np.inf #don't split 1st/last particles
    #don't split particles around special locations (TBD)

    #fdelS = 2*sim.particles.delta
    isL, = np.nonzero(xdel>fdelS*fdelXarr)

    if len(isL)>0 or len(imL)>0:
        doResample = True
    else:
        doResample = False

    if len(imL)>0:
        print('Merge happens')

    if doResample :#or len(imL)>0:
        print('Resample happens') 

        mtot = spN.mtotL
        mphy = spN.massL
        fcomp = spN.fcomp #composition fraction
        ncomp = len(fcomp[0])

        #splitting: add the locations
        addlocS = np.sqrt(loc[isL]*loc[isL+1])
        addlocM = np.sqrt(loc[imL]*loc[imL+1])

        locmidext = locmid_ext (loc)
        cummtot = np.concatenate(([0],np.cumsum(mtot)))

        #merging: remove the locations from loc (TBD)
        locn = loc.copy()
        locn[imL] = addlocM
        locn = np.delete(locn,imL+1)

        #find the index the mphy should be inserted
        idxS = np.searchsorted(locn, addlocS)

        #a bit weird
        locn = np.concatenate((locn,addlocS))
        locn.sort()


        npar = len(locn) #new number of particles

        locmidnext = locmid_ext (locn)
        locmidnext[0] = locmidext[0] #hack
        cummtotn = np.interp(locmidnext, locmidext, cummtot)
        mtotn = np.diff(cummtotn)

        #print(isL)
        mflux = sim.particles.v_r *sim.particles.sfd *sim.particles.locL
        if np.any(np.diff(mflux[1:20])<0) and False:
            print('mflux not in order')
            import pdb; pdb.set_trace()

        ##[25.01.20]lzx: now from the velocity to get the mass 
        #first get the velocity in the geometrical average point 
        addmphyS = get_addmphy(isL, addlocS, sim, loc, prop='rv') 
        addmphyM = get_addmphy(imL, addlocM, sim, loc, prop='rv')

        #insert the phycical mass into massL 
        mphyn = mphy.copy()
        mphyn[imL] = addmphyM 
        mphyn = np.delete(mphyn, imL+1)

        mphyn = np.insert(mphyn, idxS, addmphyS) 

        #if len(imL)>0:
        #    import userfun
        #    userfun.ba_resample(loc, locn, mphy, mphyn, mtot, mtotn, isL, imL,sim.time)

        #logmphyn = interp_mtot_weighted (locmidnext, locmidext, np.log(mphy))
        #mphyn = np.exp(logmphyn)
        #import pdb; pdb.set_trace()

        #mphyn[isL+1] = (0.5*(mphy[isL]**(1/3)+mphy[isL+2]**(1/3)))**3

        #if np.all(np.diff(np.log10(mphyn[:100]))<0)==False:
        #    print('physical mass not in order')

        #sfdnew = ff.sfd_special (mtotn, locn, sim.specloc)

        #[25.06.05]: don;t know how to treat softening length yet, now just crudely treat here. 
        if pars.sfdmode == 'sfd_kernel':
            hsoft = sim.particles.hsoft 
            hsoftn = hsoft.copy()
            hsoftn[imL] = np.sqrt(hsoft[imL]*hsoft[imL+1]) 
            hsoftn = np.delete(hsoftn, imL+1) 
            hsoftn = np.insert(hsoftn, idxS, np.sqrt(hsoft[isL]*hsoft[isL+1]))

        #composition... to be tested
        fcompn = np.empty((npar,ncomp))
        for k in range(ncomp):
            cummass = np.concatenate(([0], np.cumsum(mtot*fcomp[:,k])))
            cummassn = np.interp(locmidnext, locmidext, cummass)
            fcompn[:,k] = np.diff(cummassn) /mtotn

        #to make sure there are only once resample 
        #import userfun
        #userfun.ba_resample(loc, locn, mphy, mphyn, mtot, mtotn, isL, imL, sim.time)


        if full_output: 
            return locn, mtotn, mphyn, fcompn, isL, imL
        elif pars.sfdmode == 'sfd_kernel':
            return locn, mtotn, mphyn, fcompn, hsoftn
        else:
            return locn, mtotn, mphyn, fcompn
    else:
        return None





def re_sample_splitmerge (sim, spN, fdelS, fdelM=0., fdelX=1, nsampleX=0, fdelDM=0.0001, full_output=False):
    """
    [22.04.10]: mmid picked through midpoints
    [22.07.28]: also added the "simple merge" (splitsimplemerge) algorithm
                in which only the 2 core particles are affected. But this is
                not recommended!
                
    resample the distribution by (binary) merging and splitting

    the fractional spacing "xdel" b/w sp is calculated
    - when xdel>(split criterion) a new sp is erected at the midpoint
    - when xdel<(merge criterion) the two sp are combined at the midpoint

    Future:
    - add other, possibly user-defined, criteria for merging

    sp pos:         0   1   2   3   ...
    midpoints:        0   1   2

    Parameters: 
    
    fdelS: the fractionally critical distance of particles for splitting
    fdelM: the fractionally critical distance of particles for merging 
    fdelDM: the fractionally critical distance of particles for directly merging
    """
    loc = spN.locL 
    marr = spN.mtotL

    #determine no. particles involved in merging on each side
    if pars.resampleMode=='splitsimplemerge':
        nmerge = 1
    else:
        nmerge = 2


    #gather the special locations
    locspecL = []
    for k, line in enumerate(sim.icelineL+sim.planetL):
        locspecL.append(line.loc)

    #corresponding special (midpoint) indices
    #b/c midpoints we subtract one
    ixS = np.searchsorted(loc, locspecL)
    imS = ixS -1

    #logarithmic difference
    #Problem! 
    #The lcoation of particles may not be in order
    xdel = np.diff(np.log(loc))

    #take a finer sampling interval arround special locations
    fdelXarr = np.ones_like(xdel)
    for k in ixS:
        fdelXarr[k:k+nsampleX] = fdelX

    #special locations OTOH may never be resampled
    #[22.04.17] account for when planet is exterior to outer particle 
    ii = imS<len(fdelXarr) 
    fdelXarr[imS[ii]] = np.inf

    #midpoint locations elligible for splitting (isL)
    isL, = np.nonzero(xdel>fdelS *fdelXarr)

    #[22.10.25] determine the merging locations (imL)
    #freeze possible mergers of particles approaching special
    for k in ixS:
        fdelXarr[k:k+2*nsampleX] = fdelX

    #midpoint location for merging...
    #avoid neighboring indices and those in isL above (!)
    #I dont know yet how to avoid a for loop...
    idumL, = np.where(xdel<fdelM *fdelXarr)
    if len(idumL)>0:
        imL = []
        for idum in idumL:
            iarr2 = np.arange(idum-nmerge+1,idum+nmerge)
            iarr3 = np.arange(idum-nmerge,idum+nmerge+1)
            con0o = idum<len(marr)-nmerge   #outer boundaries             
            con0i = idum>=nmerge-1          #[22.07.28] inner boundary
            con1 = len(imL)==0 or idum>=imL[-1]+2*nmerge    #[22.04.16]distance 
            con2 = not(set(iarr2) & set(imS)) #much faster than instersect1d
            con3 = not(set(iarr3) & set(isL)) #much faster than instersect1d
            if con0o and con0i and con1 and con2 and con3:
                imL.append(idum)
    else:
        imL = []
    imL = np.array(imL,dtype=int)


    #[22.10.25] move this one up...
    #only set up a new structure when splitting or merging has actually  happend
    if len(imL)>0 or len(isL)>0:

        #[24.05.26]:add to message list
        #           in future ix should point to hash
        for ix in imL:
            sim.add_message('merge', f'particle {ix} has been selected for merging')
        for ix in isL:
            sim.add_message('split', f'particle {ix} has been selected for splitting')
        

        #the locations where new particles are made (geometric mean)
        loc_split  = np.sqrt(loc[isL]*loc[isL+1])


        #wdel is a bit like xdel --- but gives support 
        #width of particle in symmetric way
        wdel = np.concatenate(([np.log(loc[1]/loc[0])], 
                               (np.log(loc[2:]) -np.log(loc[:-2]))/2,
                               [np.log(loc[-1]/loc[-2])]))


        xdelNl = (np.log(loc_split)  -np.log(loc[isL]) ) /2 
        xdelNu = (np.log(loc[isL+1]) -np.log(loc_split)) /2

        #[22.04.12]
        #fraction of mass donated by upper/lower particles
        #and remaining mass fraction from existing particles
        fnewu = xdelNu /wdel[isL+1]
        fnewl = xdelNl /wdel[isL]
        frem = np.ones_like(loc)
        frem[isL]   -= fnewl 
        frem[isL+1] -= fnewu

        msup_split = fnewl*spN.mtotL[isL]   +fnewu*spN.mtotL[isL+1]
        mphy_split = (fnewl*spN.massL[isL]  +fnewu*spN.massL[isL+1]) /(fnewl+fnewu)   #weighed

        #[22.03.29]changed since order has changed
        fcom_split = fnewl[:,np.newaxis]*spN.fcomp[isL] +fnewu[:,np.newaxis]*spN.fcomp[isL+1] 
        fcom_split /= (fnewl[:,np.newaxis] +fnewu[:,np.newaxis])

        #... continue w/ mering
        loc_merge = np.sqrt(spN.locL[imL]*spN.locL[imL+1])

        #there cannot be overlap b/w merged and splitting particles
        if np.any(frem[imL]!=1.) or np.any(frem[imL+1]!=1.):
            print('ERROR!')

        #merging particles may also donate mass to their neighbors
        #[22.04.12]: this doesn't work at alll... probably a bug somewhere
        #[22.04.16]: solved. Need to care about changing particles' physical properties

        #algorithm (ii) -- donate according to surface density of the merged particles
        #[22.07.28]: doesnt happen in splitsimplemerge algorithm
        if pars.resampleMode=='splitmerge':
            xdelNl = (np.log(loc_merge)  -np.log(loc[imL-1])) /2 
            xdelNr = (np.log(loc[imL+2]) -np.log(loc_merge) ) /2
            delm1 = (wdel[imL]   -xdelNl) /wdel[imL]   *marr[imL]   #mass transfer to (i-1)
            delm2 = (wdel[imL+1] -xdelNr) /wdel[imL+1] *marr[imL+1] #mass transfer to (i+2)
        else:
            delm1 = 0
            delm2 = 0

        #mass conservation: add i+i+1, subtract donors
        mass_merge = spN.mtotL[imL] +spN.mtotL[imL+1] -delm1 -delm2

        if pars.resampleMode=='splitmerge':
            mass_merge1 = spN.mtotL[imL-1] +delm1 #particle i-1
            mass_merge2 = spN.mtotL[imL+2] +delm2 #particle i+2

        try:
            assert(np.all(mass_merge>0))
        except:
            print('[resample.splitmerge]:ERROR: negative mass after merging!')
            import pdb; pdb.set_trace()
            sys.exit()

        #merging locations will be removed
        frem[imL]   = 0.0
        frem[imL+1] = 0.0
        if pars.resampleMode=='splitmerge':
            frem[imL-1] = 0.0
            frem[imL+2] = 0.0


        #simple merging. Doesn't work well as sfd no longer smooth
        #frem[imL] = 0.0
        #frem[imL+1] = 0.0
        #mass_merge = spN.mtotL[imL] +spN.mtotL[imL+1]


        #properties of the merged particles and its new neighbors
        #this is very important!

        mphy_merge = ( spN.massL[imL] * (marr[imL]-delm1) +spN.massL[imL+1]  *(marr[imL+1]-delm2)) /mass_merge
        fcom_merge = ((spN.fcomp[imL].T*(marr[imL]-delm1) +spN.fcomp[imL+1].T*(marr[imL+1]-delm2)) /mass_merge).T

        if pars.resampleMode=='splitmerge':
            mphy_merge1 = ( spN.massL[imL-1]*marr[imL-1] +spN.massL[imL]  *delm1) /mass_merge1
            mphy_merge2 = ( spN.massL[imL+2]*marr[imL+2] +spN.massL[imL+1]*delm2) /mass_merge2

            fcom_merge1 = ((spN.fcomp[imL-1].T*marr[imL-1] +spN.fcomp[imL].T*delm1) /mass_merge1).T
            fcom_merge2 = ((spN.fcomp[imL+2].T*marr[imL+2] +spN.fcomp[imL+1].T*delm2) /mass_merge2).T

            loc_merge1 = spN.locL[imL-1]
            loc_merge2 = spN.locL[imL+2]

   ##TBR -- old branching location
   ##only set up a new structure when splitting or merging has actually  happend
   #if len(imL)>0 or len(isL)>0:

        if pars.resampleMode=='splitmerge':
            dumiiL = [isL+1,imL,imL+1,imL+2]
            dumloc = [loc_split,loc_merge1,loc_merge,loc_merge2]
            dumsup = [msup_split,mass_merge1, mass_merge, mass_merge2]
            dumphy = [mphy_split, mphy_merge1, mphy_merge, mphy_merge2]
            dumcom = [fcom_split, fcom_merge1, fcom_merge, fcom_merge2]
        else:
            dumiiL = [isL+1,imL]
            dumloc = [loc_split, loc_merge]
            dumsup = [msup_split, mass_merge]
            dumphy = [mphy_split, mphy_merge]
            dumcom = [fcom_split, fcom_merge]

        #new sp arrays
        locI =  np.insert(loc,      np.concatenate(dumiiL), np.concatenate(dumloc))
        msupI = np.insert(marr*frem,np.concatenate(dumiiL), np.concatenate(dumsup))
        mphyI = np.insert(spN.massL, np.concatenate(dumiiL), np.concatenate(dumphy))
        fcomI = np.insert(spN.fcomp,np.concatenate(dumiiL), np.concatenate(dumcom), axis=0)


        #when sp merge, remove the old locations
        iremL, = np.nonzero(msupI==0)
        if len(iremL)>0:
            locI = np.delete(locI, iremL)
            msupI = np.delete(msupI, iremL)
            mphyI = np.delete(mphyI, iremL)
            fcomI = np.delete(fcomI, iremL, axis=0)
            #print('removed', len(iremL)//4, ' particles')

        #[24.05.15]LZX: because the location of particles is not in order at early stage, 
        #so comment this out for now.
        #try:
        #    assert(np.all(np.diff(locI)>0)) #locations must be increasing
        #except:
        #    print('[resample.splitmerge]:locations not increasing')
        #    import pdb; pdb.set_trace()

        try:
            assert(np.abs(msupI.sum()/marr.sum()-1) <1e-10) #mass conservations
        except:
            print('[resample.splitmerge]:mass conservation violated')
            import pdb; pdb.set_trace()

        if full_output:
            newarr = (locI, msupI, mphyI, fcomI, isL, imL)
        else:
            newarr = (locI, msupI, mphyI, fcomI)


    #neither merging or splitting has occurred
    #[24.05.15]cwo: omit DM for the moment
    else:
        #search of "Direct Merging" (2 particles very close)
        #left index loc that potentially merge with right neigbor
        #(first 2, 2 right of special, 2 left of special, final 2)
        ixL = [0] +ixS.tolist() +(imS-1).tolist() + [len(loc)-2]   

        #remove the special locations (we cannot merge over them!)
        #[22.12.15]:also remove the final index in case any ixS is only interior to last one
        ixL = list(set(ixL) -set(imS) -set([-1]) -set([len(loc)-1, len(loc)]))

        #check if any meet the direct merging condition
        con_direct = xdel[ixL]<fdelDM

        if con_direct.any():
            #only 1 direct merge per timestep!
            ic, = con_direct.nonzero()
            try:
                im = ixL[int(ic)] #the merge index
            except:
                import pdb;pdb.set_trace()

            print('[resample.splitmerge]:direct merge particle no.', im)

            loc_merge = np.sqrt(spN.locL[im]*spN.locL[im+1])
            mass_merge = spN.mtotL[im] +spN.mtotL[im+1]
            mphy_merge =  (spN.massL[im]   *marr[im] +spN.massL[im+1]   *marr[im+1]) /mass_merge
            fcom_merge = ((spN.fcomp[im].T*marr[im] +spN.fcomp[im+1].T*marr[im+1]) /mass_merge).T

            locI = np.concatenate ((spN.locL[:im],   [loc_merge],    spN.locL[im+2:]))
            msupI = np.concatenate((marr[:im],      [mass_merge],   marr[im+2:]))
            mphyI = np.concatenate((spN.massL[:im],  [mphy_merge],   spN.massL[im+2:]))
            fcomI = np.concatenate((spN.fcomp[:im], [fcom_merge],   spN.fcomp[im+2:]))
            newarr = (locI, msupI, mphyI, fcomI)

        else:
            newarr = None

    #if newarr is not None:
        #import pdb;pdb.set_trace()
    #if len(imL>0):
    #    import pdb;pdb.set_trace()
    return newarr

def re_sample_dropmerge(sim, spN, fdelS=0.05, fdelM=0., fdelX=1, nsampleX=0, fdelDM = 0.0):
    """
    new algorithm: drop: 
    When the distance b/w too particles become too large, 
    the inner particle will drop some mass in the logarithmical 
    average point b/w them. 

    Drop: 
    delDmass = 1/2 * mtot(i-1)
    """
    loc = spN.locL 
    marr = spN.mtotL

    #determine no. particles involved in merging on each side
    if pars.resampleMode=='dropsimplemerge':
        nmerge = 1
    else:
        nmerge = 2


    #gather the special locations
    locspecL = []
    for k, line in enumerate(sim.icelineL+sim.planetL):
        locspecL.append(line.loc)

    #corresponding special (midpoint) indices
    #b/c midpoints we subtract one
    ixD = np.searchsorted(loc, locspecL)
    imD = ixD -1

    #logarithmic difference
    #Problem! 
    #The lcoation of particles may not be in order
    xdel = np.diff(np.log(loc))

    #take a finer sampling interval arround special locations
    fdelXarr = np.ones_like(xdel)
    for k in ixD:
        fdelXarr[k:k+nsampleX] = fdelS

    #special locations OTOH may never be resampled
    #[22.04.17] account for when planet is exterior to outer particle 
    ii = imD<len(fdelXarr) 
    fdelXarr[imD[ii]] = np.inf

    #midpoint locations elligible for splitting (isL)
    idL, = np.nonzero(xdel>fdelS *fdelXarr)

    #[22.10.25] determine the merging locations (imL)
    #freeze possible mergers of particles approaching special
    for k in ixD:
        fdelXarr[k:k+2*nsampleX] = fdelX

    #midpoint location for merging...
    #avoid neighboring indices and those in isL above (!)
    #I dont know yet how to avoid a for loop...
    idumL, = np.where(xdel<fdelM *fdelXarr)
    if len(idumL)>0:
        imL = []
        for idum in idumL:
            iarr2 = np.arange(idum-nmerge+1,idum+nmerge)
            iarr3 = np.arange(idum-nmerge,idum+nmerge+1)
            con0o = idum<len(marr)-nmerge   #outer boundaries             
            con0i = idum>=nmerge-1          #[22.07.28] inner boundary
            con1 = len(imL)==0 or idum>=imL[-1]+2*nmerge    #[22.04.16]distance 
            con2 = not(set(iarr2) & set(imD)) #much faster than instersect1d
            con3 = not(set(iarr3) & set(idL)) #much faster than instersect1d
            if con0o and con0i and con1 and con2 and con3:
                imL.append(idum)
    else:
        imL = []
    imL = np.array(imL,dtype=int)


    #[22.10.25] move this one up...
    #only set up a new structure when splitting or merging has actually  happend
    if len(imL)>0 or len(idL)>0:

        #[24.05.26]:add to message list
        #           in future ix should point to hash
        for ix in imL:
            sim.add_message('merge', f'particle {ix} has been selected for merging')
        for ix in idL:
            sim.add_message('drop', f'particle {ix} has been selected for dropping')
        

        #the locations where new particles are made (geometric mean)
        loc_drop  = np.sqrt(loc[idL]*loc[idL+1])


        #wdel is a bit like xdel --- but gives support 
        #width of particle in symmetric way
        wdel = np.concatenate(([np.log(loc[1]/loc[0])], 
                               (np.log(loc[2:]) -np.log(loc[:-2]))/2,
                               [np.log(loc[-1]/loc[-2])]))


        xdelNl = (np.log(loc_drop)  -np.log(loc[idL]) ) /2 
        xdelNu = (np.log(loc[idL+1]) -np.log(loc_drop)) /2

        #[22.04.12]
        #fraction of mass donated by upper/lower particles
        #and remaining mass fraction from existing particles
        #[24.10.14]Artificially define the mass dropped from the inner particles is half of its mass 
        fnewu = 0.0*np.ones(len(idL))
        fnewl = 1/2*np.ones(len(idL))
        frem = np.ones_like(loc)
        frem[idL]   -= fnewl 
        frem[idL+1] -= fnewu

        msup_drop = fnewl*spN.mtotL[idL]   +fnewu*spN.mtotL[idL+1]
        #[24.10.14]keep the physical mass being weighted from upper and lower mass 
        mphy_drop = (fnewl*spN.massL[idL]  +fnewu*spN.massL[idL+1]) /(fnewl+fnewu)   #weighed

        #[22.03.29]changed since order has changed
        fcom_drop = fnewl[:,np.newaxis]*spN.fcomp[idL] +fnewu[:,np.newaxis]*spN.fcomp[idL+1] 
        fcom_drop /= (fnewl[:,np.newaxis] +fnewu[:,np.newaxis])

        #... continue w/ mering
        loc_merge = np.sqrt(spN.locL[imL]*spN.locL[imL+1])

        #there cannot be overlap b/w merged and splitting particles
        if np.any(frem[imL]!=1.) or np.any(frem[imL+1]!=1.):
            print('ERROR!')

        #merging particles may also donate mass to their neighbors
        #[22.04.12]: this doesn't work at alll... probably a bug somewhere
        #[22.04.16]: solved. Need to care about changing particles' physical properties

        #algorithm (ii) -- donate according to surface density of the merged particles
        #[22.07.28]: doesnt happen in splitsimplemerge algorithm
        if pars.resampleMode=='dropmerge':
            xdelNl = (np.log(loc_merge)  -np.log(loc[imL-1])) /2 
            xdelNr = (np.log(loc[imL+2]) -np.log(loc_merge) ) /2
            delm1 = (wdel[imL]   -xdelNl) /wdel[imL]   *marr[imL]   #mass transfer to (i-1)
            delm2 = (wdel[imL+1] -xdelNr) /wdel[imL+1] *marr[imL+1] #mass transfer to (i+2)
        else:
            delm1 = 0
            delm2 = 0

        #mass conservation: add i+i+1, subtract donors
        mass_merge = spN.mtotL[imL] +spN.mtotL[imL+1] -delm1 -delm2

        if pars.resampleMode=='dropmerge':
            mass_merge1 = spN.mtotL[imL-1] +delm1 #particle i-1
            mass_merge2 = spN.mtotL[imL+2] +delm2 #particle i+2

        try:
            assert(np.all(mass_merge>0))
        except:
            print('[resample.splitmerge]:ERROR: negative mass after merging!')
            import pdb; pdb.set_trace()
            sys.exit()

        #merging locations will be removed
        frem[imL]   = 0.0
        frem[imL+1] = 0.0
        if pars.resampleMode=='dropmerge':
            frem[imL-1] = 0.0
            frem[imL+2] = 0.0


        #simple merging. Doesn't work well as sfd no longer smooth
        #frem[imL] = 0.0
        #frem[imL+1] = 0.0
        #mass_merge = spN.mtotL[imL] +spN.mtotL[imL+1]


        #properties of the merged particles and its new neighbors
        #this is very important!

        mphy_merge = ( spN.massL[imL] * (marr[imL]-delm1) +spN.massL[imL+1]  *(marr[imL+1]-delm2)) /mass_merge
        fcom_merge = ((spN.fcomp[imL].T*(marr[imL]-delm1) +spN.fcomp[imL+1].T*(marr[imL+1]-delm2)) /mass_merge).T

        if pars.resampleMode=='dropmerge':
            mphy_merge1 = ( spN.massL[imL-1]*marr[imL-1] +spN.massL[imL]  *delm1) /mass_merge1
            mphy_merge2 = ( spN.massL[imL+2]*marr[imL+2] +spN.massL[imL+1]*delm2) /mass_merge2

            fcom_merge1 = ((spN.fcomp[imL-1].T*marr[imL-1] +spN.fcomp[imL].T*delm1) /mass_merge1).T
            fcom_merge2 = ((spN.fcomp[imL+2].T*marr[imL+2] +spN.fcomp[imL+1].T*delm2) /mass_merge2).T

            loc_merge1 = spN.locL[imL-1]
            loc_merge2 = spN.locL[imL+2]

   ##TBR -- old branching location
   ##only set up a new structure when splitting or merging has actually  happend
   #if len(imL)>0 or len(isL)>0:

        if pars.resampleMode=='dropmerge':
            dumiiL = [idL+1,imL,imL+1,imL+2]
            dumloc = [loc_drop,loc_merge1,loc_merge,loc_merge2]
            dumsup = [msup_drop,mass_merge1, mass_merge, mass_merge2]
            dumphy = [mphy_drop, mphy_merge1, mphy_merge, mphy_merge2]
            dumcom = [fcom_drop, fcom_merge1, fcom_merge, fcom_merge2]
        else:
            dumiiL = [idL+1,imL]
            dumloc = [loc_drop, loc_merge]
            dumsup = [msup_drop, mass_merge]
            dumphy = [mphy_drop, mphy_merge]
            dumcom = [fcom_drop, fcom_merge]

        #new sp arrays
        locI =  np.insert(loc,      np.concatenate(dumiiL), np.concatenate(dumloc))
        msupI = np.insert(marr*frem,np.concatenate(dumiiL), np.concatenate(dumsup))
        mphyI = np.insert(spN.massL, np.concatenate(dumiiL), np.concatenate(dumphy))
        fcomI = np.insert(spN.fcomp,np.concatenate(dumiiL), np.concatenate(dumcom), axis=0)


        #when sp merge, remove the old locations
        iremL, = np.nonzero(msupI==0)
        if len(iremL)>0:
            locI = np.delete(locI, iremL)
            msupI = np.delete(msupI, iremL)
            mphyI = np.delete(mphyI, iremL)
            fcomI = np.delete(fcomI, iremL, axis=0)
            #print('removed', len(iremL)//4, ' particles')

        #[24.05.15]LZX: because the location of particles is not in order at early stage, 
        #so comment this out for now.
        #try:
        #    assert(np.all(np.diff(locI)>0)) #locations must be increasing
        #except:
        #    print('[resample.splitmerge]:locations not increasing')
        #    import pdb; pdb.set_trace()

        try:
            assert(np.abs(msupI.sum()/marr.sum()-1) <1e-10) #mass conservations
        except:
            print('[resample.dropmerge]:mass conservation violated')
            import pdb; pdb.set_trace()
        newarr = (locI, msupI, mphyI, fcomI)


    #neither merging or splitting has occurred
    #[24.05.15]cwo: omit DM for the moment
    else:
        #search of "Direct Merging" (2 particles very close)
        #left index loc that potentially merge with right neigbor
        #(first 2, 2 right of special, 2 left of special, final 2)
        ixL = [0] +ixD.tolist() +(imD-1).tolist() + [len(loc)-2]   

        #remove the special locations (we cannot merge over them!)
        #[22.12.15]:also remove the final index in case any ixS is only interior to last one
        ixL = list(set(ixL) -set(imD) -set([-1]) -set([len(loc)-1, len(loc)]))

        #check if any meet the direct merging condition
        #[24.10.14]: dont consider direct merge now 
        con_direct = xdel[ixL] < fdelDM

        if con_direct.any():
            #only 1 direct merge per timestep!
            ic, = con_direct.nonzero()
            try:
                im = ixL[int(ic)] #the merge index
            except:
                import pdb;pdb.set_trace()

            print('[resample.dropmerge]:direct merge particle no.', im)

            loc_merge = np.sqrt(spN.locL[im]*spN.locL[im+1])
            mass_merge = spN.mtotL[im] +spN.mtotL[im+1]
            mphy_merge =  (spN.massL[im]   *marr[im] +spN.massL[im+1]   *marr[im+1]) /mass_merge
            fcom_merge = ((spN.fcomp[im].T*marr[im] +spN.fcomp[im+1].T*marr[im+1]) /mass_merge).T

            locI = np.concatenate ((spN.locL[:im],   [loc_merge],    spN.locL[im+2:]))
            msupI = np.concatenate((marr[:im],      [mass_merge],   marr[im+2:]))
            mphyI = np.concatenate((spN.massL[:im],  [mphy_merge],   spN.massL[im+2:]))
            fcomI = np.concatenate((spN.fcomp[:im], [fcom_merge],   spN.fcomp[im+2:]))
            newarr = (locI, msupI, mphyI, fcomI)

        else:
            newarr = None

    #if newarr is not None:
        #import pdb;pdb.set_trace()
    #if len(imL>0):
    #    import pdb;pdb.set_trace()
    return newarr


def locmid_ext (loc):
    """
    provides the midpoints and extends them 
    """
    locmid = np.sqrt(loc[1:]*loc[:-1])
    locmidext = np.concatenate((
        [loc[0]*np.sqrt(loc[0]/loc[1])],
        locmid,
        [np.sqrt(loc[-1]/loc[-2])*loc[-1]]))
    #[24.0112]LZX: here we use this thing to be consistent with the sfd_simple

    return locmidext


def locmid_ext_special (loc, specL=[]):
    """
    this is a modified version of the above, in which
    midpoint locations around special locations are shifted
    to concide with the latter
    """
    locmidext = locmid_ext (loc)

    for locs in specL:
        il = np.searchsorted(loc, locs)
        ix = np.searchsorted(locmidext, locs)
        if loc[il]>loc[ix]:# p<special<ix<il
            locmidext[ix] = locs
        else:# ix-1<special<il<ix
            locmidext[ix-1] = locs

    return locmidext


def interp_mtot_weighted (xn, xmid, qarr, marr=None, mn=None, 
                          neval=0, sigval=1e-9):
    """
    interpolates a certain quantity qarr defined with respect to the 
    midpoints xmid onto new locations (midpoints) xn, 
    weighted (optionally) by masses mn

    [25.01.18]:try to detect loss of significance
    """
    if marr is None: marr = np.ones_like(qarr)
    if mn is None: mn = np.ones(len(xn)-1)

    neval += 1
    if neval==1:
        ix = 0
    else:
        ix = (xmid>xn[0]).argmax() -1

    cummass = np.concatenate(([0], np.cumsum(marr[ix:]*qarr[ix:])))
    cummassn = np.interp(xn, xmid[ix:], cummass)
    diffcumm = np.diff(cummassn)
    qn = diffcumm /(mn+1e-16) #1e-16 to prevent the zero division 

    #relative error/significance
    sig = diffcumm/cummassn[1:]
    itrust = (sig<sigval).argmax() #trust until here

    if itrust>0:
        qn[itrust:] = interp_mtot_weighted (xn[itrust:], xmid, qarr, marr, mn[itrust:], neval, sigval)

    return qn


def fixed_resample (sim, spN, specloc, fchange=0.9, Xspecial=1, **kwargs):
    """
    like global_resample, resample to fixed positions
    - finer resampling near iceline locations (TBD)
    - piecewise, but not mass-conserving

    [25.01.23]: accounting for "specials" makes the algorith very terse bookkeeping
    """

    #we need to make copies b/c of some b/c this could change...
    loc = spN.locL
    mtot = spN.mtotL.copy()
    mphy = spN.massL.copy()
    fcomp = spN.fcomp.copy() #composition fraction

    ncomp = len(fcomp[0])

    xdel = np.diff(np.log(loc))

    #splitting/changes based on fdelS
    fdelS = sim.particles.delta /fchange

    #initial locations and midpoins
    locn, locmidnext = sim.particles.loc_init(specL=sim.specloc)
    npar = len(locn)

    #areas surrounding special locations are better resolved
    #and need to have a more stringent merging criterion
    fdelM = np.ones_like(xdel) *sim.particles.delta *fchange
    conspecial = False
    for spec in list(sim.specloc):
        ii = (loc>spec*np.exp(-2*sim.particles.delta)) *\
                (loc<spec*np.exp(2*sim.particles.delta))
        fdelM[ii[:-1]] /= Xspecial #contains iceline


    conmerge = xdel<fdelM
    consplit = np.any(xdel>fdelS)

    nspec = len(specloc)
    specL = [0] +list(specloc)
    loc1 = sim.rout
    loc1mod = sim.rout #modified outer bounds

    mloss = 0
    if conmerge.any() or consplit:

        mphyn = np.zeros_like(locn)
        mtotn = np.zeros_like(locn)
        fcompn = np.zeros((npar,ncomp))

        #segmented resample, inverse order!
        #This is very annoying bookkeeping
        for kseg,loc0 in enumerate(specL[::-1]):
            ii = (loc0<loc) *(loc<loc1mod)
            iin = (loc0<locn) *(locn<=loc1)

            #make cumulative array
            locmidext = locmid_ext (loc[ii])
            cummtot = np.concatenate(([0],np.cumsum(mtot[ii])))

            #corresponding midpoint indices that span iin
            imn0 = iin.argmax()
            ss = slice(imn0,imn0 +np.sum(iin)+1)

            if locmidnext[imn0]<locmidext[0]:
                #I addded this but I don't think it's necessary...
                #locmidext[0] = locmidnext[imn0]
                problem = False
            else:
                #Here there is overflow into the next segment!
                #it would be very BAD to adjust locmidext, which results in mass
                #pileup near the special location
                problem = True


            cummd = np.interp(locmidnext[ss], locmidext, cummtot)

            #this is the mass lost at the inner boundary...
            if kseg==nspec:
                mloss += cummd[0]

            #sometimes midpoints exceeds domian of midpoint-new..
            #.. need to ensure we include all mass at the outer domain end
            cummd[-1] = cummtot[-1]
            mtotn[iin] = np.diff(cummd)

            mphyn[iin] = interp_mtot_weighted (locmidnext[ss], locmidext, mphy[ii], mtot[ii], mtotn[iin])
            for k in range(ncomp):
                fcompn[iin,k] = interp_mtot_weighted (locmidnext[ss], locmidext, fcomp[ii,k], mtot[ii], mtotn[iin])

            #check the mphy around the iceline  
            # import cgs
            # plt.loglog(loc, mphy/mphy.max(), 'x-', label ='old')
            # plt.loglog(locn, mphyn/mphy.max(), 'x-', label = str(sim.time/cgs.yr))
            # plt.loglog(loc, mtot/mtot.max(), '.-', label='oldmtot')
            # plt.loglog(locn, mtotn/mtot.max(), '.-', label='newmtot')
            # plt.axvline(sim.specloc[0], ls='--', lw=1, c = 'gray')
            # plt.xlim(sim.specloc[0]*0.9, sim.specloc[0]*1.1)
            # plt.legend()
            # plt.savefig('mphy/{:.2f}.jpg'.format(sim.time))
            # plt.close()
            #The Problem particle -- that would be created at the location interior to the special
            #This represents mass overflow.
            if imn0>0 and problem: 
                mtotp = cummd[0] #the mass not included in the above

                #other properties are taken from the "old" particle
                ii0 = ii.argmax()
                spiP = sim.particles.select_single(ii0)
                spiP.mtotL = mtotp

                #this particle should cross the special location (e.g., lose its ice; pebble accretion):
                #   spiP --> spiP'
                #and from spiP' we extract the new properties and we replace the particle
                #at ii0 with these new properties (and mtot->mtotp)

                #so stuff TBD here!!

                #here we mimic water icelines where 50% is lost
                #for illustrative purposes:
                fcomp[ii0] = np.array([1,0]) #it has crossed the iceline!
                mphy[ii0] /=2       #physical mass is reduces
                mtot[ii0] = mtotp/2 #only a fraction has crossed
                mloss += mtotp/2

                loc1mod = (1+1e-8)*loc[ii0] #make sure it will be included in the next segment
            else:
                loc1mod = loc0

            loc1 = loc0

        #check mass conservation...
        err = (mtotn.sum() +mloss)/spN.mtotL.sum()-1
        if abs(err)>2e-14:
            print('[fixed_resample]:mass loss detected')
            import pdb; pdb.set_trace()

        return locn, mtotn, mphyn, fcompn
    else:
        return None



def global_resample4 (sim, spN, fchange=0.5, fdelX=1, nsampleX=0, nn=1,**args):
    """
    similar to global_resample2... but with:
    -- segmented sampling (special locations; specL)
    -- do not adjust positions of nn particles near specL
    """
    n1 = nn-1 #so 0 or 1

    loc = spN.locL 
    mtot = spN.mtotL #total mass
    mphy = spN.massL #physical mass
    fcomp = spN.fcomp #composition fraction
    ncomp = len(fcomp[0])

    #general
    fdelS = sim.particles.delta/fchange
    fdelM = sim.particles.delta*fchange

    specL = list(sim.specloc) +[np.inf]

    #first consider need of resample per segment
    doResample = []
    loc0 = 0
    for kseg, loc1 in enumerate(specL):

        ii = (loc0<loc) *(loc<loc1)

        #we only interested in the spacing among these particles in the segment
        xdel = np.diff(np.log(loc[ii]))


        #TBD: direct merge and single split

        #normal mode
        #the first/last particles two cannot be merged...
        #so need to have at least 5 particles in the segment!
        ss = slice(n1,len(xdel)-n1)
        if (np.any(xdel[ss]<fdelM) or np.any(xdel[ss]>fdelS)) and sum(ii)>=2*nn+1:
            doResample.append(True)
        else:
            doResample.append(False)

        loc0 = loc1 #for next segment

    if np.any(doResample):
        loc0 = 0
        locnL=[]; mtotnL=[]; mphynL=[]; fcompnL=[] 
        for kseg, loc1 in enumerate(specL):
            ii = (loc0<loc) *(loc<loc1)
            if doResample[kseg]:

                #extended midpoints locations
                locmidext = locmid_ext (loc[ii])

                #cumulative mass function is defined on the midpoints
                cummtot = np.concatenate(([0],np.cumsum(mtot[ii])))

                if nn==2:
                    #number of particles to place in segment 
                    npar = int(np.log(loc[ii][-nn]/loc[ii][n1]) /sim.particles.delta)

                    #at these (new) locations
                    locn = np.concatenate(([loc[ii][0]],
                                np.exp(np.linspace(np.log(loc[ii][1]), np.log(loc[ii][-2]), npar+1)),
                                [loc[ii][-1]]))

                    npar = len(locn)

                elif nn==1:
                    npar = int(np.log(locmidext[-1]/locmidext[0]) /sim.particles.delta)
                    locmid = np.exp(np.linspace(np.log(locmidext[0]),np.log(locmidext[-1]),npar+1))
                    locn = np.sqrt(locmid[1:]*locmid[:-1])


                if False:
                    xdum = np.log(sim.rout/sim.rinn)
                    npar = int(xdum /sim.particles.delta)
                    rmid = sim.rinn *np.exp(np.linspace(0,1,npar+1)*xdum)
                    locn = np.sqrt(rmid[1:]*rmid[:-1])
                    #import pdb; pdb.set_trace()


                locmidnext = locmid_ext (locn)
                cummtotn = np.interp(locmidnext, locmidext, cummtot)
                mtotn = np.diff(cummtotn)

                ## this should give the same...
                #dum = interp_mtot_weighted (locmidnext, locmidext, mtot[ii])

                ## it would be very weird if the particles cross the boundary
                if nn==1 and kseg==0 and locn[-1]>loc1:
                    locn[-1] = (1-1e-10)*loc1

                if nn==1 and kseg==1 and locn[0]<loc0:
                    locn[0] = (1+1e-10)*loc0

                #with these rules we can also sample other (mass-weighted quantities)
                #if False:
                #cummass = np.concatenate(([0], np.cumsum(mtot[ii]*mphy[ii])))
                #cummassn = np.interp(locmidnext, locmidext, cummass)
                #dum = np.diff(cummassn) /(mtotn+1e-16) #1e-16 to prevent the zero division 

                #log-interpolation of the mass seems much better
                #logmphyn = interp_mtot_weighted (locmidnext, locmidext, np.log(mphy[ii]), mtot[ii], mtotn)
                #mphyn = np.exp(logmphyn)

                #this is very diffusive
                #mphyn = interp_mtot_weighted (locmidnext, locmidext, mphy[ii], mtot[ii], mtotn)

                #this is also quite diffusive
                pwl = 0.5
                dum = interp_mtot_weighted (locmidnext, locmidext, mphy[ii]**pwl, mtot[ii], mtotn)
                mphyn = dum**(1/pwl)

                #composition... to be tested
                fcompn = np.empty((npar,ncomp))
                for k in range(ncomp):
                    cummass = np.concatenate(([0], np.cumsum(mtot[ii]*fcomp[ii,k])))
                    cummassn = np.interp(locmidnext, locmidext, cummass)
                    fcompn[:,k] = np.diff(cummassn) /mtotn


                locnL.append(locn)
                mtotnL.append(mtotn)
                mphynL.append(mphyn)
                fcompnL.append(fcompn)

            else:#just copy stuff
                locnL.append(loc[ii])
                mtotnL.append(mtot[ii])
                mphynL.append(mphy[ii])
                fcompnL.append(fcomp[ii])

            loc0 = loc1 #for next segment

        locn = np.concatenate(locnL)
        mtotn = np.concatenate(mtotnL)
        mphyn = np.concatenate(mphynL)
        fcompn = np.concatenate(fcompnL)

        return locn, mtotn, mphyn, fcompn
    else:
        return None



def global_resample2 (sim, spN, fchange=0.9, fdelX=1, nsampleX=0, nspec=1,**args):
    """
    A variation on the below
    """

    loc = spN.locL 
    mtot = spN.mtotL #total mass
    mphy = spN.massL #physical mass
    fcomp = spN.fcomp #composition fraction
    ncomp = len(fcomp[0])

    xdel = np.diff(np.log(loc))
    npar = pars.dparticleprops['nini']

    #[24.12.30]: in global_resample, better to let fdelS, fdelM tied to the initial grid spacing
    #the disired xdel
    xdel_aim = np.log((sim.rout/sim.rinn)**(1.0/npar))
    fdelS = xdel_aim/fchange
    fdelM = xdel_aim*fchange
    fspec = (fdelS+fdelM)/2

    opL, = np.where(xdel<fdelM)
    opL = np.append(opL, np.where(xdel>fdelS)[0])

    if len(opL)>0:
        #new locations -- this follows initialization and we better save them!
        rdum = np.exp(np.linspace(np.log(sim.rinn),np.log(sim.rout),npar+1))


        locn = np.sqrt(rdum[1:]*rdum[:-1])

        #get the special locations 
        #[24/11/18] b/c we have consider the iceline location later, so here maybe we can consider the 
        #           planet location only.
        locspecL = [] 
        #try to also treat the inner edge as a special location 
        locspecL.append(sim.rinn)
        for line in sim.planetL+sim.icelineL:
            locspecL.append(line.loc)

        # not sure whether it's a good idea to add the outermost location of particles to the special location 
        # locspecL.append(loc[-1])

        #to prevent there are some overlapping special locations 
        locspecL = list(set(locspecL))


        #consider the special locations like I did before. Just first add the 1 surrounding particles 
                #find the locations that close to the special locations 
        specpar_idx = np.searchsorted(loc, locspecL) 
        #insert the surrounding particles into locn 
        #we also need to remove the locations that closer to the special locations in locn 

        #LZX[24/11/2] not sure whether we should get the 'special particles' within two particles or within a 
        #             certain location range.
        for i, idx in enumerate(specpar_idx):
            #maybe we should combine the two methods together. [24/11/27]
            #first get the special particles, and also special range, choose the minimal range to implement
            #inside of special range, we regard the particles as 'effective particles for special locations'
            if idx < nspec:
                specL_par = np.array(loc[:idx+nspec])
            elif idx>len(loc)-nspec:
                specL_par = np.array(loc[idx-nspec:])
            else:
                specL_par = np.array(loc[idx-nspec:idx+nspec]) 

            #'special range' 
            speclim = [locspecL[i]/np.exp(fspec), locspecL[i]*np.exp(fspec)]

            #'crop' the specL_par 
            specL_par = specL_par[np.where((specL_par>speclim[0]) & (specL_par<speclim[1]))[0]]


            #particles in special range 
            #specL_par = np.where((loc>speclim[0]) & (loc<speclim[1]))[0]

            if len(specL_par)>0:
                locn = np.delete(locn, np.where((locn>specL_par.min()) & (locn<specL_par.max())))
                #locn = np.delete(locn, np.where((locn>speclim[0]) & (locn<speclim[1])))
                locn = np.insert(locn, np.searchsorted(locn, specL_par), specL_par)
                #check: if the difference of the locations are too small, 
                #       just remove! 
                locn = np.delete(locn, np.where(np.diff(np.log(locn))<fdelM))



        #midpoints and extensions
        locmid = np.sqrt(loc[1:]*loc[:-1])
        locmidext = np.concatenate(([loc[0]*np.sqrt(loc[0]/loc[1])],
                                    locmid,
        #[24.0112]LZX: here we use this thing to be consistent with the sfd_simple
                                    [np.sqrt(loc[-1]/loc[-2])*loc[-1]]))

        #cumulative mass function is defined on the midpoints
        cummtot = np.concatenate(([0],np.cumsum(mtot)))

        #we also sample at the midpoints; then do a "diff" 
        #to get the new masses
        locmidn = np.sqrt(locn[1:]*locn[:-1])

        locmidnext = np.concatenate(([sim.rinn],locmidn,[sim.rout]))
        cummtotn = np.interp(locmidnext, locmidext, cummtot)
        mtotn = np.diff(cummtotn)
        print('[global_resample2]:mass losss = ', cummtotn[0], 100*cummtotn[0]/cummtotn[-1], '%')

        #maybe we just remove the 0 mtot particles 
        non0idx = np.argwhere(mtotn>0).flatten()
        locn = locn[non0idx]
        mtotn = mtotn[non0idx]
        locmidnext = locmidnext[np.append(non0idx, non0idx[-1]+1)]


        #with these rules we can also sample other (mass-weighted quantities)
        cummass = np.concatenate(([0], np.cumsum(mtot*mphy)))
        cummassn = np.interp(locmidnext, locmidext, cummass)
        mphyn = np.diff(cummassn) /(mtotn+1e-16) #1e-16 to prevent the zero division 

        #check for the mass conservation 
        try: 
            assert(np.abs(cummtotn[-1]-cummtot[-1])/cummtot[-1]<1e-10)
        except:
            print('[global_resample2]:mass conservation violated')
            import pdb;pdb.set_trace()

        #composition... to be tested
        fcompn = np.empty((len(non0idx),ncomp))
        for k in range(ncomp):
            cummass = np.concatenate(([0], np.cumsum(mtot*fcomp[:,k])))
            cummassn = np.interp(locmidnext, locmidext, cummass)
            fcompn[:,k] = np.diff(cummassn) /mtotn

        return locn, mtotn, mphyn, fcompn
    else:
        return None


def global_resample3 (sim, spN, fchange=0.9, fdelX=1, nsampleX=0, nspec=1, fdelDM= 1e-3, **args):
    """
    A variation on the below
    """

    loc = spN.locL 
    mtot = spN.mtotL #total mass
    mphy = spN.massL #physical mass
    fcomp = spN.fcomp #composition fraction
    ncomp = len(fcomp[0])


    #also try the keep two particles 
    nspec = 2
    idxspec = np.searchsorted(loc, sim.specloc)
    idxspec = np.concatenate(([0], idxspec, [len(loc)-1]))

    xdel = np.diff(np.log(loc))
    npar = pars.dparticleprops['nini']

    #[24.12.30]: in global_resample, better to let fdelS, fdelM tied to the initial grid spacing
    #the disired xdel
    xdel_aim = np.log((sim.rout/sim.rinn)**(1.0/npar))
    fdelS = xdel_aim/fchange
    fdelM = xdel_aim*fchange
    fspec = (fdelS+fdelM)/2

    opL, = np.where(xdel<fdelM)
    opL = np.append(opL, np.where(xdel>fdelS)[0])

    if len(opL)>0:
        print(opL, np.searchsorted(loc, sim.icelineL[0].loc))

        ## get the rdum according to the special location 
        locspecL = [] 
        #try to also treat the inner edge as a special location 
        for line in sim.planetL+sim.icelineL:
            locspecL.append(line.loc)

        #new locations -- this follows initialization and we better save them!
        rdum = np.exp(np.linspace(np.log(sim.rinn),np.log(sim.rout),npar+1))


        locn = np.sqrt(rdum[1:]*rdum[:-1])
        if len(locspecL) > 0:
            locfull = np.concatenate(([locn[0]], loc, [locn[-1]]))
            #get the index of special particles in the original locations
            specidx = np.searchsorted(locfull, locspecL) 
            specidx = np.concatenate(([0], specidx, [len(locfull)]))

            #maybe also need to apply the special range to specidx 

            locnn = np.array([])

            for i in range(len(specidx)-1): 
                numpar = len(locn[np.argwhere((locn>locfull[specidx[i]]) & (locn<locfull[specidx[i+1]-1])).flatten()])
                slice = np.exp(np.linspace(np.log(locfull[specidx[i]]), np.log(locfull[specidx[i+1]-1]), numpar))
                # slice = np.sqrt(slice_dum[1:]*slice_dum[:-1])
                locnn = np.append(locnn, slice)

        #to see what the locnn different from locn 
        # import cgs
        # y1 = np.ones_like(locn)
        # y2 = 2*np.ones_like(locnn) 
        # y3 = 3*np.ones_like(loc)
        # plt.ylim(-4,6)
        # plt.xlim(locspecL[0]-cgs.RJ, locspecL[0]+cgs.RJ)
        # plt.axvline(locspecL[0], ls='--', c = 'gray', label='iceline')
        # plt.axvline(locfull[specidx[1]], ls='--', c = 'gray', lw = 0.5)
        # plt.axvline(locfull[specidx[1]-1], ls='--', c = 'gray', lw = 0.5)
        # plt.scatter(locn, y1, label = 'locn', c='r', s=4)
        # plt.scatter(locnn, y2, label = 'locnn', s = 4)
        # plt.scatter(loc, y3, label = 'old', s=4)
        # plt.legend()
        # plt.savefig('resample_spec.jpg')
        # plt.close()



        #midpoints and extensions
        locmid = np.sqrt(loc[1:]*loc[:-1])
        locmidext = np.concatenate(([loc[0]*np.sqrt(loc[0]/loc[1])],
                                    locmid,
        #[24.0112]LZX: here we use this thing to be consistent with the sfd_simple
                                    [np.sqrt(loc[-1]/loc[-2])*loc[-1]]))

        #cumulative mass function is defined on the midpoints
        cummtot = np.concatenate(([0],np.cumsum(mtot)))


        locn = locnn

        #we also sample at the midpoints; then do a "diff" 
        #to get the new masses
        locmidn = np.sqrt(locn[1:]*locn[:-1])

        locmidnext = np.concatenate(([sim.rinn],locmidn,[sim.rout]))
        cummtotn = np.interp(locmidnext, locmidext, cummtot)
        mtotn = np.diff(cummtotn)
        print('[global_resample2]:mass losss = ', cummtotn[0], 100*cummtotn[0]/cummtotn[-1], '%')

        #maybe we just remove the 0 mtot particles 
        non0idx = np.argwhere(mtotn>0).flatten()
        locn = locn[non0idx]
        mtotn = mtotn[non0idx]
        locmidnext = locmidnext[np.append(non0idx, non0idx[-1]+1)]


        #with these rules we can also sample other (mass-weighted quantities)
        cummass = np.concatenate(([0], np.cumsum(mtot*mphy)))
        cummassn = np.interp(locmidnext, locmidext, cummass)
        mphyn = np.diff(cummassn) /(mtotn+1e-16) #1e-16 to prevent the zero division 

        #check for the mass conservation 
        merr = np.abs(cummtotn[-1] - cummtot[-1])/cummtot[-1]
        if merr>1e-5:
            print('mass conservation is violated')
            #here plot the cumulative mass and mass distribution to check what happens 
            import cgs
            fig, (axcu, axm) = plt.subplots(1,2, figsize=(10,5)) 
            axcu.loglog(locmidext/cgs.RJ, cummtot, 'o-', label='old')
            axcu.loglog(locmidnext/cgs.RJ, cummtotn, 'x-', label='new loss{:.3e}'.format(merr))
            axcu.set_ylim(cummtot[-1]*0.99, cummtot[-1]*1.01)
            axcu.set_xlim(locmidext[-1]/cgs.RJ*0.99, locmidext[-1]/cgs.RJ*1.01)

            axm.loglog(loc/cgs.RJ, mtot, '.-', label='old')
            axm.loglog(locn/cgs.RJ, mtotn, 'x-', label='new')
            axm.set_ylim(mtot[-1]*0.99, mtot[-1]*1.01)
            axm.set_xlim(loc[-1]/cgs.RJ*0.99, loc[-1]/cgs.RJ*1.01)
            axcu.axvline(pars.dgasgrid['rinn']/cgs.RJ, color='k', linestyle='--')
            axm.axvline(pars.dgasgrid['rinn']/cgs.RJ, color='k', linestyle='--')
            axcu.legend() 
            axm.legend()

            #make the ticks smaller in fontsize
            axcu.tick_params(axis='both', labelsize=2) 
            axm.tick_params(axis='both', labelsize=2)

            plt.show()
            import pdb;pdb.set_trace()

        #check the mphy around the iceline  
        # import cgs
        # plt.loglog(loc, mphy/mphy.max(), 'x-', label ='old')
        # plt.loglog(locn, mphyn/mphy.max(), 'x-', label = str(sim.time/cgs.yr))
        # plt.loglog(loc, mtot/mtot.max(), '.-', label='oldmtot')
        # plt.loglog(locn, mtotn/mtot.max(), '.-', label='newmtot')
        # plt.axvline(locspecL[0], ls='--', lw=1, c = 'gray')
        # plt.xlim(locspecL[0]*0.9, locspecL[0]*1.1)
        # plt.legend()
        # plt.savefig('/home/lzx/CpdPhysics/Test/zxl_test/mphy/{:.2f}.jpg'.format(sim.time))
        # plt.close()

        #composition... to be tested
        fcompn = np.empty((len(non0idx),ncomp))
        for k in range(ncomp):
            cummass = np.concatenate(([0], np.cumsum(mtot*fcomp[:,k])))
            cummassn = np.interp(locmidnext, locmidext, cummass)
            fcompn[:,k] = np.diff(cummassn) /mtotn

        return locn, mtotn, mphyn, fcompn
    else:
        return None

def global_resample (sim, spN, fchange=0.9, fdelX=1, nsampleX =0, nspec = 1,**args):
    """
    Follow the Schoonenberg. 2018 Lagrangian model. 

    Basic idea: resample all the particles to the initial configuration when conditions 
    are met.

    nspec: 
        The number of particial bond with the special locations 

    Conditions:
        Currently I would like to follow the conditions above 

    Resampling:
        1. Firstly identify the special locations (e.g. iceline, planet location), the particles 
            in these special locations will not be resampled.
        2. Resample principle:
            
            Idea 1 (rejected): 
                1. first get the locations 
                2. modify the locations: delete the locations that close to the special locations within 
                    a certain range, and insert the surrounding particles into the locations.
                3. we consider the interpolation according to the locations of icelines
                4. for total mass and physical mass, interpolate the mass of the slice.
                5. for composition, just use the composition of that region.

            Idea 2: 
                resample the total mass by the cumulative mass fraction of the particles.
    """
    loc = spN.locL 
    marr = spN.mtotL

    xdel = np.diff(np.log(loc))
    npar = pars.dparticleprops['nini']

    #[24.12.30]: in global_resample, better to let fdelS, fdelM tied to the initial grid spacing
    #the disired xdel
    xdel_aim = np.log((sim.rout/sim.rinn)**(1.0/npar))
    fdelS = xdel_aim/fchange
    fdelM = xdel_aim*fchange
    fspec = (fdelS+fdelM)/2

    opL, = np.where(xdel<fdelM)
    opL = np.append(opL, np.where(xdel>fdelS)[0])

    if len(opL)>0:
        #get new locations
        locn = sim.rinn*(sim.rout/sim.rinn)**np.linspace(1/npar, 1, npar)

        #get the special locations 
        #[24/11/18] b/c we have consider the iceline location later, so here maybe we can consider the 
        #           planet location only.
        locspecL = [] 
        #try to also treat the inner edge as a special location 
        locspecL.append(sim.rinn)
        for line in sim.planetL+sim.icelineL:
            locspecL.append(line.loc)

        # not sure whether it's a good idea to add the outermost location of particles to the special location 
        # locspecL.append(loc[-1])

        #to prevent there are some overlapping special locations 
        locspecL = list(set(locspecL))


        #consider the special locations like I did before. Just first add the 1 surrounding particles 
                #find the locations that close to the special locations 
        specpar_idx = np.searchsorted(loc, locspecL) 
        #insert the surrounding particles into locn 
        #we also need to remove the locations that closer to the special locations in locn 

        #LZX[24/11/2] not sure whether we should get the 'special particles' within two particles or within a 
        #             certain location range.
        for i, idx in enumerate(specpar_idx):
            #maybe we should combine the two methods together. [24/11/27]
            #first get the special particles, and also special range, choose the minimal range to implement
            #inside of special range, we regard the particles as 'effective particles for special locations'
            if idx < nspec:
                specL_par = np.array(loc[:idx+nspec])
            elif idx>len(loc)-nspec:
                specL_par = np.array(loc[idx-nspec:])
            else:
                specL_par = np.array(loc[idx-nspec:idx+nspec]) 

            #'special range' 
            speclim = [locspecL[i]/np.exp(fspec), locspecL[i]*np.exp(fspec)]

            #'crop' the specL_par 
            specL_par = specL_par[np.where((specL_par>speclim[0]) & (specL_par<speclim[1]))[0]]


            #particles in special range 
            #specL_par = np.where((loc>speclim[0]) & (loc<speclim[1]))[0]

            if len(specL_par)>0:
                locn = np.delete(locn, np.where((locn>specL_par.min()) & (locn<specL_par.max())))
                #locn = np.delete(locn, np.where((locn>speclim[0]) & (locn<speclim[1])))
                locn = np.insert(locn, np.searchsorted(locn, specL_par), specL_par)
                #check: if the difference of the locations are too small, 
                #       just remove! 
                locn = np.delete(locn, np.where(np.diff(np.log(locn))<fdelM))

        #[24.01.07]LZX: don't resample the particles outer than the final one 
        #locn = np.append(locn[locn<loc[-2]/np.exp(fdelM)], [loc[-2],loc[-1]])

        #get the iceline locations and get the new properties according to them 
        # iceloc = [] #[i.loc for i in sim.icelineL]
        # slice_idxo = np.searchsorted(loc, iceloc)
        # slice_idxo = np.concatenate(([0], slice_idxo, [len(loc)]))
        #
        # slice_idxn = np.searchsorted(locn, iceloc) 
        # slice_idxn = np.concatenate(([0], slice_idxn, [len(locn)]))

        mphyo = spN.massL 

        #[24/11/27]: 
        #if the innermost particle is inside of the rinn+deltar/2 
        #we add the rinn and rinn-deltar/2  to the beginning of the locn 
        #and add the rinn-deltar/2 to the beginning of the loc 

        #add_virtual = False #loc[0]<locn[0]

        # if add_virtual:
        #     r0 = sim.rinn/(locn[1]/locn[0])
        #
        #     locn = np.append(sim.rinn, locn) 
        #     locn = np.append(r0, locn)
        #
        #     loc = np.append(r0, loc) 
        #     marr = np.append(0, marr) 
        #     mphyo = np.append(mphyo[0], mphyo)
            

        #initialize the new property arrays
        mtotn = np.array([]) 
        massn = np.array([])

        #[24/11/27]resample according to iceline is not a good idea...
        cum_mtoto = np.cumsum(marr)
        cum_mtotn = np.interp(locn, loc, cum_mtoto)
        mtotn = np.append(cum_mtotn[0], np.diff(cum_mtotn)) 
        massn = np.interp(locn, loc, mphyo) 
        
        #cumdum = np.cumsum(marr*mphyo)
        #dum1 = np.interp(locn, loc, cumdum)
        #massn = np.append(dum1[0], np.diff(dum1))/mtotn

        #remove the first 'virtual' particle 
        # if add_virtual: 
        #     locn = locn[2:]
        #     mtotn = mtotn[2:]
        #     massn = massn[2:]


        fcompn = np.zeros((len(locn), len(sim.particles.fcomp[0])))
        icelineL = np.array([i.loc for i in sim.icelineL])
        idxn = np.searchsorted(locn, icelineL) 
        idxo = np.searchsorted(loc, icelineL) 

        idxn = np.append(0, idxn) 
        idxn = np.append(idxn, len(locn))
        idxo = np.append(0, idxo) 

        for i in range(len(idxn)-1):
            fcompn[idxn[i]:idxn[i+1]] = sim.particles.fcomp[idxo[i]] 


        # import cgs
        # plt.loglog(loc, mphyo/mphyo.max(), 'x-', label ='old')
        # plt.loglog(locn, massn/massn.max(), 'x-', label = str(sim.time/cgs.yr))
        # plt.loglog(loc, marr/marr.max(), '.-', label='oldmtot')
        # plt.loglog(locn, mtotn/mtotn.max(), '.-', label='newmtot')
        # plt.axvline(locspecL[1], ls='--', lw=1, c = 'gray')
        # plt.xlim(locspecL[1]*0.9, locspecL[1]*1.1)
        # plt.legend()
        # plt.savefig('/home/lzx/CpdPhysics/Test/zxl_test/mphy/{:.2f}.jpg'.format(sim.time))
        # plt.close()
        #
        #check the mass conservation every slices
        merr = np.abs(cum_mtotn[-1] - cum_mtoto[-1])
        if merr>1e-10:
            print('mass conservation is violated')
            #here plot the cumulative mass and mass distribution to check what happens 
            import cgs
            fig, (axcu, axm) = plt.subplots(1,2, figsize=(10,5)) 
            axcu.loglog(loc/cgs.RJ, cum_mtoto, 'o-', label='old')
            axcu.loglog(locn/cgs.RJ, cum_mtotn, 'x-', label='new')
            axm.loglog(loc/cgs.RJ, marr, '.-', label='old')
            axm.loglog(locn/cgs.RJ, mtotn, 'x-', label='new')
            axcu.axvline(pars.dgasgrid['rinn']/cgs.RJ, color='k', linestyle='--')
            axm.axvline(pars.dgasgrid['rinn']/cgs.RJ, color='k', linestyle='--')
            axcu.legend() 
            axm.legend()
            plt.show()
            import pdb;pdb.set_trace()       

        return locn, mtotn, massn, fcompn
    else:
        return None




if False:
    if len(opL)>0:
        #if the particles are too close or too far away, just resample 

        #locations 
        npar = pars.dparticleprops['nini']
        locn = sim.rinn*(sim.rout/sim.rinn)**np.linspace(1/npar, 1, npar)

        #find the special locations 
        locspecL = [] 
        for line in sim.icelineL+sim.planetL:
            locspecL.append(line.loc)

        #find the locations that close to the special locations 
        specpar_idx = np.searchsorted(loc, locspecL) 
        #insert the surrounding particles into locn 
        #we also need to remove the locations that closer to the special locations in locn 

        #LZX[24/11/2] not sure whether we should get the 'special particles' within two particles or within a 
        #             certain location range.
        for idx in specpar_idx:
            #specL_par = np.array(loc[idx-nspec:idx+nspec]) 
            #'special range'
            speclim = [loc[idx]*(10**fspec-1)/10**fspec, (10**fspec-1)*loc[idx]]

            #particles in special range 
            specL_par = np.where((loc>speclim[0]) & (loc<speclim[1]))[0]

            #locn = np.delete(locn, np.where((locn>specL_par.min()) & (locn<specL_par.max())))
            locn = np.delete(locn, np.where((locn>speclim[0]) & (locn<speclim[1])))
            locn = np.insert(locn, np.searchsorted(locn, specL_par), specL_par)
            #check: if the difference of the locations are too small, 
            #       just remove! 
            locn = np.delete(locn, np.where(np.diff(locn)<fdelM))

        #initialize the new property arrays 
        locI = locn 
        massI = np.array([])
        mtotI = np.array([]) 
        fcompI = np.zeros((len(locn), len(sim.particles.fcomp[0])))

        #total mass
        illocs = np.array([i.loc for i in sim.icelineL])
        #get the slice of new locations 
        slice_idxn = np.searchsorted(locn, illocs)
        locn_slice = np.split(locn, slice_idxn)
        slice_idxn = np.concatenate(([0], slice_idxn, [len(locn)]))

        n_slice = illocs+1 # n icelines, n+1 slices 
        
        #get the slices of old properties 
        slice_idx = np.searchsorted(loc, illocs)
        slice_idx = np.concatenate(([0], slice_idx, [len(loc)]))

        wdeln = np.concatenate(([np.log(locn[1]/locn[0])], 
                               (np.log(locn[2:]) -np.log(locn[:-2]))/2,
                               [np.log(locn[-1]/locn[-2])]))
        #wdelo = np.concatenate(([np.log(loc[1]/loc[0])], 
        #                       (np.log(loc[2:]) -np.log(loc[:-2]))/2,
        #                       [np.log(loc[-1]/loc[-2])]))

        for i,locl in enumerate(locn_slice):
            #get the mass of the slice 
            m_slice = marr[slice_idx[i]:slice_idx[i+1]]
            sfd_slice = sfd[slice_idx[i]:slice_idx[i+1]]
            summlice = m_slice.sum()

            loco_slice = loc[slice_idx[i]:slice_idx[i+1]]

            #get the physical mass of the slice 
            mphy_slice = spN.massL[slice_idx[i]:slice_idx[i+1]]

            #insert the slice into the new locations 
            sfd_slicen = np.interp(locl, loco_slice, sfd_slice) 

            #get the total mass of the slice
            mn_slice = sfd_slicen*2*np.pi*locl**2*wdeln[slice_idxn[i]:slice_idxn[i+1]] #mass of the slice

            #maybe normalized to the old total mass???

            #check the mass conservation 
            if np.abs(mn_slice.sum()/summlice-1)>2e-2:
                print('mass conservation is violated') 
                import pdb;pdb.set_trace()

            mtotI = np.append(mtotI, mn_slice)

            #get the physical mass of the slice 
            mphy_slicen = np.interp(locl, loco_slice, mphy_slice)
            massI = np.append(massI, mphy_slicen) 

            #get the composition of the slice 
            comp = sim.particles.fcomp[slice_idx[i]]
            fcomp_slicen = np.zeros((len(locl), len(comp))) 
            fcomp_slicen[:] = comp
            fcompI[slice_idxn[i]:slice_idxn[i+1]] = fcomp_slicen
        #return locI, mtotI, massI, fcompI
    
    #else:
#        return None


