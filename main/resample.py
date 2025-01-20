import numpy as np
import parameters as pars
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import functions as ff
import physics


def v_rad (marr,disk,rhoint, Stg=None, vadd=0):
    """
    obtain radial velocity of particles
    - radarr    :input locations
    - disk      :disk object at locations
    - rhoint    :internal densities
    """

    sarr = physics.mass_to_radius(marr,rhoint)
    St, vr = ff.Stokes_number(disk, sarr, rhoint, Sto=Stg)

    import pdb; pdb.set_trace()
    return vr + vadd


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
    #imL = np.array([],dtype=np.int64) #no merging


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

        mflux = sim.particles.v_r *sim.particles.sfd *sim.particles.locL
        if np.any(np.diff(mflux[1:20])<0) and False:
            print('mflux not in order')


        pm = 0.5
        dum = interp_mtot_weighted (locmidnext, locmidext, mphy**pm, mtot, mtotn)
        mphyn = dum**(1/pm)
        #logmphyn = interp_mtot_weighted (locmidnext, locmidext, np.log(mphy), mtot, mtotn)
        #mphyn = np.exp(logmphyn)

        #the desired velocities
        vraim = (sim.particles.v_r[isL] +sim.particles.v_r[isL+1])/2

        ## search for proper particles mass
        disk = sim.get_disk (loc=addlocS)
        rhoint = sim.particles.rhoint[isL] #this should be changed
        stgarr = sim.particles.St[isL]
        out = v_rad (mphy[isL],disk,rhoint, stgarr, vadd=0)

        import pdb; pdb.set_trace()
        #
        # def func(mass,disk,rhoint):
        #       return vr

        if loc[0]<sim.rinn:
            import pdb; pdb.set_trace()

        if np.all(np.diff(np.log10(mphyn[:100]))<0)==False and False:
            print('physical mass not in order')

        sfdnew = ff.sfd_special (mtotn, locn, sim.specloc)

        return locn, mtotn, mphyn, fcompn
    else:
        return None

def new_splitmerge_zxl (sim, spN, fdelS, fdelM=0., fdelX=1, nsampleX=0, fdelDM=0.0001):
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

    #midpoint locations elligible for splitting (isL)
    #don't split 1st/last particles (i.e., around special locations; TBD)
    #fdelXarr[0] = np.inf; fdelXarr[-1] = np.inf


    #now the masses
    #ydel = np.diff(np.log(mphy))
    #isL2, = np.nonzero(np.abs(ydel)>0.5*fdelXarr)
    #isL = np.union1d(isL1, isL2)

    #merging
    #fdelXarr[0] = 0; fdelXarr[-1] = 0
    fdelXarr = np.ones_like(xdel)
    imL, = np.nonzero(xdel<fdelXarr*sim.particles.delta*2/3)


    fdelXarr[imL] = np.inf #dont split where we merge
    fdelXarr[imL+1] = np.inf
    fdelXarr[imL-1] = np.inf

    fdelS = 2*sim.particles.delta
    isL, = np.nonzero(xdel>fdelS*fdelXarr)

    #splitting: add the locations
    addlocS = np.sqrt(loc[isL]*loc[isL+1])
    addlocM = np.sqrt(loc[imL]*loc[imL+1])

    if len(isL)>0 or len(imL)>0:
        doResample = True
    else:
        doResample = False


    if False and len(imL)>0:
        addloc = np.sqrt(loc[imL]*loc[imL+1])
        addmtot = mtot[imL]+mtot[imL+1]
        addmphy = np.sqrt(mphy[imL]*mphy[imL+1])

        #addmphy = (mphy[imL]*mtot[imL] +mphy[imL+1]*mtot[imL+1]) /addmtot
        #addmloc = (loc[imL]*mtot[imL] +loc[imL+1]*mtot[imL+1]) /addmtot

        locn = loc.copy()
        locn[imL] = addloc
        locn = np.delete(locn,imL+1)
        npar = len(locn)

        mtotn = mtot.copy()
        mtotn[imL] = addmtot
        mtotn = np.delete(mtotn,imL+1)

        mphyn = mphy.copy()
        mphyn[imL] = addmphy
        mphyn = np.delete(mphyn,imL+1)

        fcompn = np.ones((npar,ncomp))

        sfdnew = ff.sfd_special (mtotn, locn, sim.specloc)
        return locn, mtotn, mphyn, fcompn


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

        #print(isL)
        mflux = sim.particles.v_r *sim.particles.sfd *sim.particles.locL
        if np.any(np.diff(mflux[1:20])<0) and False:
            print('mflux not in order')
            import pdb; pdb.set_trace()

        pm = 0.5
        dum = interp_mtot_weighted (locmidnext, locmidext, mphy**pm, mtot, mtotn)
        mphyn = dum**(1/pm)

        ##[25.01.20]lzx: now from the velocity to get the mass 
        #first get the velocity in the geometrical average point 
        addloc = np.append(addlocS, addlocM)
        for i, idx in enumerate(np.append(isL, imL)):
            v1 = sim.particles.v_r[idx] 
            v2 = sim.particles.v_r[idx+1]
            r1 = loc[idx]
            r2 = loc[idx+1]
            vnew = (r1*v2+r2*v1)/2/addloc[i]

            #get the mass from the velocity

        import pdb;pdb.set_trace()
        #logmphyn = interp_mtot_weighted (locmidnext, locmidext, np.log(mphy))
        #mphyn = np.exp(logmphyn)
        #import pdb; pdb.set_trace()

        #mphyn[isL+1] = (0.5*(mphy[isL]**(1/3)+mphy[isL+2]**(1/3)))**3

       #if np.all(np.diff(np.log10(mphyn[:100]))<0)==False:
       #    print('physical mass not in order')

        sfdnew = ff.sfd_special (mtotn, locn, sim.specloc)

        #if len(imL)>0: import pdb; pdb.set_trace()

        #composition... to be tested
        fcompn = np.empty((npar,ncomp))
        for k in range(ncomp):
            cummass = np.concatenate(([0], np.cumsum(mtot*fcomp[:,k])))
            cummassn = np.interp(locmidnext, locmidext, cummass)
            fcompn[:,k] = np.diff(cummassn) /mtotn

        return locn, mtotn, mphyn, fcompn
    else:
        return None


def v2mass (vr, loc, sim):
    """
    [25.01.20]: get the mass from the velocity
    """
    from scipy.optimize import fsolve 
    pass



def resamplr_spiltmerge_zxl(sim, fchange, *args):
    """
    do the split and merge by the cumulative total mass distrbution
    """
    loco = sim.particles.locL 
    mtoto = sim.particles.mtotL 
    masso = sim.particles.massL 
    fcompo = sim.particles.fcomp 

    #let's forget about the special locations for now [25.01.18]lzx 
    xdel = np.log(loco[1:]/loco[:-1])

    fdelS = sim.particles.delta/fchange
    fdelM = sim.particles.delta*fchange

    isL = np.where(xdel>fdelS)[0] 
    imL = np.where(xdel<fdelM)[0]

    if len(isL)>0: 
        for id in isL:
            if id>0 and id<len(loco)-2:
                locn = np.insert(loco, id+1, np.sqrt(loco[id]*loco[id+1]))

                #do the cumulative mass distribution to find the total mass 
                cummasso = np.cumsum(mtoto)
                cummassn = np.interp(locn, loco, cummasso) 
                mtotn = np.diff(cummassn)




def re_sample_splitmerge (sim, spN, fdelS, fdelM=0., fdelX=1, nsampleX=0, fdelDM=0.0001):
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


def interp_mtot_weighted (xn, xmid, qarr, marr=None, mn=None, neval=0):
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
    qn = np.diff(cummassn) /(mn+1e-16) #1e-16 to prevent the zero division 

    #print(neval, xn[0], len(xn))

    #relative error/significance
    sig = np.diff(cummassn)/cummassn[1:]
    itrust = (sig<1e-8).argmax() #trust until here
    if itrust>0:
        qn[itrust:] = interp_mtot_weighted (xn[itrust:], xmid, qarr, marr, mn[itrust:], neval)

    #if neval==1: import pdb; pdb.set_trace()
    return qn


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
                logmphyn = interp_mtot_weighted (locmidnext, locmidext, np.log(mphy[ii]), mtot[ii], mtotn)
                mphyn = np.exp(logmphyn)
                #mphyn = interp_mtot_weighted (locmidnext, locmidext, mphy[ii])

                #this is quite diffusive
                #pwl = 0.5
                #dum = interp_mtot_weighted (locmidnext, locmidext, mphy[ii]**pwl, mtot[ii], mtotn)
                #mphyn = dum**(1/pwl)

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


