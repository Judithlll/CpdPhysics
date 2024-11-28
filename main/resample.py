import numpy as np
import parameters as pars
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


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


def global_resample(sim, spN, fdelS, fdelM, fdelX=1, nsampleX =0, nspec = 1, fspec= 0.002,**args):
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


    opL, = np.where(xdel<fdelM)
    opL = np.append(opL, np.where(xdel>fdelS)[0])

    npar = pars.dparticleprops['nini']

    if len(opL)>0:
        #get new locations
        locn = sim.rinn*(sim.rout/sim.rinn)**np.linspace(1/npar, 1, npar)
        
        #get the special locations 
        #[24/11/18] b/c we have consider the iceline location later, so here maybe we can consider the 
        #           planet location only.
        locspecL = [] 
        for line in sim.planetL+sim.icelineL:
            locspecL.append(line.loc)

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
            specL_par = np.array(loc[idx-nspec:idx+nspec]) 

            #'special range' 
            speclim = [locspecL[i]/10**fspec, locspecL[i]*10**fspec]

            #'crop' the specL_par 
            specL_par = specL_par[np.where((specL_par>speclim[0]) & (specL_par<speclim[1]))[0]]

            #particles in special range 
            #specL_par = np.where((loc>speclim[0]) & (loc<speclim[1]))[0]

            locn = np.delete(locn, np.where((locn>specL_par.min()) & (locn<specL_par.max())))
            #locn = np.delete(locn, np.where((locn>speclim[0]) & (locn<speclim[1])))
            locn = np.insert(locn, np.searchsorted(locn, specL_par), specL_par)
            #check: if the difference of the locations are too small, 
            #       just remove! 
            locn = np.delete(locn, np.where(np.diff(np.log(locn))<fdelM))

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

        add_virtual = loc[0]<locn[0]

        if add_virtual:
            r0 = sim.rinn/(locn[1]/locn[0])

            locn = np.append(sim.rinn, locn) 
            locn = np.append(r0, locn)

            loc = np.append(r0, loc) 
            marr = np.append(0, marr) 
            mphyo = np.append(mphyo[0], mphyo)
            

        #initialize the new property arrays
        mtotn = np.array([]) 
        massn = np.array([])

        #[24/11/27]resample according to iceline is not a good idea...
        cum_mtoto = np.cumsum(marr)
        cum_mtotn = np.interp(locn, loc, cum_mtoto)
        mtotn = np.append(cum_mtotn[0], np.diff(cum_mtotn)) 
        
        massn = np.interp(locn, loc, mphyo) 

        #remove the first 'virtual' particle 
        if add_virtual: 
            locn = locn[2:]
            mtotn = mtotn[2:]
            massn = massn[2:]


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


