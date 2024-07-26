import numpy as np
import parameters as pars
import sys


def re_sample_splitmerge (sim, spN, fdelS, fdelM=0., fdelX=1, nsampleX=0, fdelDM=0.001):
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


