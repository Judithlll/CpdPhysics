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
    disk.add_auxiliary ()
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
    vt1 = bt1*(planetMass/mcp)*(disk.sigmaG*planetLoc**2/mcp)*(disk.vK/disk.cs)**2 *disk.vK

    v_mig=vt1+vd
    return v_mig


