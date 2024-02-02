import numpy as np
import cgs
import functions as f
import disk_properties as dp
import physics

def planet_migration (gas,planetLoc,planetMass,time,rhopl):
    #aerodynamic drag migration
    out = gas.get_key_disk_properties (planetLoc, time)

    mcp = dp.Mcp_t(time)

    disk = physics.DISK (*out, planetLoc, time, mcp) #pro
    disk.add_auxiliary ()
    disk.add_uservar (dp.user_add_var())    #variables
    disk.add_userfun (dp.user_add_fun())
    disk.add_user_eval (dp.user_add_eval())
    #disk = core.DISK (*out, planetLoc, time)
    #disk.add_auxiliary ()
    #disk.user_difined ()

    # eta=disk.eta
    # v_K=disk.vK
    # v_th=disk.vth
    # lmfp=disk.lmfp
    # rho_g=disk.rhog
    # Omega_K=disk.OmegaK
    # dotMg=disk.dotMg
    # Mcp=disk.Mcp
    # mg=disk.mu*cgs.mp
    # Sigmag=disk.sigmaG
    # cs=disk.cs
    
    # Rpl=(planetMass/(4/3*np.pi*rhopl))**(1/3)

    # Stpl,vd=f.St_iterate(eta,v_K,v_th,lmfp,rho_g,Omega_K,Rpl)

    #Type I migration
    qr=-0.14493*disk.dot_Mg(disk.time)*cgs.gC*disk.mcp*(-0.206349*planetLoc**(5.5)+planetLoc**4.5*disk.rout)/disk.alpha/planetLoc**(8.5)/disk.rout/np.sqrt(cgs.kB*(disk.dot_Mg(time)*cgs.gC*disk.mcp/planetLoc**3/cgs.sigmaSB)**(1/4)/disk.mg) #pressure gradient

    CI=0.1
    bt1=CI*(2.7+1.1*qr)   #a constant Ogihara 2014
    vt1=bt1*(planetMass/disk.mcp)*(disk.sigmaG*planetLoc**2/disk.mcp)*(disk.vK/disk.cs)**2*disk.vK

    # v_mig=vt1+vd

    return vt1 # v_mig

        
