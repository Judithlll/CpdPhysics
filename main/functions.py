import numpy as np
error=1e-8

def St_iterate(eta,v_K,v_th,lmfp,rho_g,
                Omega_K,R_d,rhoint=1.4,Sto=0.001):
        """
        obtain Stokes number and drift velocity by fixed-point iteration
        """
        #put in physics.py
    
        St=Sto
        nither=1
        while(nither<40):
            nither += 1
            v_r=-2*St/(St**2+1)*eta*v_K
            v_dg=np.abs(v_r)
            Rep=4*R_d*v_dg/v_th/lmfp
            #following Shibaike, eq....
            CD=24/Rep*(1+0.27*Rep)**0.43+0.47*(1-np.exp(-0.04*Rep**0.38))
            Stn=8/3/CD*rhoint*R_d/rho_g/v_dg*Omega_K
            delta_St=abs(St-Stn)
            
            if np.max(delta_St)<error:
                break
            else:
                St=Stn

        return Stn, v_r


def determine_type (val):
    """
    this determines the type of val

    hist
        [23.12.03]:copied from /NewLagrange by Chris
    """

    if val=='None':
        return None
    #boolean value
    elif val=='True':
        return True
    elif val=='False':
        return False

    #integer
    elif val.isdigit():
        return int(val)

    else:
        #float/string value
        try:
            out = float(val)
        except:
            out = val

        return out
