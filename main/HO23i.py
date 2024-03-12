import numpy as np

def al_j(j):
    j1 = j+1
    return (j1/j)**(-2/3)


def calc_laplace (al, s, m, p=1, q=1):
    from pylaplace import LaplaceCoefficient
    """
    calculate laplace coefficients for given j
    assuming j+1:j
    """

    #use Brute mode in order to prevent spurious sign change
    laplace = LaplaceCoefficient(method='Brute')

    lp = laplace (al,s,m,p,q)
    lp_d = laplace.derivative(al,s,m,p,q)

    return lp, lp_d



def calc_f1 (j):
    """
    definition follows eq.7 of Terquem & Papaloizou
    """
    j1 = j+1
    al = al_j(j)

    lp, lp_d = calc_laplace(al, 0.5, j1)

    f1 = -0.5 *(2*j1*lp +al*lp_d)
    return f1


def calc_f2 (j):
    """
    definition follows eq.8 of Terquem & Papaloizou
    """
    j1 = j+1
    al = al_j(j)

    lp, lp_d = calc_laplace(al, 0.5, j)

    f2 = 0.5 *((2*j+1)*lp +al*lp_d)
    return f2


def get_f1 (j):
    """
    a calculation of the above
    """
    df1 = {1: -1.1904936978494907,
            2: -2.0252226899386003,
            3: -2.8404318567218776,
            4: -3.649618244115698,
            5: -4.4561427850589315,
            6: -5.261253831888332,
            7: -6.065523627134226,
            8: -6.869251916654784,
            9: -7.672611034804858}

    return df1[j]


def get_f2p (j):
    """
    get the modified f2 coefficient
    (this equals f2 except for j=1)
    """
    df2p = {1: 0.428389834143891,
            2: 2.4840051833039474,
            3: 3.2832567218222226,
            4: 4.083705371796215,
            5: 4.884706297500403,
            6: 5.6860074115090455,
            7: 6.487489727978894,
            8: 7.289089771449692,
            9: 8.09077059747424}

    return df2p[j]


def ta_crit (j, qp, innerperturber=True, tau_e=None, te_over_ta=None, tau_wave=None,
                haspect=0.1, gamI=1.9, Ce=0.11, tPer=None):
    """
    Calculates the critical migration time, above which a migrating planet will get 
    trapped in resonance due to the resonance forcing of a pertuber planet.

    This follows eq.26 of Huang & Ormel (2023), which we generalized to also cover
    the case of an inner pertuber.

    !NOTE  
        default Ce parameter has been updated to 0.11, after
        fixing a mistake in the pusblished paper. See the erratum.

        in comments here any "tau" refers to dimensionless time
        (times multiplied by orbital frequency), "t"-s are real times.

    INPUT 
        j:          resonance index (j+1:j resonance)
        qp:         mass-to-central (stellar) mass for the perturber
                    (the more massive planet, assumed to move on near-circular orbit)

        innerperturber: 
                    whether the perturber planet is an inner or outer one
                    in case of an inner perturber all quantities other than qp
                    (tau_e... tPer) concern the outer planet. Vice versa for
                    an outer perturber (innerperturber=False)
                    HO23 considered an outer perturber

       [tau_e]:     the dimensionless eccentricity damping parameter
       [tau_wave]:  the dimensionless wave damping parameter
                    (see eq. 8 of HO23)
       [te_over_ta]:if provided, use this value instead of calculating 
                    it using the isothermal disk profile
       [haspect]:   aspect ratio at the location of the migrating planet
       [tPer]:      orbital period at the location of the migrating planet


    RETURNS 
        the critical migration threshold ta_crit, either in its dimensionles 
        form (multiplied by the orbital frequence) or dimensional when the 
        orbital period tPer is provided.
    """
    j1 = j+1
    al = al_j(j)


    if tau_wave is not None:
        tau_e = Ce*tau_wave/0.78
    elif te_over_ta is None:
        #eq. 9 of HO23
        te_over_ta = Ce/(0.78*gamI) *haspect**2

    #after this we either have an eccentricity damping timescale (te)
    #or the damping ratio te/ta

    if innerperturber:
        if j<=9:
            f2p = get_f2p(j)
        else:
            f2p = calc_f2(j)
        ffac = f2p
        jj = j
    else:
        if j<=9:
            f1 = get_f1(j)
        else:
            f1 = calc_f1(j)
        ffac = al*np.abs(f1)
        jj = j1

    #eq.26 of HO23
    if tau_e is not None:
        tau_crit = 1/(2*jj*ffac**2*qp**2*tau_e)
    else:
        tau_crit = 1/np.sqrt(2*j1*te_over_ta) /(ffac*qp)

    if tPer is not None:
        return tau_crit/np.sqrt(2*np.pi) *tPer
    else:
        return tau_crit


#just to test it out
if __name__=='__main__':
    for j in range(1,15):
        f2p = calc_f2 (j)
        if j==1:
            f2p = f2p -2**(1/3)

        tacrit1 = ta_crit(j, 1e-3, tau_wave=1e-2)
        tacrit1 = ta_crit(j, 1e-3, te_over_ta=1e-4)
        tacrit2 = ta_crit(j, 1e-3, Ce=0.1)
        print('{:4d} {:10.2e} {:10.2e}'.format(j, tacrit1, tacrit2))
