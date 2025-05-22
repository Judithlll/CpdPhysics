import numpy as np

eta_default = 1.8 # default smoothing length


#NOTE: maybe the temperature feedback also matter? 
# 1. to use a higher order kernel filling, 
# 2. try to wisely ease the tension between this and resample

def kernel_fill (qM):

    if np.any(qM<0):
        print('[newkernel.py]warning: qM<0')
        import pdb; pdb.set_trace()
    i0 = qM<1
    i1 = qM<2

    wM = np.zeros_like(qM)
    wpM = np.zeros_like(qM)

    #kernel matrix
    wM[i1]  = 1/4 *(2-qM[i1])**3 
    wM[i0] -=      (1-qM[i0])**3

    #kernal derivative matrix
    wpM[i1]  = -3/4 *(2-qM[i1])**2
    wpM[i0] +=  3   *(1-qM[i0])**2

    sig = 2/3 #1D
    wM *= sig
    wpM *= sig
    return wM, wpM

def kernel_fill_M5 (qM):
    """
    try the higer order filling
    """
    if np.any(qM<0):
        print('[newkernel.py]warning: qM<0')
        import pdb; pdb.set_trace()

    i0 = qM<0.5
    i1 = qM<1.5 
    i2 = qM<2.5 

    wM = np.zeros_like(qM) 
    wpM = np.zeros_like(qM) 

    #kernel matrix wM 
    wM[i2]  = (5/2 - qM[i2])**4 
    wM[i1] -= 5*(3/2 - qM[i1])**4
    wM[i0] += 10*(1/2 - qM[i0])**4

    #kernel derivative matrix wpM 
    wpM[i2]  = -4*(5/2 - qM[i2])**3
    wpM[i1] += 20*(3/2 - qM[i1])**3
    wpM[i0] -= 40*(1/2 - qM[i0])**3


    sig = 1/24 #1D for M5 quartic 
    wM *= sig
    wpM *= sig 

    return wM, wpM

def kernel_fill_M6 (qM):
    if np.any(qM<0):
        print('[newkernel.py]warning: qM<0')
        import pdb; pdb.set_trace()

    i0 = qM<1
    i1 = qM<2
    i2 = qM<3

    wM = np.zeros_like(qM) 
    wpM = np.zeros_like(qM) 

    #kernel matrix wM 
    wM[i2]  = (3 - qM[i2])**5 
    wM[i1] -= 6*(2 - qM[i1])**5
    wM[i0] += 15*(1 - qM[i0])**5

    #kernel derivative matrix wpM 
    wpM[i2]  = -5*(3 - qM[i2])**4 
    wpM[i1] += 30*(2 - qM[i1])**4 
    wpM[i0] -= 75*(1 - qM[i0])**4


    sig = 1/120 #1D for M6 quintic 
    wM *= sig
    wpM *= sig 

    return wM, wpM



def get_weight (loc, face, hguess=None, eta=eta_default):
    """
    we use an initial guess h0 and solve for x
        h = h0*x
        q = q0/x;   q0 = d_ij /h0

    matrix q0[i,:] gives the distance array for particle i
    such that matrix

        wM[i,:] = w(q)

    is the (normalized) kernel barring the softening length.

    We find the root of the equation: 

        y_i = sum(wM[i,:]) -eta

    through a Newton-Rhapson method to solve for x_i

    Hence we obtain:
        - the softenings lengths h
        - the kernel matrix wM

    from which the density is obtained
    """

    #initial guess (if not given)
    if hguess is None or len(hguess)!=len(loc):
        h0 = np.diff(face)
        # h0 = np.concatenate(( [2*(loc[1]-loc[0])], 
        #                         loc[2:] -loc[:-2],
        #                     [2*(loc[-1]-loc[-2])] ))
    else:
        h0 = hguess.copy()


    #initial guess for q distance matrix
    #I swap the axes such that q0[i,:] are the
    #distances to particle i
    q0 = (np.abs(loc -loc[:,np.newaxis]) /h0).T

    xarr = np.zeros_like(h0)
    yarr = 1e99 *np.ones_like(h0)
    ii = np.ones_like(xarr, dtype=bool)

    iloop = 0
    iloop_max = 20
    while iloop<iloop_max:

        if iloop==0:
            qM = q0.copy()
            wM, wM_q = kernel_fill_M6(qM)
        else:
            #decrease the q, to account for further particles.
            qM = q0[ii,:] /(1+xarr[ii][:,np.newaxis])

            wM[ii,:], wM_q = kernel_fill_M6(qM)
            #set all the minus value to zero 
            wM[ii,:] = np.maximum(0,wM[ii,:])


        yval =    (wM[ii,:]).sum(1) -eta
        yprime = -(wM_q *qM/(1+xarr[ii][:,np.newaxis])).sum(1)

        if iloop>5:
            print(iloop, yval[0], xarr[0])


        #limit the change somewhat
        # xarr[ii] = xarr[ii] -np.maximum(-0.4,np.minimum(0.4, yval/yprime))

        delx_mch = 0.1 +0.8*np.abs(xarr[ii])
        delx = np.minimum(delx_mch, np.maximum(-delx_mch,-yval/yprime))
        xarr[ii] += delx

        if np.any(yprime==0):
            import pdb;pdb.set_trace()

        yarr[ii] = yval

        ii = abs(yarr)>1e-8 #indices still involved
        if ii.sum()==0: 
            break
        else:
            iloop += 1

    if iloop>=20:
        print('[newkernel.py]warning:iloop=10 reached!')
        import pdb; pdb.set_trace()

    #finally, the softening lengths..
    hsoftarr = h0 *(1+xarr)
    return hsoftarr, wM

