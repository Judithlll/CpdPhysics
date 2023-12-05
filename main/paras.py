import numpy as np 
import matplotlib.pyplot as plt
import cgs
"""
some test in disk properties with space and time 
"""
def r_Cavity(Bcp,Rcp,Mcp,dotMg):
    '''
    have been assumed to be a constant: 5.89RJ

    Bcp in Guass
    Rcp in Jupiter mass
    Mcp in Jupiter radius
    dotMg in Jupiter mass per million years
    '''
    rCavity=1.07*(Bcp/40)**(4/7)*(Rcp)**(12/7)*(Mcp/0.4)**(-1/7)*(dotMg**2/0.2)**(-2/7) *cgs.RJ
    return rCavity

def dot_M_g(t,dotMggap,tgap,tdep):
    '''
    get the gas accration rate, which is equal to the inflow of mass from CSD to CJD
    
    Parameters:
    t: simulation time from the tgap;
    dotMggap: $/dot{M}_{g,gap}$ is initial gas accration rate (constant: 0.2 MJ/Myr)
    tdep: depletion tiescale of gas (constant: 3e6 yr)
    '''
    dotMg=dotMggap*np.exp(-(t-tgap)/tdep)   #Shibaike 2019 equation (1)
    return dotMg

def Omega_K(Mcp,r):
    '''
    Keplian angular velocity

    Parameters:
    Mcp: mass of the central planet [0.4MJ-MJ]
    r: distance from central planet
    '''
    OmegaK= np.sqrt(cgs.gC*Mcp/r**3)
    return OmegaK

def Sigma_g(dotMg,OmegaK,alpha,cs):
    '''
    surface density of disk 
    
    Parameters:
    dotMg: gas accration rate
    OmegaK: Keplian angular velocity
    alpha: strength of turbulence
    cs: sound speed
    '''
    Sigmag=dotMg*OmegaK/3/np.pi/alpha/cs**2
    return Sigmag

def Tempurature(Mcp,dotMg,r,rgg,Sigmag):
    '''
    Assume that the disk is viscously heated, calculate the midplane tempurature

    Parameters:
    Mcp: mass of the central planet
    dotMg: gas accration rate
    r: distance from the central planets
    rgg: the ratio of the surface dencity of grains that affect the temparature to gas surface dencity
    '''
    Ti=100
    n=1
    while n<20:

        if 0<Ti<160:
            kapa=450*(Ti/160)**2*rgg
        elif Ti<0:
            print('Something wrong! T=',Ti,'<0')
            break
        else:
            kapa=450*rgg
    
        tau=kapa*Sigmag

        g=(3/8*tau+1/4.8/tau)**(1/4)
        Td=(3*cgs.gC*Mcp*dotMg/8/np.pi/cgs.sigmaSB/r**3)**(1/4)*g
        diff= abs(Ti-Td)
        Ti=Td
     #   print(n,'diff_Td=',diff)
        n+=1

    return Td

def Stokes_number(p,q,x,alpha,Td,Mcp,r):
    '''
    Stokes number considering radial drift without fragmentation

    Paremeters:
    p,q: exponents of Sigmag and Td
    x: the ratio between the pebble accration rate and the gas accration rate
    alpha: strenght of the turbulence
    Td: midplane tempurature [K]
    Mcp: central planetary mass [Jupiter mass]
    r: distance from central planet [Jupiter radius]
    '''
    St=0.23*(2/(3+2*p+q))**(4/5) *(10/(18-39*q))**(2/5) *(x/0.003)**(2/5) *(alpha/1e-4)**(1/5) *(Td/160)**(-2/5) *(Mcp/1)**(2/5) *(r/10)**(-2/5)
    return St

# def relative_velocity(Mcp,St,r,Td):
#     vK=np.sqrt(cgs.gC*Mcp/r**3)
#     eta=
#     vr=-2*St/(St**2+1)*eta*vK
#     return vr

dotMggap=0.2*cgs.MJ/1e6/cgs.yr2s
alpha= 1e-4
tdep= 3e6*cgs.yr2s
tgap= 1e6*cgs.yr2s

rcav= 5.89*cgs.RJ
rgg=1.7e-7

plt.figure(figsize=(12,20),dpi=192)
plt.subplot(211)
plt.title('Temperature [K]')
plt.xlabel('Distance from central planets [RJ]')
plt.plot(np.linspace(1,100,1000),np.linspace(160,160,1000),color='black',label='T=160K',linestyle='dashed')
plt.legend()
plt.subplot(212)
plt.title('Surface density [$g/cm^2$]')
plt.xlabel('Distance from central planets [RJ]')

for t in [0,3e6*cgs.yr2s,10e6*cgs.yr2s,30e6*cgs.yr2s]:
    dotMg=dot_M_g(t,dotMggap,tgap,tdep) 

    r_span=np.linspace(cgs.RJ,100*cgs.RJ,1000)  #unsure
    Mcp=0.4*cgs.MJ           #initial condition

    Sigmag_r=[]
    Td_r=[]
    for r in r_span:
        csi=3e4  #sound speed in cm/s
        ni=1
        while ni<30:
            OmegaK=Omega_K(Mcp,r)
            Sigmag=Sigma_g(dotMg,OmegaK,alpha,csi)  
            Td=Tempurature(Mcp,dotMg,r,rgg,Sigmag)
            cs=np.sqrt(cgs.kB*Td/cgs.mp/2.34)

            diff_cs=abs(cs-csi)

            csi=cs
            ni+=1
        #    print('diff_cs=',diff_cs,ni)
        Sigmag_r.append(Sigmag)
        Td_r.append(Td)

    plt.subplot(211)
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(100,3000)
    plt.plot(r_span/cgs.RJ,Td_r)
    plt.subplot(212)
    plt.xscale('log')
    plt.yscale('log') 
    plt.ylim(1,1e5)
    plt.plot(r_span/cgs.RJ,Sigmag_r)

