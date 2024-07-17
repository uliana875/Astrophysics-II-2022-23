import numpy as np
import matplotlib.pyplot as plt
from numpy import pi 
#%%
Ms, Rs, Rhos, DRho, DM = [], [], [], [], []

d = 20
G = 6.67e-8
c = 3e10
rhoo = 1e15*(0.5 + d/70)
rhoc = rhoo

def proces(m0,rho0,r0,dm_dr0,drho_dr0):
    while rho0>rhoo:
        m1 = m0 + dm_dr0*d_r
        rho1 = rho0 + drho_dr0*d_r
        r1 = r0 + d_r
        
        dm_dr1 = 4*pi*r1**2*rho1
        drho_dr1 = - ((G*m1)/(c**2*r1**2))*(4*rho1 - rhoo)*(1 + (4*pi*r1**3*(rho1-rhoo))/(3*m1))*(1-(2*G*m1)/(c**2*r1))**(-1)
        
        m0 = m1
        rho0 = rho1
        dm_dr0 = dm_dr1
        drho_dr0 = drho_dr1
        r0 = r1

    return m1, rho1, r1, dm_dr1, drho_dr1
    
while rhoc*1e-16<5.5:
    
    if rhoc*1e-16<0.28:
        rhoc += 2e13        
    else:
        rhoc += 5e14

    Rhos.append(rhoc)
    
    m0 = 0
    rho0 = rhoc
    r0 = 1e-10
    dm_dr0 = 4*pi*r0**2*rho0
    drho_dr0 = - ((4*pi*G*r0**2)/(3*c**2))*(4*rho0**2-5*rho0*rhoo+rhoo**2)
    d_r = 1000
    m1, rho1, r1, dm_dr1, drho_dr1 = proces(m0,rho0,r0,dm_dr0,drho_dr0)

    DRho.append(drho_dr1)
    DM.append(dm_dr1)
    Ms.append(m1)
    Rs.append(r1)

    
#%%
plt.plot(np.array(Rhos)*1e-16, np.array(Rs)*1e-5,'o', markersize=1, color='indigo')
plt.ylabel(r'M [$M_{\odot}$]')
plt.xlabel(r'$\rho_c$'+' '+r'$[10^{-16}g/cm^3]$')
#plt.savefig('/home/dell/Astrofizyka II/R_rhoc.png',dpi=150)
plt.show()

#%%
plt.plot(np.array(Rs)*1e-5,np.array(Ms)*0.5*1e-33, 'o', markersize=1, color='indigo')
plt.ylabel(r'M [$M_{\odot}$]')
plt.xlabel('R [km]')
# plt.savefig('/home/dell/Astrofizyka II/M_R.png',dpi=150)
plt.show()

#%%
maxMs, Ds, maxRhos, maxRs = [], [], [],[]


G = 6.67e-8
c = 3e10


d=0
while d<120:
    rhoo = 1e15*(0.5 + d/70)
    rhoc = rhoo
    
    Ms, Rs, Rhos = [],[], []
    
    while rhoc*1e-16<15:
        
        if rhoc*1e-16<0.28:
            rhoc += 2e13        
        else:
            rhoc += 1e14
            
        Rhos.append(rhoc)

        m0 = 0
        rho0 = rhoc
        r0 = 1e-10
        dm_dr0 = 4*pi*r0**2*rho0
        drho_dr0 = - ((4*pi*G*r0**2)/(3*c**2))*(4*rho0**2-5*rho0*rhoo+rhoo**2)
        d_r = 1000
        m1, rho1, r1, dm_dr1, drho_dr1 = proces(m0,rho0,r0,dm_dr0,drho_dr0)
        
        Ms.append(m1)
        Rs.append(r1)
        

    plt.plot(np.array(Rs)*1e-5,np.array(Ms)*0.5*1e-33, markersize=1,label='d='+str(d))

    maxMs.append(max(Ms))
    ind = Ms.index(max(Ms))
    maxRhos.append(Rhos[ind])
    maxRs.append(Rs[ind])
    Ds.append(d)
    d += 10
    
plt.legend()
plt.ylabel(r'M [$M_{\odot}$]')
plt.xlabel('R [km]')    
#plt.savefig('/home/dell/Astrofizyka II/.png',dpi=150)
plt.show()
#%%
plt.plot(Ds,np.array(maxMs)*0.5*1e-33,'o',markersize=2,color='indigo')
plt.ylabel(r'$M_{max}$ [$M_{\odot}$]')
plt.xlabel('d')
#plt.savefig('/home/dell/Astrofizyka II/maxM_d.png',dpi=150)
plt.show()
#%%
plt.plot(Ds,np.array(maxRhos)*1e-16,'o',markersize=2,color='indigo')
plt.ylabel(r'$\rho_c$'+' '+r'$[10^{-16}g/cm^3]$')
plt.xlabel('d')
#plt.savefig('/home/dell/Astrofizyka II/maxrho_d.png',dpi=150)
plt.show()
#%%
plt.plot(Ds,np.array(maxRs)*1e-5,'o',markersize=2,color='indigo')
plt.ylabel('R [km]')
plt.xlabel('d')
#plt.savefig('/home/dell/Astrofizyka II/maxR_d.png',dpi=150)
plt.show()