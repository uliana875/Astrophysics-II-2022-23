import numpy as np
from numpy import sqrt,pi,sin,cos
import matplotlib.pyplot as plt
from scipy.integrate import quad
#%%
def d_ksi(r,b):
    return r**(-2)*(b**(-2) -r**(-2)*(1 - rs/r))**(-0.5)

def get_b(rmin):
    return rmin*sqrt(rmin/(rmin-rs))
#%%

# Images for big rmin (light rays are just bent)

rs = 1
rmin = 5.7*rs
d_obs = 200.

phis = np.linspace(0,2*pi,500)

z_obs_plus, y_obs_plus = [],[]
z_obs_minus, y_obs_minus = [],[]
y1=-1

print('Wartości rmin:')
while y1<=2:
    print(rmin)
    
    b = get_b(rmin)
    
    r0 = d_obs
    del_r = 0.01
    
    ksi_min = quad(d_ksi, rmin, np.inf, args=b)[0]
    xs,ys = [],[]
    angles = []
    x = 0
    
    while x>=-20:
        
        r1 = r0 - del_r
        
        if ((r1<rmin) and (r1>-rmin)):
            r0 = r1
            continue
        
        if r1>0:
            angle = quad(d_ksi, abs(r1), np.inf, args=b)[0]
            
        else:
            angle = 2*ksi_min-quad(d_ksi, abs(r1), np.inf, args=b)[0]
            
        
        x = abs(r1)*cos(angle)
        y = abs(r1)*sin(angle)
              
        if (x>d_obs-0.5) or (x<-20+5*del_r):
            xs.append(x)
            ys.append(y)
            
        r0 = r1
        

    x0, x1 = xs[0],xs[-1]
    y0, y1 = ys[0],ys[-1]
    
    

    if (y1>=-2) and (y1<=2):
        
        R = abs(y1)
        z_kolo = R*cos(phis)
        y_kolo = R*sin(phis)
        
        
        for (z,y,phi) in zip(z_kolo,y_kolo,phis):
            if ((z>-1.5) and (z<-0.5) and (y>-1.5) and (y<1.5)) or ((z>0.5) and (z<1.5) and (y>-1.5) and (y<1.5)) or ((z>-0.5) and (z<0.5) and (y>-1.5) and (y<-0.5)):
                
                if y1>0:
                    z_obs = y0*cos(phi)
                    y_obs = y0*sin(phi)
                    z_obs_plus.append(z_obs)
                    y_obs_plus.append(y_obs)
                else:
                    z_obs = y0*cos(phi+pi)
                    y_obs = y0*sin(phi+pi)
                    z_obs_minus.append(z_obs)
                    y_obs_minus.append(y_obs)
        
        
    rmin += 0.05


########################################################################

# Images for smaller rmin (light rays make a loop around a black hole)

rmin = 1.54*rs


z_obs_plus2, y_obs_plus2 = [],[]
z_obs_minus2, y_obs_minus2 = [],[]

y1=-1
while y1<=2:
    print(rmin)
    
    b = get_b(rmin)
    
   
    
    r0 = d_obs
    del_r = 0.01
    
    ksi_min = quad(d_ksi, rmin, np.inf, args=b)[0]
    xs,ys = [],[]
    angles = []
    x = 0
    
    while x>=-20:
        
        r1 = r0 - del_r
        
        if ((r1<rmin) and (r1>-rmin)):
            r0 = r1
            continue
        
        if r1>0:
            angle = quad(d_ksi, abs(r1), np.inf, args=b)[0]
            
        else:
            angle = 2*ksi_min-quad(d_ksi, abs(r1), np.inf, args=b)[0]
            
        
        x = abs(r1)*cos(angle)
        y = abs(r1)*sin(angle)
              
        if (x>d_obs-0.5) or (x<-20+5*del_r):
            xs.append(x)
            ys.append(y)
            
        r0 = r1
        

    x0, x1 = xs[0],xs[-1]
    y0, y1 = ys[0],ys[-1]
    
    z_obs_list, y_obs_list = [],[]

    if (y1>=-2) and (y1<=2):
        
        R = abs(y1)
        z_kolo = R*cos(phis)
        y_kolo = R*sin(phis)
        
        
        for (z,y,phi) in zip(z_kolo,y_kolo,phis):
            if ((z>-1.5) and (z<-0.5) and (y>-1.5) and (y<1.5)) or ((z>0.5) and (z<1.5) and (y>-1.5) and (y<1.5)) or ((z>-0.5) and (z<0.5) and (y>-1.5) and (y<-0.5)):
                
                if y1>0:
                    z_obs = y0*cos(phi)
                    y_obs = y0*sin(phi)
                    z_obs_plus2.append(z_obs)
                    y_obs_plus2.append(y_obs)
                else:
                    z_obs = y0*cos(phi+pi)
                    y_obs = y0*sin(phi+pi)
                    z_obs_minus2.append(z_obs)
                    y_obs_minus2.append(y_obs)
    
    rmin += 0.0005
    
print('Koniec')


# Plots

fig, ax = plt.subplots(figsize=(5,5))

ax.plot(z_obs_plus,y_obs_plus,'.',label=r'$6.4<r_{min}<7.4$',color='indigo')
ax.plot(z_obs_minus,y_obs_minus,'.',label=r'$5.7<r_{min}<6.4$',color='indigo')
ax.plot(z_obs_plus2,y_obs_plus2,'.',label=r'$1.542<r_{min}<1.545$',color='indigo')
ax.plot(z_obs_minus2,y_obs_minus2,'.',label=r'$1.540<r_{min}<1.542$',color='indigo')


# Making a letter
dx, dy = 0.05, 0.05
size = 2*rs
x = np.arange(-size,size, dx)
y = np.arange(-size,size, dy)

X,Y = np.meshgrid(x,y)

x00, y00 = [],[]
for i in range(len(X)):
    for (x,y) in zip(X[i],Y[i]):
        if ((x>-1.5) and (x<-0.5) and (y>-1.5) and (y<1.5)) or ((x>0.5) and (x<1.5) and (y>-1.5) and (y<1.5)) or ((x>-1) and (x<1) and (y>-1.5) and (y<-0.5)):
            x00.append(x)
            y00.append(y)

#ax.plot(x00,y00,'.',color='black') # uncomment to plot the original letter


#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.grid()
#plt.savefig('4.png',dpi=150)
plt.show()
#%%

# Trajectories of the light rays in the vicinity of a black hole

rs = 1
rmin = 1.543*rs

b = get_b(rmin)

r0 = 30
del_r = 0.001

ksi_min = quad(d_ksi, rmin, np.inf, args=b)[0]
xs,ys = [],[]
angles = []
x = 0

for i in range(50000):
    r1 = r0 - del_r
    
    if ((r1<rmin) and (r1>-rmin)):
        r0 = r1
        continue
    
    if r1>0:
        angle = quad(d_ksi, abs(r1), np.inf, args=b)[0]
    else:
        angle = 2*ksi_min-quad(d_ksi, abs(r1), np.inf, args=b)[0]
    
    x = abs(r1)*cos(angle)
    xs.append(abs(r1)*cos(angle))
    ys.append(abs(r1)*sin(angle))
        
    r0 = r1
    angles.append(angle)
   
    

fig, ax = plt.subplots(figsize=(5,5))
ax.plot(xs,ys,'-',color='darkorange')
ax.plot([-19.9,-19.9],[-2,2],color='red', lw=4)
circle = plt.Circle((0,0),rs,color='black')
ax.add_patch(circle)
ax.set_xlim(-20,20)
ax.set_ylim(-20,20)
ax.set_title(r'$r_{min}=$'+str(rmin))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.grid()
#plt.savefig('/home/dell/coding_subjects/Astrofizyka II/Problem 3/'+str(rmin)+'.png',dpi=150)
plt.show()