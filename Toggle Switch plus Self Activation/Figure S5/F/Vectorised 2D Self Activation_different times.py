# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 16:15:07 2022

@author: N.V
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import os
# note that this must be executed before 'import numba'
os.environ['NUMBA_DISABLE_INTEL_SVML'] = '1'
from numba import njit
import time as process_time
import seaborn as sns
sns.set_theme(context='notebook', style='ticks',font='arial',font_scale=1.2)
 #Degradtion rate
gamma_A=0.5; gamma_B=0.5;
#Transcription rate
g_A=10 ; g_B=10;


# @njit
def hill(X,X0,l,n):
    #H+
    H1=(X**n)/(X0**n+X**n)
    #H-
    H2=(X0**n)/(X**n+X0**n)
    #Hill function =H-(+) lambda*H+
    H=H2+l*H1
    return H







def interactions_sa(X):
    dXdt=np.zeros(np.shape(X))
    X=np.asarray(X)
   
    #Degradtion rate
    gamma_A=0.5; gamma_B=0.5;
    #Transcription rate
    g_A=10 ; g_B=10;
    #Hills function threshold 
    A0B=10 ; B0A=10;
    A0A=10;B0B=10;
    
    #Co-operativity 
    nAtoB=2;
    nBtoA=2;
    nAtoA=3;
    nBtoB=3;
    
    
    #fold change
    lambda_AtoB =0.1
    lambda_BtoA=0.1
    lambda_AtoA=3
    lambda_BtoB=3
    
    #Equation=
    dXdt[0,:,:]=g_A*hill(X[0,:,:],A0A,lambda_AtoA,nAtoA)*hill(X[1,:,:],B0A,lambda_BtoA,nBtoA) -gamma_A*X[0,:,:]
    dXdt[1,:,:]=g_B*hill(X[1,:,:],B0B,lambda_BtoB,nBtoB)*hill(X[0,:,:],A0B,lambda_AtoB,nAtoB) -gamma_B*X[1,:,:]
    
    return dXdt





# dimsension size, mm
w = h = 10.
# intervals in x-, y- directions, mm
dx = dy = 0.1
#Diffusivity
D = 0.0001
alpha=10; gamma=0.5

u_initial, v_initial = 50, 50
T_max=1500;
nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = 0.01; #dx2 * dy2 / (2 * D * (dx2 + dy2))
resolution=300/dt;
#Setting the initial conditions

u0 = u_initial * np.ones((nx, ny))+ 10*np.random.random((nx, ny))
v0 = v_initial * np.ones((nx, ny))+ 10*np.random.random((nx, ny))


u = u0.copy();v = v0.copy()






# @njit
def do_timestep(u0, u, v0, v):
    # Propagate with forward-difference in time, central-difference in space
     X0= [ u0[1:-1,1:-1] , v0[1:-1,1:-1] ]
    
     F=interactions_sa(X0)
    
     u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * ((u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2+ (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )+ dt*F[0,:,:] 
    
     v[1:-1, 1:-1] = v0[1:-1, 1:-1] + D * dt * ((v0[2:, 1:-1] - 2*v0[1:-1, 1:-1] + v0[:-2, 1:-1])/dx2+ (v0[1:-1, 2:] - 2*v0[1:-1, 1:-1] + v0[1:-1, :-2])/dy2 ) +dt*F[1,:,:]
    
    ## Boundary Condition along x axis:
     u[0,:],u[-1,:]=u[1,:],u[-2,:]
     v[0,:],v[-1,:]=v[1,:],v[-2,:]
    
    ## Boundary Condition along y axis:
     u[:,0],u[:,-1]=u[:,1],u[:,-2]
     v[:,0],v[:,-1]=v[:,1],v[:,-2]
    
    
     u0 = u.copy();v0 = v.copy()
     return u0, u, v0, v
 
    
def plot_heatmap(u_k, k, species):
    # Clear the current plot figure
    f, ax = plt.subplots(figsize=(5,4)) 
    cm=ax.pcolormesh(u_k/(g_A/gamma_A), cmap=plt.cm.jet, vmin=0, vmax=0.8 ,shading='gouraud')
    clb=f.colorbar(cm,ax=ax)
    
    ax.set_aspect("equal")  
    ax.set_xlabel("x",weight='bold')
    ax.set_ylabel("y",weight='bold')
     
    #clb=f.colorbar(sc,ax=ax,shrink=0.8)
    #clb.ax.set_title(r"$\arctan(\Delta P)$ ")
    f.savefig(species+"_at_t= "+str(k)+".png",dpi=800)
    plt.close()
    

#%%


# Number of timesteps
nsteps = int(T_max/dt)

t_start=process_time.time()
uk = np.empty((nx, ny));vk=u.copy()
for m in range(nsteps):
    u0, u, v0, v = do_timestep(u0, u, v0, v)
    
    if (m%resolution)==0:
        print("At time = ",m*dt)
        vk=np.dstack((vk,v));
        uk=np.dstack((uk,u));
        plot_heatmap(u,m*dt,"A");
    
uk=np.moveaxis(uk, -1, 0)
vk=np.moveaxis(vk, -1, 0)
print("Simulation time (in min)= ",(process_time.time()-t_start)/60)


#%%



# for i in range(2,np.shape(uk)[0]):
#      plt.clf()
#      plt.rcParams.update({'font.size':14})
#     # plt.title(f"A at t = {k:.3f} min",weight='bold')
#      plt.xlabel("x",weight='bold')
#      plt.ylabel("y",weight='bold')
#      u_k=uk[i,:,:]
#      plt.pcolormesh(u_k/(g_A/gamma_A), cmap=plt.cm.jet, vmin=0, vmax=0.8 ,shading='gouraud')
#      plt.colorbar()
#      plt.show()
