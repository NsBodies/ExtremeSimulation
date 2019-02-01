# -*- coding: utf-8 -*-

from Integrador import R4, Euler

import numpy as np
from matplotlib import pyplot as plt



MT=5.972e24  #Masa Tierra
ML=7.349e22  #Masa Luna
G=6.67e-11   #Constante Gravitatoria

t0 = 0.0                     #Instante inicial
U0 =np.array([1, 0, 0, 1])   #Condiciones iniciales
T=100000000                #Tiempo de simulación 
dt=360                     #Paso temporal 




"Condiciones iniciales de posición y velocidad del satélite"
xsat = -3.6e8                
ysat = 0
zsat =0
vxsat =0
vysat = 800
vzsat =  0


[x0T,y0T,z0T,x0L,y0L,z0L,x0Sat,y0Sat,z0Sat] = np.array([4.789e6,0,0,-3.7961e8,0,0,xsat,ysat,zsat])
[vx0T,vy0T,vz0T,vx0L,vy0L,vz0L,vx0Sat,vy0Sat,vz0Sat,] = np.array([0,-12.67,0,0,1030,0,vxsat,vysat,vzsat])

U0= np.concatenate(([x0T,y0T,z0T,x0L,y0L,z0L,x0Sat,y0Sat,z0Sat],[vx0T,vy0T,vz0T,vx0L,vy0L,vz0L,vx0Sat,vy0Sat,vz0Sat,]),axis=None)


def F(t,U):
    
    muL = G*ML
    muT = G*MT
    
    rT = np.array([U[0],U[1],U[2]])
    rL = np.array([U[3],U[4],U[5]])
    rSat = np.array([U[6],U[7],U[8]])
    rTL = rT- rL
    rSatT = rSat - rT
    rSatL = rSat - rL
    
    vT =  np.array([U[9:12]])
    vL = np.array([U[12:15]])
    vSat = np.array([U[15:18]])
    v = np.concatenate((vT,vL,vSat),axis=None)
    
    aT  = -muL/np.linalg.norm(rTL)**3*rTL
    aL = - aT*MT/ML
    aSat = -muT/np.linalg.norm(rSatT)**3*rSatT -muL/np.linalg.norm(rSatL)**3*rSatL
    a = np.concatenate((aT,aL,aSat),axis=None)
    
    return np.concatenate((v,a),axis=None)

[U,t] = R4(F,t0, U0, T, dt)         #Integrador Runge Kutta 4
#[U,t] = Euler(F,t0, U0, T, dt)      #Integrador Euler

[xT,yT,zT,xL,yL,zL,xSat,ySat,zSat] = U[0:9]


fig,ax = plt.subplots(figsize=(10, 10))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18) 


plt.plot(xT,yT,color = 'b')         
plt.plot(xSat,ySat,'o',color = 'g') 
plt.plot(xL,yL,color = 'r') 

plt.savefig('Orbita_Satelite_Luna_Tierra_PlanoXY.jpg', dpi=300)

asign = np.sign(ySat)
signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)

index = np.where(signchange == 1)[0]


xSat_ = xSat[index]
ySat_ = 0*ySat[index]

fig1,ax = plt.subplots(figsize=(10, 10))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18) 
plt.xlim([-4e8,-3.5e8])


    
plt.plot(xSat_,ySat_,'o',color = 'g') 

#Save
plt.savefig('Plano_Poincaré.jpg', dpi=300)




