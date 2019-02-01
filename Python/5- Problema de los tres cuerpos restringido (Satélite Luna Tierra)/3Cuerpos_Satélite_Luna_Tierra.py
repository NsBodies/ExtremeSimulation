# -*- coding: utf-8 -*-

from Integrador import R4, Euler

import numpy as np
from matplotlib import pyplot as plt



MT=5.972e24  #Masa Tierra
ML=7.349e22  #Masa Luna
G=6.67e-11   #Constante Gravitatoria

t0 = 0.0                     #Instante inicial
U0 =np.array([1, 0, 0, 1])   #Condiciones iniciales
T=2500000                    #Tiempo de simulación (28 dias)
dt=60                      #Paso temporal de 1min  




"Condiciones iniciales de posición y velocidad del satélite"
xsat = -3.6e8                
ysat = 1e6
zsat =-1e6
vxsat = 300
vysat = 700
vzsat =  -100


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


N  = int(T/dt/1000)  #Se pintan 1000 puntos


xT = xT[0:-1:N]
yT = yT[0:-1:N]
zT = zT[0:-1:N]

xL = xL[0:-1:N]
yL = yL[0:-1:N]
zL = zL[0:-1:N]

xSat = xSat[0:-1:N]
ySat = ySat[0:-1:N]
zSat = zSat[0:-1:N]


fig1,ax = plt.subplots(figsize=(10, 10))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18)
plt.plot(xT,yT,color = 'b')      
plt.plot(xL,yL,color = 'r')  
plt.plot(xSat,ySat,color = 'g') 

#Save
plt.savefig('Orbita_Satelite_Luna_Tierra_PlanoXY.jpg', dpi=300)


fig = plt.figure()
plt.figure(figsize=(10,10))
ax = fig.gca(projection='3d')
ax.set_xlim(-0.45e9,0.45e9)
ax.set_ylim(-0.45e9,0.45e9)
ax.set_zlim(-7e6, 7e6)

ax.plot(xT,yT,zT, color = 'b')
ax.plot(xL,yL,zL, color = 'r')
ax.plot(xSat,ySat,zSat, color = 'g')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

fig.savefig('Orbita_Satelite_Luna_Tierra_3D.jpg', dpi=300)


