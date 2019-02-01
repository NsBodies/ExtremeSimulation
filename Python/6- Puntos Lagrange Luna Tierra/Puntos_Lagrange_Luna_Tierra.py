# -*- coding: utf-8 -*-

from Integrador import R4, Euler

import numpy as np
from matplotlib import pyplot as plt



MT=5.972e24  #Masa Tierra
ML=7.349e22  #Masa Luna
G=6.67e-11   #Constante Gravitatoria

t0 = 0.0                     #Instante inicial
U0 =np.array([1, 0, 0, 1])   #Condiciones iniciales
T=4500000                    #Tiempo de simulación (28 dias)
dt=60                        #Paso temporal de 1min

L_TL = 3.844e8

"Punto L1"

xL1 = -3.26e8
yL1 = 0
zL1 = 0
vxL1 = 0
vyL1 = 2*np.pi*abs(xL1)/2358176
vzL1 = 0

"Punto L2"

xL2 = -4.44e8
yL2 = 0
zL2 = 0
vxL2 = 0
vyL2 = 2*np.pi*abs(xL2)/2358176
vzL2 = 0

"Punto L3"

xL3 = 3.8647e8
yL3 = 0
zL3 = 0
vxL3 = 0
vyL3 = -2*np.pi*abs(xL3)/2358176
vzL3 = 0

"Punto L4"

xL4 = -3.84e8*np.cos(60*np.pi/180)
yL4 = 3.84e8*np.sin(60*np.pi/180)
zL4 = 0
vxL4 = 2*np.pi*abs(L_TL)/2358176*np.cos(30*np.pi/180)
vyL4 = 2*np.pi*abs(L_TL)/2358176*np.sin(30*np.pi/180)
vzL4 = 0

"Punto L5"

xL5 = -3.84e8*np.cos(60*np.pi/180) + 1e6
yL5 = -3.84e8*np.sin(60*np.pi/180)
zL5 = 0
vxL5 = -2*np.pi*abs(L_TL)/2358176*np.cos(30*np.pi/180)
vyL5 = 2*np.pi*abs(L_TL)/2358176*np.sin(30*np.pi/180)
vzL5 = 0




"Condiciones iniciales de posición y velocidad del satélite"
xsat = xL1           
ysat = yL1
zsat = zL1
vxsat = vxL1
vysat = vyL1
vzsat = vzL1


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




