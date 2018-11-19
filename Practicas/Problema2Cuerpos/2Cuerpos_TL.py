# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt

MT=5.972e24  #Masa Tierra
ML=7.349e22  #Masa Luna
G=6.67e-11   #Constante Gravitatoria

N=10000   #Número de pasos de integración
T=2500000 #28 dias
dt=T/N    

U=np.zeros((N,8))
rTL=np.zeros(2)
rLT=np.zeros(2)
aLT=np.zeros(2)
aTL=np.zeros(2)
U[0,:]=[4.789e6,0,-3.7961e8,0,0,-12.67,0,1030]

def F(U, G, ML, MT,i):
    rTL = U[i,2:4]-U[i,0:2]   
    aTL = G*ML/(np.linalg.norm(rTL)**3)*rTL
    aLT = -aTL*MT/ML
    return np.concatenate((aTL,aLT), axis=None)
    
#Leap Frog Method
for i in range (N-1):
    aux = U[i,4:8] + 0.5*dt*F(U,G, ML, MT,i)[0:4]
    U[i+1,0:4] = U[i,0:4] + dt*aux
    U[i+1,4:8] = aux + 0.5*dt*F(U,G,ML,MT,i+1)[0:4]

fig1,ax = plt.subplots(figsize=(8, 8))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18)
plt.plot(U[:,0],U[:,1],color = 'b')      
plt.plot(U[:,2],U[:,3],color = 'r')  

#Save
plt.savefig('Orbita_Luna_Tierra.jpg', dpi=300)