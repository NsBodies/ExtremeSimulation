# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt

MT=5.972e24  #Masa Tierra
ML=7.349e22  #Masa Luna
G=6.67e-11   #Constante Gravitatoria

N=50000 #Número de pasos de integración
T=25000000 #290 dias
dt=T/N    

U=np.zeros((N,12))
rTL=np.zeros(2)
rLT=np.zeros(2)
rSL = np.zeros(2)
aLT=np.zeros(2)
aTL=np.zeros(2)
aSL = np.zeros(2)

U[0,:]=[4.789e6,0,-3.7961e8,0,-3.5e8,0,0,-12.67,0,1030,0,500]

def F(U, G, ML, MT,i):
    rTL = U[i,2:4]-U[i,0:2]
    rST = U[i,0:2]-U[i,4:6]
    rSL = U[i,2:4]-U[i,4:6]
    aTL = G*ML/(np.linalg.norm(rTL)**3)*rTL
    aLT = -aTL*MT/ML
    aS = G*(MT/(np.linalg.norm(rST)**3)*rST+ML/(np.linalg.norm(rSL)**3)*rSL)
    return np.concatenate((aTL,aLT,aS), axis=None)
    
#Leap Frog Method
for i in range (N-1):
    aux = U[i,6:12] + 0.5*dt*F(U,G, ML, MT,i)[0:6]
    U[i+1,0:6] = U[i,0:6] + dt*aux
    U[i+1,6:12] = aux + 0.5*dt*F(U,G,ML,MT,i+1)[0:6]

fig1,ax = plt.subplots(figsize=(8, 8))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18)
plt.plot(U[:,0],U[:,1],color = 'b')      
plt.plot(U[:,2],U[:,3],color = 'r')  
plt.plot(U[:,4],U[:,5],color = 'g') 
#Save
plt.savefig('Orbita_Satelite_Luna_Tierra.jpg', dpi=300)