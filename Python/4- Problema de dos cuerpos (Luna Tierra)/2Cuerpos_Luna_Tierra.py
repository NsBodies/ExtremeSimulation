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
dt=3600                      #Paso temporal de 1h  


U0=[4.789e6,0,-3.7961e8,0,0,-12.67,0,1030]

def F(t,U):
    
    muL = G*ML
    
    rT = np.array([U[0],U[1]])
    rL = np.array([U[2],U[3]])
    
    rTL = rT- rL
    
    vT =  np.array([U[4:6]])
    vL = np.array([U[6:8]])
    v = np.concatenate((vT,vL),axis=None)
    
    aT  = np.array([-muL/(np.linalg.norm(rTL)**3)*rTL[0],-muL/(np.linalg.norm(rTL)**3)*rTL[1]])
    aL = - aT*MT/ML
    a = np.concatenate((aT,aL))
    
    return np.concatenate((v,a),axis=None)

[[xT,yT,xL,yL,vxT,vyT,vxL,vyL],t] = R4(F,t0, U0, T, dt)       #Integrador Runge Kutta 4
#[[xT,yT,xL,yL,vxT,vyT,vxL,vyL],t] = Euler(F,t0, U0, T, dt)   #Integrador Euler

#Plot

fig1,ax = plt.subplots(figsize=(8, 8))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18)
plt.plot(xT[:],yT[:],color = 'b')      
plt.plot(xL[:],yL[:],color = 'r')  


#Save
plt.savefig('Órbita_Luna_Tierra.jpg', dpi=300)