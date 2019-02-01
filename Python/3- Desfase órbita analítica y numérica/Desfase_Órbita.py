# -*- coding: utf-8 -*-

from Integrador import R4, Euler

import numpy as np
from matplotlib import pyplot as plt

t0 = 0.0                     #Instante inicial
U0 =np.array([1, 0, 0, 1])   #Condiciones iniciales
T = 10*2*np.pi               #Tiempo de simulación (10 vueltas)
dt = 0.1                    #Paso temporal

mu = 1.0                     #Parámetro gravitacional 
                

def F(t,U): 
    return np.array([U[2],U[3] ,-mu/(np.linalg.norm(U[0:2])**3)*U[0],-mu/(np.linalg.norm(U[0:2])**3)*U[1]]) 
   

[[x,y,vx,vy],t] = R4(F,t0, U0, T, dt)     #Integrador Runge Kutta 4
#[[x,y,vx,vy],t] = Euler(F,t0, U0, T, dt) #Integrador Euler
 

#Plot
fig1,ax = plt.subplots(figsize=(8, 8))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18)
plt.plot(x[:],y[:])
plt.plot(x[0],y[0],'o',color = 'r')      #Solución analítca de la posición en el instante final
plt.plot(x[-1],y[-1],'o',color = 'g')    #Solución numérica de la posición en el instante final

#Save
plt.savefig('Desfase_Órbita.jpg', dpi=300)