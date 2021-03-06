# -*- coding: utf-8 -*-
import numpy as np

T = 100                                            #Tiempo de simulación
dt = 0.01                                          #Paso temporal
N = int(T/dt) + 1                                  #Número de pasos
U = np.zeros([N,2])                                #Matriz solución
U[0,0:2] = [1,0]                                   #Condiciones iniciales

for i  in range(1,N):   
     U[i,:] = U[i-1,:] + dt*np.array([U[i-1,1],-U[i-1,0]])   #Euler explícito
    
print(U[N-1,1]**2+U[N-1,0]**2)

#Plot

from matplotlib import pyplot as plt

t = np.linspace(0,T,N,endpoint = True)
fig, ax = plt.subplots(figsize=(8,8))
plt.grid()
plt.xlabel("x",fontsize=18)
plt.ylabel(r'$\dot{x}$',fontsize=18)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='x',labelsize=14)
plt.plot(U[:,0],U[:,1],'b')


#Save
plt.savefig('Órbita_Euler_Explícito.jpg', dpi=300)


