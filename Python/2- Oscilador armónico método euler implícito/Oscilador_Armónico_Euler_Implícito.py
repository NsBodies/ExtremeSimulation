# -*- coding: utf-8 -*-
import numpy as np

T = 20*np.pi                                       #Tiempo de simulación
dt = 0.01                                          #Paso temporal
N = int(T/dt) + 1                                  #Número de pasos
U = np.zeros([N,2])                                #Matriz solución
U[0,0:2] = [1,0]                                   #Condiciones iniciales
 
B = np.linalg.inv(np.identity(2) - dt*np.matrix('0 1; -1 0'))
 
for i  in range(1,N):   
     U[i,:] = np.matmul(B,U[i-1,:]  )

print(U[N-1,1]**2+U[N-1,0]**2)
    
#Plot
t = np.linspace(0,T,N,endpoint = True)
from matplotlib import pyplot as plt

fig, ax = plt.subplots(figsize=(8,8))
plt.grid()
plt.xlabel("x",fontsize=18)
plt.ylabel(r'$\dot{x}$',fontsize=18)
ax.tick_params(axis='y',labelsize=14)
ax.tick_params(axis='x',labelsize=14)
plt.plot(U[:,0],U[:,1],'b')


#Save
plt.savefig('Órbita_Euler_Implícito.jpg', dpi=300)

