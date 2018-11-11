# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

T = 10*2*np.pi                                #Tiempo de simulación (10 vueltas)
N = 1000                                      #Número de pasos de integración
dt = T/N
t = np.linspace(0,T,num=N,endpoint=True)
                                              #mu = 1
U = np.zeros([N,4])                           #Matriz solución
U[0,:] =np.matrix('1 0 0 1')                  #Condiciones iniciales

def F(U,i): 
    return np.array([U[i,2:4], -1/(np.linalg.norm(U[i,0:2])**3)*U[i,0:2]]) 
   
#Leap Frog Method
for i in range(N-1):
    aux = U[i,2:4] + 0.5*dt*F(U,i)[1]
    U[i+1,0:2] = U[i,0:2] + dt*aux
    U[i+1,2:4] = aux + 0.5*dt*F(U,i+1)[1]

r = np.sqrt(U[:,0]**2+U[:,1]**2)  
v = np.sqrt(U[:,2]**2+U[:,3]**2)  

#Plot
fig1,ax = plt.subplots(figsize=(8, 8))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("x",fontsize=18)
plt.ylabel("y",fontsize=18)
plt.plot(U[:,0],U[:,1])
plt.plot(U[0,0],U[0,1],'o',color = 'r')      #Solución analítca de la posición en el instante final
plt.plot(U[N-1,0],U[N-1,1],'o',color = 'g')  #Solución numérica de la posición en el instante final

#Save
plt.savefig('Orbita.jpg', dpi=300)