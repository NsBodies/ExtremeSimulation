# -*- coding: utf-8 -*-

from Integrador import R4, Euler

import numpy as np
from matplotlib import pyplot as plt


t0 = 0.0                     #Instante inicial
U0 =np.array([1,0])   #Condiciones iniciales
T=100                #Tiempo de simulación 
dt=0.01                    #Paso temporal 

epsilon0 = 1e-9      #Perturbación inicial

def F(t,U):
    v = U[1]
    x = U[0]
    return  np.array([v,-x])

[U1,t] = R4(F,t0, U0, T, dt)         #Integrador Runge Kutta 4
[U2,t] = R4(F,t0, U0, T, dt/2)
q = 4
#[U1,t] = Euler(F,t0, U0, T, dt)      #Integrador Euler
#[U2,t] = Euler(F,t0, U0, T, dt/2)
# q = 1


k = np.linalg.norm(U2[:,-1] - U1[:,-1])/(dt**q*(1-0.5**q))
E1 = dt**q*k
E2 = (dt/2)**q*k

Lyapunov_coef = np.log(E1/epsilon0)/T

t = np.linspace(0,10*T,1000,endpoint = True)
E = np.zeros(len(t))
for i in range(len(t)):
    
    E[i] = epsilon0*np.exp(Lyapunov_coef*t[i])

fig,ax = plt.subplots(figsize=(10, 10))
plt.grid()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("t",fontsize=18)
plt.ylabel("Error",fontsize=18) 
plt.title(r"$\ddot{x}+x = 0 \ ; \ x(0)=1 \ , \ \dot{x}(0)=0 \ y \ dt = 0.01$", fontsize = 20)
plt.plot(t,E,color = 'b')         


plt.savefig('Error_Tiempo.jpg', dpi=300)




