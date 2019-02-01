# -*- coding: utf-8 -*-

def Euler(F,t0, U0, tmax, dt): 
    
    import numpy as np
    
    
    if type(U0) == float:
    
        n = int((tmax - t0)/dt) + 1
        t = np.zeros(n)
        U = np.zeros(n)
        U[0] = U0
        t = np.linspace(t[0], tmax, n, endpoint = True)
        
        for i in range(0, n-1):     
         
           U[i+1] = U[i] + dt*F(t0, U[i]) 
           t0 = t0 + dt 
    
    else:
        s = len(U0)
        n = int((tmax - t0)/dt) + 1
        t = np.zeros(n)
        U = np.zeros((s,n))
        U[:,0] = U0
        t = np.linspace(t[0], tmax, n, endpoint = True)
    
        for i in range(0, n-1): 
            
            U[:,i+1] = U[:,i] + dt*F(t0, U[:,i]) 
            t0 = t0 + dt 
      
    return U,t


def R4(F,t0, U0, tmax, dt): 
    
    import numpy as np
    
    
    if type(U0) == float:
    
        n = int((tmax - t0)/dt) + 1
        t = np.zeros(n)
        U = np.zeros(n)
        U[0] = U0
        t = np.linspace(t[0], tmax, n, endpoint = True)
        
        for i in range(0, n-1):     
          
           k1 = dt*F(t0, U[i]) 
           k2 = dt*F(t0 + 0.5*dt, U[i] + 0.5*k1) 
           k3 = dt*F(t0 + 0.5*dt, U[i] + 0.5*k2) 
           k4 = dt*F(t0 + dt, U[i] + k3) 
           
           U[i+1] = U[i] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4) 
           t0 = t0 + dt 
    
    else:
        s = len(U0)
        n = int((tmax - t0)/dt) + 1
        t = np.zeros(n)
        U = np.zeros((s,n))
        U[:,0] = U0
        t = np.linspace(t[0], tmax, n, endpoint = True)
    
        for i in range(0, n-1): 
        
            k1 = dt*F(t0, U[:,i]) 
            k2 = dt*F(t0 + 0.5*dt, U[:,i] + 0.5*k1) 
            k3 = dt*F(t0 + 0.5*dt, U[:,i] + 0.5*k2) 
            k4 = dt*F(t0 + dt, U[:,i] + k3) 
               
            U[:,i+1] = U[:,i] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4) 
            t0 = t0 + dt 
      
    return U,t