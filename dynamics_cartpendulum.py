#!/usr/bin/env python
# coding: utf-8

# In[60]:


import  numpy as np
import math as mt

mc = 1
mp = 0.3
l  = 0.5 
g  = 9.806
pi = 3.14159
dt = 0.001





def dynamics(x_seed,u_seed):
    
    x = float(x_seed[0,0])
    theta = float(x_seed[0,1])
    x_dot = float(x_seed[0,2])
    theta_dot = float(x_seed[0,3])
    u = u_seed


    H = np.array([[mc+mp ,mp*l*mt.cos(theta)],
                       [mp*l*mt.cos(theta), mp*(l**2)]],np.float32)

    C = np.array([[0 , -mp*l*theta_dot*mt.sin(theta)],
                      [0, 0]],np.float32)

    G = np.array([[0,mp*g*l*mt.sin(theta)]],np.float32)

    G = G.transpose()

    B = np.array([[1],[0]],np.float32)

    inv_H = np.linalg.inv(H)

    qdot = np.array([[x_dot],[theta_dot]])
    q2dot = np.matmul(inv_H,(np.multiply(B,u)-np.matmul(C,qdot)-G))

    X_dot = np.concatenate((qdot,q2dot))

    return(X_dot)

    
    
# x_seed = np.array([[-0.00716441, -2.9711456 , -0.01667254 , 0.2032204 ]],np.float32)
# u_seed = np.array([[-0.46402387523690786]],np.float32)   

#X = dynamics(x_seed,u_seed)

    
    
    



