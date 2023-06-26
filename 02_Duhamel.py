# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 18:44:59 2023

@author: shvet
"""

import math
import numpy as np
import matplotlib.pyplot as plt

m = 1000
xi = 0.05
f = 1.5
wn = 2*math.pi*f
wd = wn*math.sqrt(1-xi**2)

tMax = 20
delT = 0.01
time = np.arange(0,tMax+delT, delT)

#[Hz] Forcing frequency
f_Force = 1
#[rad/s] Angular forcing frequency
wf = 2*math.pi*f_Force
#Force amplitude
P = 100
#Force vector
force = P*np.sin(wf*time)

def duhamel(T,F):
    U = np.zeros(len(T))
    
    # Initialise values for the cumulative sum used 
    # to calculate A and B at each time-step.
    ACum_i = 0
    BCum_i = 0
    
    #Cycle through the time vector and evaluate the responce at each time point
    for i, t in enumerate(T):
        if i > 0:
            #Calculate A[i]
            #Value if integrand at current time-step
            y_i = math.e**(xi*wn*T[i])*F[i]*math.cos(wd*T[i])
            #Value of integrand at previous time-step
            y_im1 = math.e**(xi*wn*T[i-1])*F[i-1]*math.cos(wd*T[i-1])
            #Area of current trapezoid
            Area_i = 0.5*delT*(y_i + y_im1)
            #Cumulative area from t=0 to current time
            ACum_i +=  Area_i
            #Value of A for the current time-step
            A_i = (1/(m*wd))*ACum_i
            
            #Calculate B[i]
            #Value if integrand at current time-step
            y_i = math.e**(xi*wn*T[i])*F[i]*math.sin(wd*T[i])
            #Value of integrand at previous time-step
            y_im1 = math.e**(xi*wn*T[i-1])*F[i-1]*math.sin(wd*T[i-1])
            #Area of current trapezoid
            Area_i = 0.5*delT*(y_i + y_im1)
            #Cumulative area from t=0 to current time
            BCum_i +=  Area_i
            #Value of B for the current time-step
            B_i = (1/(m*wd))*BCum_i
            
            #Calculate the responce
            U[i] = A_i*math.e**(-xi*wn*T[i])*math.sin(wd*T[i]) - \
                B_i*math.e**(-xi*wn*T[i])*math.cos(wd*T[i])
    return U

responce = duhamel(time, force)

fig = plt.figure()
axes = fig.add_axes([0.1,0.1,2,1])
axes.plot(time,responce)

#Housekeeping
axes.set_xlabel('time (s')
axes.set_ylabel('Disp (m')
axes.set_title('SDoF system responce to harmonic loading')
axes.set_xlim([0,tMax])
plt.grid()
plt.show()






