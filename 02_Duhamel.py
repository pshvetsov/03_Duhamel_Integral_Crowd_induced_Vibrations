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
    
    #Cycle through the time vector and evaluate the response at each time point
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
            
            #Calculate the response
            U[i] = A_i*math.e**(-xi*wn*T[i])*math.sin(wd*T[i]) - \
                B_i*math.e**(-xi*wn*T[i])*math.cos(wd*T[i])
    return U

response = duhamel(time, force)

#-------------------------- Plot Duhamel response -----------------------------
Plot = False
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])
    axes.plot(time,response)
    #Housekeeping
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Disp (m)')
    axes.set_title('SDoF system response to harmonic loading')
    axes.set_xlim([0,tMax])
    plt.grid()
    plt.show()


#Validating the numerical solution
beta = wf/wd
k = m*wd**2

O = (P/k)*(1/((1-beta**2)**2 + (2*xi*beta)**2))
response_cf = O*((1-beta**2)*np.sin(wf*time) - 2*xi*beta*np.cos(wf*time))

#------------------------- Plot analytical response ---------------------------
Plot = False
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])
    axes.plot(time, response, label='Duhamel')
    axes.plot(time, response_cf, label='Analytical')
    #Housekeeping
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Disp (m)')
    axes.set_title('SDoF system response to harmonic loading')
    axes.legend(loc='lower right')
    axes.set_xlim([0,tMax])
    plt.grid()
    plt.show()
    
#------------------------ Dynamic crowd loading -------------------------------
L = 60 #[m] Bridge span
vp = 1.3 #[m/s] Pedestrian walking velocity
tMax = L/vp #[s] Crossing time
m = 80 #[kg] Pedestrian mass
G = m*9.81 #[N] Static weight of pedestrian

fv = 0.35*vp**3 - 1.59*vp**2 + 2.93*vp
DLF = 0.41*(fv-0.95)

print(f"- The DLF = {round(DLF,3)} and the pacing frequency is {round(fv,2)} Hz"\
      " ({round(fv,2)} steps per second)")
print(f"- Duration of a single step is {round(1/fv,2)} seconds")

delT = 0.005
time = np.arange(0, tMax+delT, delT)
Fv = G + abs(G*DLF*np.sin(2*math.pi*(fv/2)*time))

#---------------------- Plot ground reaction force ----------------------------
Plot = False
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])
    axes.plot(time,Fv,label='GRF')
    #Housekeeping  
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Force (N)')
    axes.set_title('Vertical ground reaction force')
    axes.legend(loc='lower right')
    axes.set_xlim([0,5])
    plt.grid()
    plt.show()

#---------------- Dynamic analysis: Bridge +1 pedestrian ----------------------
xp = vp*time
phi = np.sin(math.pi*xp/L)
Fn = Fv*phi
#----------------------------- Plot modal force -------------------------------
Plot = False
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])
    axes.plot(time, Fn, label='Modal load')
    #Housekeeping
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Force (N)')
    axes.set_title('Modal Force')
    axes.set_xlim([0,tMax])
    plt.grid()
    plt.show()
    
M = 2000 #[kg/m] Mass per unit length
m = 0.5*M*L #[kg] Modal mass of mode 1
xi = 0.025 #[-] Damping ration
fn = 2.5 #[Hz] bridge modal frequency
wn = 2*math.pi*2.5 #[rads/s] Angular modal frequency
wd = wn*math.sqrt(1-xi**2) #[rads/s] Damped angular modal frequency
response = duhamel(time, Fn) #Response calculated using Duhamel integral function

#--------------------------- Plot Modal responce ------------------------------
Plot = False
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])
    axes.plot(time, -response)
    #Housekeeping
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Disp (m)')
    axes.set_title('Modal response (static+dynamic)')
    axes.set_xlim([0,tMax])
    plt.grid()
    plt.show()

#Static component of GRF
Fn_static = G*phi 
#Dynamic component of GRF
Fn_dynamic = abs(G*DLF*np.sin(2*math.pi*(fv/2)*time))*phi

response_static = duhamel(time, Fn_static)
response_dynamic = duhamel(time, Fn_dynamic)

#----------------- Plot separated dynamic and static responses-----------------
#--------------------- as well as total combined response ---------------------
Plot = False
if Plot:
    fig, axes = plt.subplots(figsize=(14,10), nrows=2,ncols=1)
    axes[0].plot(time,-response_dynamic,label='Dynamic')
    axes[0].plot(time,-response_static,label='Static')
    axes[0].set_xlabel('time (s)')
    axes[0].set_ylabel('Disp (m)')
    axes[0].set_title('Modal response (separate components)')
    axes[0].legend(loc='lower right')
    axes[0].set_xlim([0,tMax])
    axes[0].grid()

    axes[1].plot(time,-response,label='Total')

#Function to calculate min values of oscillation
def Peaks(disp, time):
    peaks = np.empty([1,0])
    times = np.empty([1,0])
    
    #Calculate slopes for each data point
    slopes = np.zeros(len(disp))
    for i, u in enumerate(disp):
        if (i < len(disp)-1):
            slopes[i] = disp[i+1] - disp[i]
            
    #Cycle through all slopes and pick out peaks
    for i, s in enumerate(slopes):
        if (i<len(slopes)-1):
            if (slopes[i+1]<0 and slopes[i]>0):
                peaks = np.append(peaks, disp[i])
                times = np.append(times, time[i])
    return [peaks, times]

#--------------------------- Plot reponse envelope ----------------------------
peaks, times = Peaks(response, time)
Plot = True
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])
    axes.plot(time, -response, label='Response')
    axes.plot(times, -peaks, label='Response envelope')
    #Housekeeping
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Disp (m)')
    axes.set_title('Modal response (static+dynamic)')
    axes.set_xlim([0,tMax])
    plt.grid()
    plt.show()

k = m*wn**2 #[N/m] Original system stiffness
Masses = [1750, 2000, 2250] #[kg/m] Masses per unit length to test

#------------------- Plot reponse envelope for var. systems -------------------
Plot = True
if Plot:
    fig = plt.figure()
    axes = fig.add_axes([0.1,0.1,2,1])

    for M in Masses:
        m = 0.5*M*L #[kg] Modal mass of mode 1
        wn = math.sqrt(k/m)
        wd = wn*math.sqrt(1-xi**2) #[rad/s] Damped angular modal frequncy
        
        response = duhamel(time, Fn)
        peaks, times = Peaks(response, time)
        
        axes.plot(times, -peaks, label=f"Response envelope M={M}")
        
    #Housekeeping
    axes.set_xlabel('time (s)')
    axes.set_ylabel('Disp (m)')
    axes.set_title('Modal response (varying bridge mass)')
    axes.legend(loc='lower right')
    axes.set_xlim([0,tMax])
    plt.grid()
    plt.show()

















































