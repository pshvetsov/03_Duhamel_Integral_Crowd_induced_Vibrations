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
Plot = False
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
Plot = False
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

#----------------- Dynamic analysis bridge +N pedestrian crowd-----------------
#Number of pedestrians that coss the bridge in the time window
N = 100 
# [s] Simulation window
window = 10*60 
# [s] Additional seconds to allow simualtion of responce beyond window length (late finishers)
buffer = 200 
# [kg] Pedestrian mass
mp = 80
# [N] static weight of pedestrian
G = 9.81*mp

#Random variables
# Uniformly distributed start times
tStart = np.random.uniform(low=0.0, high=window, size=N)
tStart[0] = 0.0
# Normally distributed walking velocities
Vp = np.random.normal(loc=1.3, scale=0.125, size=N)

tMax = window + buffer
time = np.arange(0, tMax+delT, delT)

crowdForce = np.zeros([N, len(time)])
crowdResponse = np.zeros([N, len(time)])

# For each pedestrian...
for i, n in enumerate(np.arange(N)):
    vp = Vp[i] # [m/s] Walking velocity
    startTime = tStart[i] #[s] Start time
    tCross = L/vp #[s] Crossing time
    tEnd = startTime + tCross #[s] Finish time
    
    fv = 0.35*vp**2 - 1.59*vp**2 + 2.93*vp # [Hz] Pacing frequency
    DLF = 0.41*(fv-0.95)
    # Time vector for this pedestrian
    timeVector = np.arange(0, tCross+delT, delT)
    # Static + Dynamic GRF (ignore static component)
    Fv = G + abs(G*DLF*np.sin(2*math.pi*(fv/2)*timeVector))
    
    xp = vp*timeVector
    phi = np.sin(math.pi*xp/L)
    Fn = Fv*phi
    # Responce calculated using the Duhamel integral function
    response = duhamel(timeVector, Fn)
    
    # Save the GRF and response for this pedestrian at the 
    # correct position in the overal simlation records
    iStart = round(startTime/delT) # Index for start time
    crowdForce[i,iStart:iStart+len(Fn)] = Fn
    crowdResponse[i,iStart:iStart+len(Fn)] = response

#------- Plot idividual modal forces and responses for each pedestrian --------
Plot = False
if Plot:
    fig, axes = plt.subplots(figsize=[14,10], nrows=2, ncols=1)
    for i in np.arange(len(crowdForce)):
        axes[0].plot(time,crowdForce[i,:])
        axes[1].plot(time,-crowdResponse[i,:])
        
    #Housekeeping
    axes[0].plot([window,window],[0,np.max(crowdForce)],'r--')
    axes[0].plot([window+buffer,window+buffer],[0,np.max(crowdForce)],'r--')
    axes[0].set_xlabel('time (s)')
    axes[0].set_ylabel('Force (N)')
    axes[0].set_title('Individual Modal forces')
    axes[0].set_xlim([0,tMax])
    #axes[0].set_xlim([startTime,startTime+tCross])
    axes[0].grid()

    axes[1].plot([window,window],[0,-np.max(crowdResponse)],'r--')
    axes[1].plot([window+buffer,window+buffer],[0,-np.max(crowdResponse)],'r--')
    axes[1].set_xlabel('time (s)') 
    axes[1].set_ylabel('Disp (m)')
    axes[1].set_title('Individual Modal responses')
    axes[1].set_xlim([0,tMax])
    #axes[1].set_xlim([startTime,startTime+tCross])
    axes[1].grid()

#Sum across rows of crowdForce and crowdResponse
F_Crowd = sum(crowdForce)
Res_crowd = sum(crowdResponse)

peaks, times = Peaks(Res_crowd, time)

#------- Plot cummulative modal forces and responses for all pedestrian -------
Plot = False
if Plot:
    fig, axes = plt.subplots(figsize=[14,10], nrows=2, ncols=1)

    axes[0].plot(time, F_Crowd)
    axes[1].plot(time, -Res_crowd)
    axes[1].plot(times, -peaks, 'r-')
    
    axes[0].plot([window,window],[0,max(F_Crowd)],'r--')
    axes[0].plot([window+buffer,window+buffer],[0,max(F_Crowd)],'r--')
    axes[0].set_xlabel('time (s)')
    axes[0].set_ylabel('Force (N)')
    axes[0].set_title('Cummulative Modal forces')
    axes[0].set_xlim([0,tMax])
    axes[0].grid()

    axes[1].set_xlabel('time (s)') 
    axes[1].set_ylabel('Disp (m)')
    axes[1].set_title('Cummulative Modal responses')
    axes[1].set_xlim([0,tMax])
    axes[1].grid()
    
    plt.show()

# -------------- Animating bridge response & traffic flow ---------------------
from matplotlib.animation import FuncAnimation
import matplotlib.gridspec as gridspec
#%matplotlib notebook

#Animation parameters
animLength = 100 #[s]
frameRate = 12 #[-] Frames per second
plotInterval = 1/frameRate #[s] time betweeen frame plots
dataInterval = int(plotInterval/delT) #Plot moving elements every 'dataInterval-th' point
defScale = 500 #Scale factor on bridge deflection (for visibility)

fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10,5))
gs = gridspec.GridSpec(2,1,height_ratios=(1,1))
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

ax1.set_aspect('equal', adjustable='box') #Set equal scale for axes top subplot

#Set axis limits
ax1.set_xlim([0,L])
ax1.set_ylim([-3,3])
yLim2 = defScale*max(Res_crowd)
ax2.set_xlim([0,L])
ax2.set_ylim([-yLim2, yLim2])

#Housekeeping
ax1.set_title('Bridge plan view')
ax1.set_xlabel('(m)')
ax1.set_ylabel('(m)')
ax1.grid()
ax2.set_title('Bridge oscillation')
ax2.set_xlabel('(m)')
ax2.set_ylabel(f'Scaled displacement\n x{defScale} (m)')
ax2.grid()

plt.show()

#Define initial state of pedestrians in top plot
topPedList = [] # Initialise an empty list to hold markers representing pedestrians
for i in np.arange(N):
    yPos = np.random.uniform(low=-2.5, high=2.5, size=1)
    pedTop, = ax1.plot(0,yPos,'o',markeredgecolor='k', markersize=10)
    topPedList.append(pedTop)

#Define initial state of pedestrians in bottom plot
btmPedList = []
for i in np.arange(N):
    ped, = ax2.plot([0,0], [0,0.6*yLim2])
    btmPedList.append(ped)

#Define the initial state of the beam, in the bottom plot
xVals = np.arange(0,L+1,1) #An array of x-values along the beam
phiVals = np.sin(math.pi*xVals/L) #Corresponding phi-values
beamDisp = 0*phiVals #Initial array of displacements along the beam

axisLine, = ax2.plot(xVals,beamDisp, 'k')
defLine, = ax2.plot(xVals,beamDisp, 'r')

# Function to animate plot objects
def animate(i):
    frm = i*dataInterval # Index of data for this frame
    simTime = time[frm] # Simulation time for this animation frame
    
    # Update the pedestrian positions (top plot) for the current frame
    for i in np.arange(N):
        if (simTime>=tStart[i] and simTime<tStart[i] + L/Vp[i]):
            Pt = topPedList[i]
            pos = (simTime - tStart[i])*Vp[i]
            Pt.set_xdata([pos,pos])
            
    # Update the beam deflected shape for the current frame
    defLine.set_data(xVals, -defScale*phiVals*Res_crowd[frm])
    
    # Update the pedestrian positions (bottom plot) for the current frame
    for i in np.arange(N):
        if (simTime>=tStart[i] and simTime<tStart[i] + L/Vp[i]):
            Pb = btmPedList[i]
            pos = (simTime - tStart[i])*Vp[i]
            h = 0.1 + 0.1*crowdForce[i,frm]/max(crowdForce[i,:])
            Pb.set_data([pos,pos], [0,h])
            
    return Pt, Pb, defLine,

# Function to generate the animation
myAnimation = FuncAnimation(fig,
                            animate, 
                            frames = int(1+(animLength/plotInterval)),
                            interval = plotInterval*1000,
                            blit=True,
                            repeat=True)
plt.show()
myAnimation.save('Brindge_response.gif')










































