# -*- coding: utf-8 -*-
"""
Created on Sun Jun 25 10:12:02 2023

@author: shvet
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from sympy.printing import pprint as pp

sym.init_printing()

#Define symbols
m = sym.Symbol('m')
w = sym.Symbol('w')
P0 = sym.Symbol('P0')
t1 = sym.Symbol('t1')
tau = sym.Symbol('tau')
t = sym.Symbol('t')

#Construct the function to integrate
f = tau * sym.sin(w*t - w*tau)

# Use SymPy to get the definite integral with respect to tau between 0 and t
defInt = sym.integrate(f, (tau,0,t))
# Simplify the integral
sym.simplify(defInt)

P0 = 1000 #[N] Max load
t1 = 10 #(sec) Rise time
delT = 0.1 #(sec) Time-step
t = np.arange(0,t1+delT,delT) #Time vector

m = 20
periodRange = [0.3, 0.4, 0.5] #Range of system periods as a propertion of rise time

#Initialise a figure to plot onto
fig = plt.figure()
axes = fig.add_axes([0.1,0.1,2,1])
for pr in periodRange:
    T = pr*t1
    wn = 2*math.pi/T
    k = m*wn**2
    u = (P0/k)*((t/t1) - ((np.sin(wn*t))/(wn*t1)))
    axes.plot(t/t1,u/(P0/k), label=f'T = {pr}t1')

#Housekeeping
axes.set_xlabel('t/t1')
axes.set_ylabel('Displacement ratio')
axes.set_title('SDoF system response to ramp loading')
axes.legend(loc='upper left')
axes.set_xlim([0,1])
plt.grid()
plt.show()