#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 10:22:25 2020

@author: shawvey
"""

import matplotlib.pyplot as plt
import numpy as np
import random as rand
import math

def Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N):
    S = [0]*N
    timeSteps = int(t/delta_t)+1
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    for i in range(1, timeSteps):
        for r in range(0, N):
            S[r]=S[r]-S[r]*delta_t/Tau_s
            randnum = rand.uniform(0,1)
            if (randnum < fire_rate * delta_t):
                S[r]+=delta_s
        sums = np.sum(S)
        Vol[i] = Vol[i-1]+(El-Vol[i-1]+Rm*Ie+Rm*g*sums*(Es-Vol[i-1]))*delta_t/Tau_m
        if (Vol[i] > Vth):
            Vol[i] = Vrest
    return Vol


######### Part B
    
#### Question 1

t = 1000
Tau_m = 10
El = -65
Vrest = -65
Vth = -50
Rm = math.pow(10, 11)  
Ie = 0
N = 40
g = 4*math.pow(10, -12)
delta_s = 0.5
Tau_s = 2
Es = 0
delta_t = 0.25
fire_rate = 15*math.pow(10, -3)


times= np.arange(0, t+delta_t, delta_t)
Vol= Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N)
plt.plot(times, Vol, 'steelblue')
plt.legend(loc='best')
plt.title('The simulation of model with conductance = 4nS and firing rate = 15Hz')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mv)')
plt.savefig('partB_Q1')
plt.show()


