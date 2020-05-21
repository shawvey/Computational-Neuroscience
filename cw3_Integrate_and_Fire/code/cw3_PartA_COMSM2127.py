#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 20:34:49 2020

@author: shawvey
"""

import matplotlib.pyplot as plt
import math
import numpy as np
plt.style.use("ggplot") 
from matplotlib import rcParams
rcParams['font.sans-serif']=['Tahoma']

def Simulation(t,Vrest,El,Vth,Rm,Ie,delta_t,Tau_m):
    timeSteps = int(t/delta_t)
    Vol = [0]*(timeSteps + 1)
    Vol[0] = El
    for i in range (1, timeSteps+1):
        Vol[i] = Vol[i-1] + (El-Vol[i-1]+Rm*Ie)*delta_t/Tau_m
        if(Vol[i] > Vth):
            Vol[i] = Vrest
    return Vol



def firing_rate(t,Vrest,El,Vth,Rm,Ie,delta_t,Tau_m):
    timeSteps = int(t/delta_t)
    Vol = [0]*(timeSteps + 1)
    Vol[0] = El
    count = 0
    for i in range (1, timeSteps+1):
        Vol[i] = Vol[i-1] + (El-Vol[i-1]+Rm*Ie)*delta_t/Tau_m
        if(Vol[i] > Vth):
            Vol[i] = Vrest
            count+=1
    return count


######### COMS2127
    
#### Question 2

t = 1000
Tau_m = 10
El = -70
Vrest = -70
Vth = -40
Rm = math.pow(10, 10)
delta_t = 0.25
Ie = 2.9*math.pow(10, -9)
times= np.arange(0, t + delta_t, delta_t)
Vol = Simulation(t,Vrest,El,Vth,Rm,Ie,delta_t,Tau_m)
plt.plot(times, Vol, color='steelblue') 
plt.xlabel("Time (ms)")
plt.ylabel("Voltage (mv)")
plt.title("The simulation of an integrate and fire model with Ie = 2.9nA")
plt.savefig('PartA_COMS2127_2')
plt.show()

#### Question 3

Ies = np.arange(2,5,0.1) * math.pow(10, -9)
spikes = []
for i in Ies:
    Ie = i
    spike = firing_rate(t,Vrest,El,Vth,Rm,Ie,delta_t,Tau_m)
    spikes.append(spike)

plt.plot(Ies, spikes,color='olivedrab')
plt.xlabel("Current (A)")
plt.ylabel("Firing Rate (Hz)")
plt.title("Firing rate for different input currents")
plt.savefig('PartA_COMS2127_3')
plt.show()