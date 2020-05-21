#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 15:26:19 2020

@author: shawvey
"""

import matplotlib.pyplot as plt
import math
import numpy as np
plt.style.use("ggplot") 

def Simulation(t,Vrest,El,Vth,Rm,Ie,delta_t,Tau_m):
    timeSteps = int(t/delta_t)
    Vol = [0]*(timeSteps + 1)
    Vol[0] = El
    for i in range (1, timeSteps+1):
        Vol[i] = Vol[i-1] + (El-Vol[i-1]+Rm*Ie)*delta_t/Tau_m
        if(Vol[i] > Vth):
            Vol[i] = Vrest
    return Vol

######### Part A COMSM2127

t = 1000
Tau_m = 10
El = -70
Vrest = -70
Vth = -40
Rm = math.pow(10, 10)
Ie = 3.1*math.pow(10, -9)
delta_t = 0.25

times= np.arange(0, t + delta_t, delta_t)
Vol = Simulation(t,Vrest,El,Vth,Rm,Ie,delta_t,Tau_m)
plt.plot(times, Vol, color='steelblue') 
plt.xlabel("Time (ms)")
plt.ylabel("Voltage (mv)")
plt.title("The simulation of an integrate and fire model")
plt.savefig('PartA_Q1')
plt.show()



