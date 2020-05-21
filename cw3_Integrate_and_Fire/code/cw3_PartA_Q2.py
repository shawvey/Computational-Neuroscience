#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 12 17:50:19 2020

@author: shawvey
"""

import matplotlib.pyplot as plt
import numpy as np
import random as rand


def Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, RmIe, Rmgs, Es, delta_s, delta_t):
    timeSteps = int(t/delta_t) + 1
    firstVol = [0]*timeSteps
    secondVol= [0]*timeSteps
    S1 = [0]*timeSteps
    S2 = [0]*timeSteps
    firstVol[0] = rand.randint(Vrest, Vth)
    secondVol[0] = rand.randint(Vrest, Vth)
    for i in range(1, timeSteps):
        S1[i] = S1[i-1]-S1[i-1]*delta_t/Tau_s
        S2[i] = S2[i-1]-S2[i-1]*delta_t/Tau_s
        firstVol[i] = firstVol[i-1]+(El-firstVol[i-1]+RmIe+Rmgs*S1[i-1]*(Es-firstVol[i-1]))*delta_t/Tau_m
        secondVol[i] = secondVol[i-1]+(El-secondVol[i-1]+RmIe+Rmgs*S2[i-1]*(Es-secondVol[i-1]))*delta_t/Tau_m
        if (firstVol[i]>Vth):
            firstVol[i] = Vrest
            S2[i] += delta_s
        if (secondVol[i]>Vth):
            secondVol[i] = Vrest
            S1[i] += delta_s
    return firstVol,secondVol


######### Part A
    
#### Question 2

t = 1000
Tau_m = 20
El = -70
Vrest = -80
Vth = -54
RmIe = 18
Rmgs = 0.15
delta_s = 0.5
Tau_s = 10
Es = -80
delta_t = 0.25


##### The synapses are inhibitory with Es = -80 mv
times= np.arange(0, t + delta_t, delta_t)
firstVol,secondVol = Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, RmIe, Rmgs, Es, delta_s, delta_t)
plt.plot(times, firstVol, 'indianred', label = 'Neuron A')
plt.plot(times, secondVol, 'steelblue', label = 'Neuron B')
plt.legend(loc='best', bbox_to_anchor=(1.01,1.0))
plt.title('The synapses are inhibitory with Es = -80 mv')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mv)')
plt.savefig('partA_Q2_inhibitory_-80mv')
plt.show()

#### The synapses are excitatory with Es = 0 mv
Es = 0
firstVol,secondVol = Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, RmIe, Rmgs, Es, delta_s, delta_t)
plt.plot(times, firstVol, 'indianred', label = 'Neuron A')
plt.plot(times, secondVol, 'steelblue', label = 'Neuron B')
plt.legend(loc='best', bbox_to_anchor=(1.01,1.0))
plt.title('The synapses are excitatory with Es = 0 mv')
plt.xlabel('Time (ms)')
plt.ylabel('Voltage (mv)')
plt.savefig('partA_Q2_excitatory_0mv')
plt.show()


