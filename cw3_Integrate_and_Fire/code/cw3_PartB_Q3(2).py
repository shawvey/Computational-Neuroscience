#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 18:11:31 2020

@author: shawvey
"""

# Plot the steady-state synaptic strength distribution.
import matplotlib.pyplot as plt
import numpy as np
import random as rand
import math
plt.style.use('ggplot')
def Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag):
    S = [0]*N
    pre_st = [0]*N
    g_bar = [g]*N
    post_st = -1500
    timeSteps = int(t/delta_t)
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    g_ = g
    for i in range(0, timeSteps-1):
        for r in range(0, N):
            S[r]=S[r]-S[r]*delta_t/Tau_s
            randnum = rand.uniform(0,1)
            if (randnum < fire_rate * delta_t): # pre_s
                S[r]+=delta_s
                time_diff = post_st - i*delta_t
                g_bar[r] = g_bar[r] - A_minus * math.exp(-1*abs(time_diff)/Tau_minus)
                if(g_bar[r] < 0):
                    g_bar[r] = 0
                pre_st[r] = i*delta_t
                
        if(Vol[i] > Vth):  # post
            Vol[i] = Vrest
            for j in range(0,N):
                time_diff = i*delta_t - pre_st[j]
                g_bar[j] = g_bar[j] + A_plus * math.exp(-1*abs(time_diff)/Tau_plus)
                if(g_bar[j] > g_):
                    g_bar[j] = g_
            post_st = i*delta_t
        gsst = 0
        for q in range(0,N):
            gsst += S[q]*g_bar[q]
        Vol[i+1] = Vol[i]+(El-Vol[i]+Rm*Ie+Rm*gsst*(Es-Vol[i]))*delta_t/Tau_m
    return g_bar
######### Part B

#### Question 3 (2)

t = 300*1000
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
A_plus = 0.2*math.pow(10, -12)
A_minus = 0.25*math.pow(10, -12)
Tau_plus = 20
Tau_minus = 20
flag='on'
fire_rate1 = 10*math.pow(10, -3)
fire_rate2 = 20*math.pow(10, -3)

data1 = Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate1, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)
data2 = Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate2, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)

plt.hist(data1, bins=8, normed=0, facecolor="firebrick", alpha = 1)
plt.xlabel("Synaptic weight (A/mV)")
plt.ylabel("Frequency")
plt.title("The histogram of the steady-state synaptic weights for <r>=10HZ ")
plt.show()

plt.hist(data2, bins=8, normed=0, facecolor="olive", alpha = 1)
plt.xlabel("Synaptic weight (A/mV)")
plt.ylabel("Frequency")
plt.title("The histogram of the steady-state synaptic weights for <r>=20HZ")
plt.show()
