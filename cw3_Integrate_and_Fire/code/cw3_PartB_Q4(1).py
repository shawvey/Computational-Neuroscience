#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 19:32:09 2020

@author: shawvey
"""
# Plot the mean and standard deviation affect the steady-state synaptic weights
import matplotlib.pyplot as plt
import numpy as np
import random as rand
import math
plt.style.use('ggplot')
def Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, N, A_plus, A_minus, Tau_plus, Tau_minus, flag, r_zero, f, B):
    S = [0]*N
    pre_st = [0]*N
    g_bar = [g]*N
    post_st = -1500
    timeSteps = int(t/delta_t)
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    g_ = g
    for i in range(0, timeSteps-1):
        fire_rate = r_zero + B * math.sin(2*math.pi*f*i*delta_t)
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
    mean = np.mean(g_bar)
    sd = np.std(g_bar,ddof = 1)
    return mean,sd
######### Part B

#### Question 4 (1)

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
r_zero = 20*math.pow(10, -3)
f = 10*math.pow(10, -3)
Bs=[0,5,10,15,20]
means = []
sds = []
for i in range(0,len(Bs)):
    B = Bs[i]*math.pow(10, -3)
    mean,sd = Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, N, A_plus, A_minus, Tau_plus, 
                  Tau_minus, flag, r_zero, f, B)
    means.append(mean)
    sds.append(sd)

x_ticks = np.linspace(0, 20, 5)
plt.xticks(x_ticks)
plt.plot(Bs, means, color='indianred')
plt.xlabel("B(Hz)")
plt.ylabel("Mean value")
plt.title("The mean of the steady-state synaptic strengths as a function of B",y=1.1)
plt.savefig('PartB_Q4_mean')
plt.show()

x_ticks = np.linspace(0, 20, 5)
plt.xticks(x_ticks)
plt.plot(Bs, sds, color='steelblue')
plt.xlabel("B(Hz)")
plt.ylabel("Standard deviation")
plt.title("The standard deviation of the steady-state synaptic strengths as a function of B",y=1.1)
plt.savefig('PartB_Q4_standard_deviation')
plt.show()
