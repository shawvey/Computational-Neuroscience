    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 19:45:01 2020

@author: shawvey
"""
#  plot the average firing rate of the postsynaptic neuron over time
import matplotlib.pyplot as plt
import numpy as np
import random as rand
import math

def Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag):
    S = [0]*N
    pre_st = [0]*N
    g_bar = [g]*N
    post_st = -1000
    timeSteps = int(t/delta_t)
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    g_ = g
    count=0
    spikes = []
    spikes.append(0)
    
    for i in range(0, timeSteps-1):
            
        for r in range(0, N):
            S[r]=S[r]-S[r]*delta_t/Tau_s
            randnum = rand.uniform(0,1)
            if (randnum < fire_rate * delta_t and flag == 'on'): # pre_s
                S[r]+=delta_s
                time_diff = post_st - i*delta_t
                if(time_diff <= 0):
                    g_bar[r] = g_bar[r] - A_minus * math.exp(-1*abs(time_diff)/Tau_minus)
                if(g_bar[r] < 0):
                    g_bar[r] = 0
                pre_st[r] = i*delta_t
                
        if(Vol[i] > Vth and flag == 'on'):  # post
            count+=1
            Vol[i] = Vrest
            for j in range(0,N):
                time_diff = i*delta_t - pre_st[j]
                if(time_diff > 0):
                    g_bar[j] = g_bar[j] + A_plus * math.exp(-1*abs(time_diff)/Tau_plus)
                if(g_bar[j] > g_):
                    g_bar[j] = g_
            post_st = i*delta_t
            
        gsst = 0
        
        for q in range(0,N):
            gsst += S[q]*g_bar[q]
            
        Vol[i+1] = Vol[i]+(El-Vol[i]+Rm*Ie+Rm*gsst*(Es-Vol[i]))*delta_t/Tau_m
        if((i+1)%40000 == 0):
            spikes.append(count/10)
            count=0
    spikes.append(count/10)
    return spikes

######### Part B

#### Question 2

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
fire_rate = 15*math.pow(10, -3)
A_plus = 0.2*math.pow(10, -12)
A_minus = 0.25*math.pow(10, -12)
Tau_plus = 20
Tau_minus = 20
flag = 'on'
spikes=Simulation(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)
times = np.arange(0,300+10,10)
plt.plot(times, spikes,color='olivedrab')
plt.xlabel("Time(s)")
plt.ylabel("Firing Rate(Hz)")
plt.title("The average firing rate of the postsynaptic neuron over time")
plt.savefig('PartB_Q2(3)')
plt.show()



