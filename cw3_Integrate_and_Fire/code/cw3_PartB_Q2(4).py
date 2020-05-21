#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 14 20:58:35 2020

@author: shawvey
"""

# calculate the estimate firing rate for 'on' and 'off'
import matplotlib.pyplot as plt
import numpy as np
import random as rand
import math

def Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag):
    S = [0]*N
    pre_st = [0]*N
    g_bar = [g]*N
    post_st = -1500
    timeSteps = int(t/delta_t)
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    g_ = g
    count=0
    for i in range(0, timeSteps-1):
        for r in range(0, N):
            S[r]=S[r]-S[r]*delta_t/Tau_s
            randnum = rand.uniform(0,1)
            if (randnum < fire_rate * delta_t and flag == 'on'): # pre_s
                S[r]+=delta_s
                time_diff = post_st - i*delta_t
                if(time_diff > 0):
                    g_bar[r] = g_bar[r] + A_plus * math.exp(-1*abs(time_diff)/Tau_plus)
                else:
                    g_bar[r] = g_bar[r] - A_minus * math.exp(-1*abs(time_diff)/Tau_minus)
                if(g_bar[r] < 0):
                    g_bar[r] = 0
                pre_st[r] = i*delta_t
                
        if(Vol[i] > Vth and flag == 'on'):  # post
            if(i>1079999):
                count+=1
            Vol[i] = Vrest
            for j in range(0,N):
                time_diff = i*delta_t - pre_st[j]
                if(time_diff > 0):
                    g_bar[j] = g_bar[j] + A_plus * math.exp(-1*abs(time_diff)/Tau_plus)
                else:
                    g_bar[j] = g_bar[j] - A_minus * math.exp(-1*abs(time_diff)/Tau_minus)
                if(g_bar[j] > g_):
                    g_bar[j] = g_
            post_st = i*delta_t
            
        gsst = 0
        
        for q in range(0,N):
            gsst += S[q]*g_bar[q]
            
        Vol[i+1] = Vol[i]+(El-Vol[i]+Rm*Ie+Rm*gsst*(Es-Vol[i]))*delta_t/Tau_m
    g_mean = np.mean(g_bar)
    return count/30,g_mean

def Simulation_off(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag):
    S = [0]*N
    timeSteps = int(t/delta_t)
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    count=0
    for i in range(0, timeSteps-1):
        for r in range(0, N):
            S[r]=S[r]-S[r]*delta_t/Tau_s
            randnum = rand.uniform(0,1)
            if (randnum < fire_rate * delta_t): # pre_s
                S[r]+=delta_s
        if(Vol[i] > Vth):  # post
            if(i>1079999): # last 30seconds
                count+=1
            Vol[i] = Vrest
        sums = np.sum(S)
        Vol[i+1] = Vol[i]+(El-Vol[i]+Rm*Ie+Rm*g*sums*(Es-Vol[i]))*delta_t/Tau_m
    return count/30


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

### flag = on
flag = 'on'
print('------ flag = on ------')
fire_rate_on = 0
g_off = 0

for x in range(0,4):
    fire_rate1, g_mean = Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)
    fire_rate_on += fire_rate1
    g_off += g_mean
print(fire_rate_on/4)
print(g_off/4)
g_off = g_off/4
### flag = off
flag = 'off'
print('------ flag = off ------')
fire_rate2 = 0
for x in range(0,4):
    fire_rate2 += Simulation_off(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g_off , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)
print(fire_rate2/4)
