#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:06:33 2020

@author: shawvey
"""

# steady-state output firing rate depend on the input firing rate
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
    count=0
    for i in range(0, timeSteps-1):
        for r in range(0, N):
            S[r]=S[r]-S[r]*delta_t/Tau_s
            randnum = rand.uniform(0,1)
            if (randnum < fire_rate * delta_t): # pre_s
                S[r]+=delta_s
                time_diff = post_st - i*delta_t
                if(time_diff > 0):
                    g_bar[r] = g_bar[r] + A_plus * math.exp(-1*abs(time_diff)/Tau_plus)
                else:
                    g_bar[r] = g_bar[r] - A_minus * math.exp(-1*abs(time_diff)/Tau_minus)
                if(g_bar[r] < 0):
                    g_bar[r] = 0
                pre_st[r] = i*delta_t
                
        if(Vol[i] > Vth):  # post
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
    return count/30


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

#### Question 3

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

#fire_rate = 15*math.pow(10, -3)
input_fire_rate = np.arange(10,20+2,2)
print(input_fire_rate)

### flag = on
print('------ flag = on ------')
flag = 'on'
output_fire_rate_on=[]
for i in range(0,len(input_fire_rate)):
    fire_rate = input_fire_rate[i]*math.pow(10, -3)
    outputfr= Simulation_on(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)
    output_fire_rate_on.append(outputfr)
print(output_fire_rate_on)
### flag = off
plt.plot(input_fire_rate, output_fire_rate_on,color='indianred')
plt.xlabel("Input Firing Rate(Hz)")
plt.ylabel("Output Firing Rate(Hz)")
plt.title("The mean output firing rate as a function of the input firing rates when STDP=on")
plt.savefig('PartB_Q3(1)_STDP_on')
plt.show()


flag = 'off'
print('------ flag = off ------')
output_fire_rate_off=[]
for i in range(0,len(input_fire_rate)):
    fire_rate = input_fire_rate[i]*math.pow(10, -3)
    outputfr= Simulation_off(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, fire_rate, N, A_plus, A_minus, Tau_plus, Tau_minus, flag)
    output_fire_rate_off.append(outputfr)
print(output_fire_rate_off)


plt.plot(input_fire_rate, output_fire_rate_off,color='steelblue')
plt.xlabel("Input Firing Rate(Hz)")
plt.ylabel("Output Firing Rate(Hz)")
plt.title("The mean output firing rate as a function of the input firing rates when STDP=off")
plt.savefig('PartB_Q3(1)_STDP_off')
plt.show()