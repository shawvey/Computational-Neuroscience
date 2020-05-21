#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 21:33:25 2020

@author: shawvey
"""
# Plot the cross-correlogram for both cases


import matplotlib.pyplot as plt
import numpy as np
import random as rand
import math
plt.style.use('ggplot')


def Get_spikeTrains(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, N, A_plus, A_minus, Tau_plus, Tau_minus, flag, r_zero, f, B):
    S = [0]*N
    pre_st = [0]*N
    g_bar = [g]*N
    post_st = -1500
    timeSteps = int(t/delta_t)
    Vol=[0]*timeSteps
    Vol[0] = Vrest
    g_ = g
    newtimestep = int((t/delta_t)/4)
    spike_train_post = [0]*newtimestep
    spike_train_pre = np.zeros([N,newtimestep])
    No=0
    for i in range(0, timeSteps-1, 4):
        flag_post=0
        flag_pre = [0]*N
        for j in range(1,4):
            for r in range(0, N):           
                S[r]=S[r]-S[r]*delta_t/Tau_s
                randnum = rand.uniform(0,1)
                fire_rate = r_zero + B * math.sin(2*math.pi*f*(i+j)*delta_t)
                if (randnum < fire_rate * delta_t): # pre_s
                    flag_pre[r] = 1
                    S[r]+=delta_s
                    time_diff = post_st - (i+j)*delta_t
                    g_bar[r] = g_bar[r] - A_minus * math.exp(-1*abs(time_diff)/Tau_minus)
                    if(g_bar[r] < 0):
                        g_bar[r] = 0
                    pre_st[r] = (i+j)*delta_t
            
            if(Vol[i+j] > Vth):  # post
                flag_post=1
                Vol[i+j] = Vrest
                for s in range(0,N):
                    time_diff = (i+j)*delta_t - pre_st[s]
                    g_bar[s] = g_bar[s] + A_plus * math.exp(-1*abs(time_diff)/Tau_plus)
                    if(g_bar[s] > g_):
                        g_bar[s] = g_
                post_st = (i+j)*delta_t   

            gsst = 0
            for q in range(0,N):
                gsst += S[q]*g_bar[q]
            if(i+j+1<timeSteps):
                Vol[i+j+1] = Vol[i+j]+(El-Vol[i+j]+Rm*Ie+Rm*gsst*(Es-Vol[i+j]))*delta_t/Tau_m
                
        for r in range(0,N):
            spike_train_pre[r][No]=flag_pre[r]
            
        spike_train_post[No]=flag_post
        
        No+=1
    return spike_train_post,spike_train_pre


def first20seconds(spike_train_post,spike_train_pre,N):
    corr = [0]*101
    timesteps = int(20/(1*math.pow(10, -3)))
    spikes_post = np.nonzero(spike_train_post[0:timesteps])[0]     
    stNum = len(spikes_post)
    for i in range(-50,51,1):
        count_a = [0]*40
        for j in range(0,timesteps):
            for r in range(0,N):
                if(j+i>=0 and spike_train_post[j]==1):
                    count_a[r]+=spike_train_pre[r][j+i]
        corr_each = [0]*40
        for r in range(0,N):
            corr_each[r]= count_a[r]/(stNum)
        corr[i]=np.mean(corr_each)
    print(corr)
    return corr

def last100seconds(spike_train_post,spike_train_pre,N):
    corr = [0]*101
    timesteps_start = int(200/(1*math.pow(10, -3)))
    timesteps_end = int(300/(1*math.pow(10, -3)))
    spikes_post = np.nonzero(spike_train_post[timesteps_start:timesteps_end])[0]
    stNum = len(spikes_post)
    for i in range(-50,51,1):
        count_a = [0]*40
        for j in range(timesteps_start,timesteps_end):
            for r in range(0,N):
                if(j+i<timesteps_end and spike_train_post[j]==1):
                    count_a[r]+=spike_train_pre[r][j+i]
        corr_each = [0]*40
        for r in range(0,N):
            corr_each[r]= count_a[r]/(stNum)
        corr[i]=np.mean(corr_each)
    print(corr)
    return corr


######### Part B

#### COMSM2127

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
r_zero = 15*math.pow(10, -3)
f = 10*math.pow(10, -3)
B = 0
interval = np.linspace(-50, 50, 101)
spike_train_post,spike_train_pre = Get_spikeTrains(t, Tau_m, Tau_s, El, Vrest, Vth, Rm, Ie, g , Es, delta_s, delta_t, N, A_plus, A_minus, Tau_plus, Tau_minus, flag, r_zero, f, B)
corrfirst20 = first20seconds(spike_train_post,spike_train_pre,N)
corrlast100 = last100seconds(spike_train_post,spike_train_pre,N)
plt.plot(interval, corrfirst20,color='indianred',label = 'First 20s')
plt.plot(interval, corrlast100,color='steelblue',label = 'Last 100s')


plt.legend(loc='best', bbox_to_anchor=(1.01,1.0))
plt.xlabel("Interval (ms)")
plt.ylabel("Cross-Correlation")
plt.title("The cross-correlogram for first 20s and last 100s")
plt.savefig('PartB_COMSM2127')
plt.show()


