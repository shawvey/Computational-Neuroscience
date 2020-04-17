import random as rnd
import numpy as np
import load
from load import load_data
import matplotlib.pyplot as plt


def get_spike_train(rate,big_t,tau_ref):

    if 1<=rate*tau_ref:
        print("firing rate not possible given refractory period f/p")
        return []


    exp_rate=rate/(1-tau_ref*rate)

    spike_train=[]

    t=rnd.expovariate(exp_rate)

    while t< big_t:
        spike_train.append(t)
        t+=tau_ref+rnd.expovariate(exp_rate)

    return spike_train


def Cal_Fanofactor(spike_train,big_t,width):
    spikeCount = int(big_t/width) * [0]
    for i in range(len(spike_train)):
        j = int(spike_train[i]/width)
        spikeCount[j] += 1
    Fano_factor = np.var(spikeCount)/np.mean(spikeCount)
    return Fano_factor


def Cal_Cov(spike_train):
    inter_spike_interval = (len(spike_train)-1) * [0]
    for i in range(len(inter_spike_interval)):
        inter_spike_interval[i] = spike_train[i+1] - spike_train[i]
    Cov = np.std(inter_spike_interval)/np.mean(inter_spike_interval)
    return Cov


def Cal_corr(rho, width, sampled_time): 
    timestep = int(width/sampled_time)
    staz = np.zeros(int(timestep/2)+1)
    staf = np.zeros(int(timestep/2)+1)
    spikeTimes = np.nonzero(rho)[0]     
    stNum = len(spikeTimes)
    for i in range(0, int(timestep/2)+1):
        z = []
        f = []
        count=0
        for j in range(0,stNum):
            if spikeTimes[j]-i<0:
                count+=1
            else:
                z.append(rho[spikeTimes[j]+i])
                f.append(rho[spikeTimes[j]-i])
        staz[i] = sum(z)/(stNum-count)
        staf[i] = sum(f)/(stNum-count)
    staf = np.flip(staf)
    staf = list(staf)
    staz = list(staz)
    staf.extend(staz[1:])
    return staf


def Cal_sta(stimulus, rho, width, sampled_time):
    timestep = int(width/sampled_time)
    sta = np.zeros(timestep+1)
    spikeTimes = np.nonzero(rho)[0]     
    stNum = len(spikeTimes)
    for i in range(0, timestep+1):
        stiValue = []
        count = 0
        for j in range(0,stNum):
            if spikeTimes[j]-i < 0:
                count += 1
            else:
                stiValue.append(stimulus[spikeTimes[j]-i])
        sta[i] = sum(stiValue)/(stNum-count)
    sta = np.flip(sta)
    return sta


def Cal_sta_notadjacent(stimulus,rho,width,interval,sampled_time):
    timestep = int(width/sampled_time)
    interval = int(interval/2)
    sta = np.zeros(timestep+1)
    spikeTimes = np.nonzero(rho)[0]
    newSpikeTimes = []
    for m in range(0,len(spikeTimes)):
        if rho[int(spikeTimes[m]+interval)] == 1:
            newSpikeTimes.append(spikeTimes[m])
    stNum = len(newSpikeTimes)
    for i in range(0, timestep+1):
        stiValue = []
        count=0
        for j in range(0,stNum):
            if newSpikeTimes[j]<i:
                count += 1
            stiValue.append(stimulus[newSpikeTimes[j]-i])
        sta[i] = sum(stiValue)/(stNum-count)
    return sta


def Cal_sta_adjacent(stimulus,rho,width,interval,sampled_time):
    timestep = int(width/sampled_time)
    interval = int(interval/2)
    sta = np.zeros(timestep+1)
    spikeTimes = np.nonzero(rho)[0]
    newSpikeTimes = []
    for m in range(0,len(spikeTimes)):
        if rho[int(spikeTimes[m]+interval)]==1 and sum(rho[int(spikeTimes[m])+1:int(spikeTimes[m]+interval)])==0:
            newSpikeTimes.append(spikeTimes[m])
    stNum = len(newSpikeTimes)
    for i in range(0, timestep+1):
        stiValue = []
        count = 0
        for j in range(0,stNum):
            if newSpikeTimes[j]<i:
                count += 1
            stiValue.append(stimulus[newSpikeTimes[j]-i])
        sta[i] = sum(stiValue)/(stNum-count)
    return sta



Hz=1.0
sec=1.0
ms=0.001

### question 1
print("---Question One---")
print("~ With no refractory period ~")
rate=35.0 *Hz
big_t=1000*sec
spike_train = get_spike_train(rate,big_t,0)
print("Fano factor(10ms):",Cal_Fanofactor(spike_train,big_t,10*ms))
print("Fano factor(50ms):",Cal_Fanofactor(spike_train,big_t,50*ms))
print("Fano factor(100ms):",Cal_Fanofactor(spike_train,big_t,100*ms))
print("Coefficient of variation:",Cal_Cov(spike_train))
print("")
print("~ With 5ms refractory period ~")
tau_ref = 5*ms
spike_train = get_spike_train(rate,big_t,tau_ref)
print("Fano factor(10ms):",Cal_Fanofactor(spike_train,big_t,10*ms))
print("Fano factor(50ms):",Cal_Fanofactor(spike_train,big_t,50*ms))
print("Fano factor(100ms):",Cal_Fanofactor(spike_train,big_t,100*ms))
print("Coefficient of variation:",Cal_Cov(spike_train))
print("")


### question 2

print("---Question Two---")

big_t = 20*60*sec
rho = load_data("rho.dat",int)
spike_train = []
for i in range(len(rho)):
	if (rho[i] == 1):
		spike_train.append(i*2*ms)

print("Fano factor(10ms):",Cal_Fanofactor(spike_train,big_t,10*ms))
print("Fano factor(50ms):",Cal_Fanofactor(spike_train,big_t,50*ms))
print("Fano factor(100ms):",Cal_Fanofactor(spike_train,big_t,100*ms))
print("Coefficient of variation:",Cal_Cov(spike_train))
print("")


sampled_time = 2*ms
width = 200*ms

### question 3

rho = load_data("rho.dat",int)
sta = Cal_corr(rho, width, sampled_time)
plt.style.use('ggplot')
plt.plot(np.linspace(-100, 100,101,  endpoint=True), sta, color='teal')
plt.title('Autocorrelogram over the range -100ms to +100ms')
plt.xlabel('Interval(ms)')
plt.ylabel('correlation')
plt.savefig('autocorrelogram.png')
plt.show()


width = 100*ms

### question 4

print("---Question Four---")
stimulus = load_data("stim.dat",float)
rho = load_data("rho.dat",int)
sta = Cal_sta(stimulus, rho, width, sampled_time)
plt.style.use('ggplot')
plt.plot(np.linspace(0, 100, 51,  endpoint=True), sta, color='slateblue')
plt.title('The spike triggered average over a 100ms window')
plt.xlabel('Interval(ms)')
plt.ylabel('The spike triggered average')
plt.savefig('sta.png')
plt.show()


### COMSM2127

print("---COMSM2127---")
## spikes are not adjacent
sta1 = Cal_sta_notadjacent(stimulus, rho, width, 2, 2*ms)
sta2 = Cal_sta_notadjacent(stimulus, rho, width, 10, 2*ms)
sta3 = Cal_sta_notadjacent(stimulus, rho, width, 20, 2*ms)
sta4 = Cal_sta_notadjacent(stimulus, rho, width, 50, 2*ms)
time = np.linspace(0, 100, 51,  endpoint=True)
plt.style.use('ggplot')
plt.plot(time, sta1,label='2 ms',color='tomato',)
plt.plot(time, sta2,label='10 ms',color='orange')
plt.plot(time, sta3,label='20 ms',color='olivedrab')
plt.plot(time, sta4,label='50 ms',color='darkcyan')
plt.legend()
plt.title('The stimulus triggered by pairs of spikes(not adjacent)')
plt.xlabel('Interval(ms)')
plt.ylabel('The averge stimulus')
plt.savefig('notadjacent.png')
plt.show()


## spikes are adjacent
sta1 = Cal_sta_adjacent(stimulus, rho, width, 2, 2*ms)
sta2 = Cal_sta_adjacent(stimulus, rho, width, 10, 2*ms)
sta3 = Cal_sta_adjacent(stimulus, rho, width, 20, 2*ms)
sta4 = Cal_sta_adjacent(stimulus, rho, width, 50, 2*ms)
time = np.linspace(0, 100, 51,  endpoint=True)
plt.style.use('ggplot')
plt.plot(time, sta1,label='2 ms',color='tomato',)
plt.plot(time, sta2,label='10 ms',color='orange')
plt.plot(time, sta3,label='20 ms',color='olivedrab')
plt.plot(time, sta4,label='50 ms',color='darkcyan')
plt.legend()
plt.title('The stimulus triggered by pairs of spikes(adjacent)')
plt.xlabel('Interval(ms)')
plt.ylabel('The averge stimulus')
plt.savefig('adjacent.png')
plt.show()
