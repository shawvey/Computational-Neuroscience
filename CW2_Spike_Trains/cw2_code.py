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


def Cal_Fanofactor_Q2(rho,width,sampled_time):
    timestep = int(width/sampled_time)
    spike_count = []
    for i in range(0, len(rho), timestep):
        count = 0
        for j in range(0, timestep):
            if rho[i+j] == 1:
                count += 1
        spike_count.append(count)
    var = np.var(spike_count)
    avg = np.mean(spike_count)
    fano_factor=var/avg
    return fano_factor


def Cal_Cov_Q2(rho,sampled_time):
    intervals = []
    spikes = np.nonzero(rho)[0] 
    for i in range(0,len(spikes)-1):
        intervals.append((spikes[i+1] - spikes[i])*sampled_time)
    sd = np.std(intervals)
    avg = np.mean(intervals)
    Cov = sd/avg
    return Cov


def Cal_corr(rho, width, sampled_time): 
    timestep = int(width/sampled_time)
    staz = np.zeros(int(timestep/2)+1)
    staf = np.zeros(int(timestep/2)+1)
    spikes = np.nonzero(rho)[0]     
    stNum = len(spikes)
    for i in range(0, int(timestep/2)+1):
        z = []
        f = []
        count=0
        for j in range(0,stNum):
            if spikes[j]-i<0:
                count+=1
            else:
                z.append(rho[spikes[j]+i])
                f.append(rho[spikes[j]-i])
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
    spikes = np.nonzero(rho)[0]     
    stNum = len(spikes)
    for i in range(0, timestep+1):
        stiValue = []
        count = 0
        for j in range(0,stNum):
            if spikes[j]-i<0:
                count += 1
            else:
                stiValue.append(stimulus[spikes[j]-i])
        sta[i] = sum(stiValue)/(stNum-count)
    sta = np.flip(sta)
    return sta


def Cal_sta_notadjacent(stimulus,rho,width,interval,sampled_time):
    timestep = int(width/sampled_time)
    interval = int(interval/sampled_time)
    sta = np.zeros(timestep+1)
    spikes = np.nonzero(rho)[0]
    newSpikes = []
    for m in range(0,len(spikes)):
        if(rho[int(spikes[m]+interval)]==1):
            flag=0
            for n in range(1,interval):
                if(rho[int(spikes[m]+n)] == 1):
                    flag=1
            if(flag==1):
                newSpikes.append(spikes[m])
    stNum = len(newSpikes)
    for i in range(0, timestep+1):
        stiValue = []
        count=0
        for j in range(0,stNum):
            if newSpikes[j]<i:
                count += 1
            stiValue.append(stimulus[newSpikes[j]-i])
        sta[i] = sum(stiValue)/(stNum-count)
    sta=np.flip(sta)
    return sta


def Cal_sta_adjacent(stimulus,rho,width,interval,sampled_time):
    timestep = int(width/sampled_time)
    interval = int(interval/sampled_time)
    sta = np.zeros(timestep+1)
    spikes = np.nonzero(rho)[0]
    newSpikes = []
    for m in range(0,len(spikes)):
        if(rho[int(spikes[m]+interval)]==1):
            flag=1
            for n in range(1,interval):
                if(rho[int(spikes[m]+n)] == 1):
                    flag=0
            if(flag==1):
                newSpikes.append(spikes[m])
    stNum = len(newSpikes)
    for i in range(0, timestep+1):
        stiValue = []
        count=0
        for j in range(0,stNum):
            if newSpikes[j]<i:
                count += 1
            stiValue.append(stimulus[newSpikes[j]-i])
        sta[i] = sum(stiValue)/(stNum-count)
    sta=np.flip(sta)
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

sampled_time = 2*ms
rho = load_data("rho.dat",int)
print("Fano factor(10ms):",Cal_Fanofactor_Q2(rho,10*ms,sampled_time))
print("Fano factor(50ms):",Cal_Fanofactor_Q2(rho,50*ms,sampled_time))
print("Fano factor(100ms):",Cal_Fanofactor_Q2(rho,100*ms,sampled_time))
print("Coefficient of variation:",Cal_Cov_Q2(rho,sampled_time))
print("")




width = 200*ms

### question 3

rho = load_data("rho.dat",int)
sta = Cal_corr(rho, width, sampled_time)
plt.style.use('ggplot')
plt.plot(np.linspace(-100, 100,101,  endpoint=True), sta, color='teal')
plt.title('Autocorrelogram over the range -100ms to +100ms')
plt.xlabel('Window(ms)')
plt.ylabel('correlation')
plt.savefig('autocorrelogram.png')
plt.show()




### question 4
width = 100*ms
print("---Question Four---")
stimulus = load_data("stim.dat",float)
rho = load_data("rho.dat",int)
sta = Cal_sta(stimulus, rho, width, sampled_time)
plt.style.use('ggplot')
plt.plot(np.linspace(-100, 0, 51,  endpoint=True), sta, color='slateblue')
plt.title('The spike triggered average over a 100ms window')
plt.xlabel('Window(ms)')
plt.ylabel('The spike triggered average')
plt.savefig('sta.png')
plt.show()

### COMSM2127

print("---COMSM2127---")
## spikes are not adjacent
width = 100*ms
stimulus = load_data("stim.dat",float)
rho = load_data("rho.dat",int)
sta2 = Cal_sta_notadjacent(stimulus, rho, width,10*ms, 2*ms)
sta3 = Cal_sta_notadjacent(stimulus, rho, width,20*ms, 2*ms)
sta4 = Cal_sta_notadjacent(stimulus, rho, width,50*ms, 2*ms)
time = np.linspace(-100, 0, 51,  endpoint=True)
plt.style.use('ggplot')

plt.plot(time, sta2,label='10 ms',color='orange')
plt.plot(time, sta3,label='20 ms',color='olivedrab')
plt.plot(time, sta4,label='50 ms',color='darkcyan')
plt.legend()
plt.title('The stimulus triggered by pairs of spikes(not adjacent)')
plt.xlabel('Window(ms)')
plt.ylabel('The averge stimulus')
plt.savefig('notadjacent.png')
plt.show()


## spikes are adjacent
sta1 = Cal_sta_adjacent(stimulus, rho, width,2*ms, 2*ms)
sta2 = Cal_sta_adjacent(stimulus, rho, width,10*ms, 2*ms)
sta3 = Cal_sta_adjacent(stimulus, rho, width,20*ms, 2*ms)
sta4 = Cal_sta_adjacent(stimulus, rho, width,50*ms, 2*ms)
time = np.linspace(-100, 0, 51,  endpoint=True)
plt.style.use('ggplot')
plt.plot(time, sta1,label='2 ms',color='tomato',)
plt.plot(time, sta2,label='10 ms',color='orange')
plt.plot(time, sta3,label='20 ms',color='olivedrab')
plt.plot(time, sta4,label='50 ms',color='darkcyan')
plt.legend()
plt.title('The stimulus triggered by pairs of spikes(adjacent)')
plt.xlabel('Window(ms)')
plt.ylabel('The averge stimulus')
plt.savefig('adjacent.png')
plt.show()

