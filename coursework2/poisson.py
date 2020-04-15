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
    spike_count = int(big_t/width) * [0]
    for i in range(len(spike_train)):
        j = int(spike_train[i]/width)
        spike_count[j] += 1
    Fano_factor = np.var(spike_count)/np.mean(spike_count)
    return Fano_factor


def Cal_Cov(spike_train):
    inter_spike_interval = (len(spike_train)-1) * [0]
    for i in range(len(inter_spike_interval)):
        inter_spike_interval[i] = spike_train[i+1] - spike_train[i]
    Cov = np.std(inter_spike_interval)/np.mean(inter_spike_interval)
    return Cov

def Cal_autocorrelation(rho1, rho2, width, interval):
    timestep = int(width/interval)
    ac = [0]*timestep
    for i in range(timestep):
	    for j in range(len(rho1)):
		    if (rho1[j] == 1 and j > timestep):
			    ac[i] += rho2[j-(timestep-i)]
	    ac[i] = ac[i]/sum(rho1)
    return ac

def Cal_sta(stimulus, rho, width, interval):
    timestep = int(width/interval)
    sta = [0]*timestep
    for i in range(timestep):
	    for j in range(len(rho)):
		    if (rho[j] == 1 and j > timestep):
			    sta[i] += stimulus[j-(timestep-i)]
	    sta[i] = sta[i]/sum(rho)
    return sta

def Cal_sta_notadjacent(stimulus,rho,width,interval):
    timestep = int(width/interval)
    sta = zeros(timestep)
    spike_T = np.nonzero(rho)[0]
    spike_times=[]
    m=0
    s = zeros(len(stimulus))
    while m < len(spike_T):
        s=spike_T[m]+interval
        if rho[int(s)]!=0:
            spike_times.append(spike_T[m])
        m+=1
    num = len(spike_times)
    for i in range(0, timestep):
        x=0
        window=[]
        for j in range(0, num):
            if spike_times[j] - i < 0:
                x += 1
            window.append(stimulus[spike_times[j] - i])
        sta[i] = sum(window) / (num - x)
    return sta

def Cal_sta_adjacent(stimulus,rho,width,interval):
    posi_delta = int(width/interval)
    sta = [0]*50
    for i in range(50):
        for j in range(len(rho)):
            if(rho[j] == 1 and rho[j+posi_delta] == 1 and sum(rho[j+1:j+posi_delta-1]) == 0 and j > 50):
                sta[i] += stimulus[j-(50-i)]
        sta[i] = sta[i]/sum(rho)
    return sta
    
Hz=1.0
sec=1.0
ms=0.001

### question one
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
tau_ref=5*ms
spike_train = get_spike_train(rate,big_t,tau_ref)
print("Fano factor(10ms):",Cal_Fanofactor(spike_train,big_t,10*ms))
print("Fano factor(50ms):",Cal_Fanofactor(spike_train,big_t,50*ms))
print("Fano factor(100ms):",Cal_Fanofactor(spike_train,big_t,100*ms))
print("Coefficient of variation:",Cal_Cov(spike_train))
print("")


### question 2
print("---Question Two---")

big_t=20*60*sec
rho=load_data("rho.dat",int)
spike_train = []
for i in range(len(rho)):
	if (rho[i] == 1):
		spike_train.append(i*2*ms)

print("Fano factor(10ms):",Cal_Fanofactor(spike_train,big_t,10*ms))
print("Fano factor(50ms):",Cal_Fanofactor(spike_train,big_t,50*ms))
print("Fano factor(100ms):",Cal_Fanofactor(spike_train,big_t,100*ms))
print("Coefficient of variation:",Cal_Cov(spike_train))
print("")


### question 3 待解决
print("---Question Three---")
rho = load_data("rho.dat",int)
ac = Cal_autocorrelation(rho, rho, 200*ms, 2*ms)
time = np.linspace(-98, 100, 100,  endpoint=True)
#print(len(time))
#print(time)
#plt.plot(time, ac, color='green')
#plt.title('Spike triggered average over a 100 ms window')
#plt.xlabel('Time /ms')
#plt.ylabel('Spike Triggered Average')
#plt.savefig('sta.png')
#plt.show()




### question 4
print("---Question Four---")
stimulus = load_data("stim.dat",float)
rho = load_data("rho.dat",int)
sta = Cal_sta(stimulus, rho, 100*ms, 2*ms)
time = np.linspace(0, 98, 50,  endpoint=True)
plt.plot(time, sta, color='purple')
plt.title('Spike triggered average over a 100 ms window')
plt.xlabel('Time /ms')
plt.ylabel('Spike Triggered Average')
plt.savefig('sta.png')
plt.show()


### COMSM2127
print("---COMSM2127---")
sta1 = Cal_sta_notadjacent(stimulus, rho, 2*ms, 2*ms)
print(sta1)
print(len(sta1))
sta2 = Cal_sta_notadjacent(stimulus, rho, 10*ms, 2*ms)
sta3 = Cal_sta_notadjacent(stimulus, rho, 20*ms, 2*ms)
sta4 = Cal_sta_notadjacent(stimulus, rho, 50*ms, 2*ms)
time = np.linspace(0, 98, 50,  endpoint=True)
plt.plot(time, sta1,color='red',label='2ms')
plt.plot(time, sta2,color='green',label='10ms')
plt.plot(time, sta3,color='black',label='20ms')
plt.plot(time, sta4,color='yellow',label='50ms')
plt.legend()
plt.xlabel('Time (ms)')
plt.ylabel('Sta')
plt.title('spikes are not neccesarily adjacent')
plt.savefig('not neccesarily adjacent')
plt.show()

