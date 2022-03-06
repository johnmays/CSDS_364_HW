import numpy as np
import matplotlib.pyplot as plt
from A3a_jkm100 import *
from A1b_jkm100 import *
purples = ["#0a0612", "#392249", "#482980", "#673ab7",
           "#7a52aa", "#9779bd", "#b59fd0", "#d3c5e3"]

def movingavg(x, lam=0.5, centered=False):
    if centered:
        y = np.copy(x.astype(float))
        for i in range(np.size(y)):
            if i <= 0:
                y[i] = (1.0-lam)*y[i] + 0
            else:
                y[i] = (1.0-lam)*y[i] + lam*y[i-1]
        delay = int(np.ceil(lam/(1-lam))) -1
        y= y[delay:]
    else:
        y = np.copy(x.astype(float))
        for i in range(np.size(y)):
            if i <= 0:
                y[i] = (1.0-lam)*y[i] + 0
            else:
                y[i] = (1.0-lam)*y[i] + lam*y[i-1]
    return y

def randprocess(N, sigma=1.0):
    x = np.zeros(N).astype(float, copy=False)
    x[0] = 0
    for i in range(1, np.size(x)):
        x[i] = x[i-1] + np.random.normal(loc=0, scale=sigma)
    return x

def plotsignalandfilter(signal, filtered_signal, title = "Original Signal with Filtered Signal"):
    # t = np.linspace(0, np.size(signal), np.size(signal))
    plt.figure(figsize=(8, 5), dpi=72)
    plt.plot(signal, c=purples[6], label="signal")
    plt.plot(filtered_signal, c=purples[3], label="filtered signal")
    plt.xlabel("Time (ms)")
    plt.ylabel("Amplitude")
    plt.title(title)
    plt.xlim(0, np.size(signal))
    plt.legend()
    plt.show()
    return None

def filterIIR(x_old, a, b):
    assert type(a) == list or type(a) == np.ndarray
    assert type(b) == list or type(b) == np.ndarray
    x = np.copy(x_old)  # x and y should now be numpy arrays no matter what
    y = np.zeros(np.size(x))
    for i in range(np.size(x_old)):
        for j in range(np.size(b)):
            if (i-j) >= 0:
                y[i] += b[j]*x[i-j]
        for k in range(np.size(a)):
            if (i-1-k) >= 0:
                y[i] -= a[k]*y[i-1-k]
    return y

def uniformnoise(N, alim = (-1,1)):
    y = np.zeros(N).astype(float)
    for i in range(np.size(y)):
        y[i] = np.random.uniform(low=alim[0], high=alim[1])
    return y

def plotfilterresponsegraph(a=[-1.702, 0.81], b=[0.063], duration=0.100, fs=2000):
    frequencies = [10, 50, 100, 250]
    noise_sigmas = [0.75, 0.50, 0.25, 0.1]
    
    fig, axs = plt.subplots(4, 4, sharex=True, sharey=True)
    fig.suptitle("Filter with a = {alis} and b = {blis}".format(alis=str(a), blis=str(b)))
    fig.set_size_inches(12, 12)
    t = np.linspace(0.0, duration, int((duration)*fs)+1)
    for i in range(4):
        f = frequencies[i]
        for j in range(4):
            sigma = noise_sigmas[j]
            s = noisysignal(t, g=lambda t: sinewave(t, f=f), tau=0.0,T=0.10,sigma=sigma)
            filtered_s = filterIIR(s, a=a, b=b)
            axs[i,j].plot(t,s, c=purples[6], label='signal')
            axs[i,j].plot(t,filtered_s, c=purples[3], label='filtered signal')
            axs[i,j].set_title("{freq}Hz sinewave with $\sigma$ = {sig} noise".format(freq=f, sig=sigma), fontdict={'fontsize': 8})
    plt.legend()
    plt.show()
    return None

def plotpowervsfrequency(f, filter_func, fs = 2000, title = "Filter Response Function"):
    t = np.linspace(0, 1, int((1-0)*fs)+1)
    nyquist = fs/2
    frequencies = np.linspace(0, nyquist, 101)
    powers = np.zeros(np.size(frequencies))
    for i in range(np.size(frequencies)):
        f = frequencies[i]
        s = sinewave(t, f)
        filtered_s = filter_func(s)
        powers[i] = power(filtered_s)   
    plt.figure(figsize=(8, 6), dpi=72)
    plt.scatter(frequencies, powers, c=purples[4])
    plt.title(title)
    plt.xlabel("Frequency(Hz)")
    plt.ylabel("Power(dB)")
    plt.xscale('log')
    plt.xlim(10,f)
    plt.show()
    return None

def impulseresponse(g):
    signal = np.zeros(1000)
    signal[0] = 1
    h = g(signal)
    return h

def plotimpulseresponse(a=[-0.9], b=[0.1], fs=2000, duration=0.1, title="Impulse Response"):
    # creating an impulse:
    signal = np.zeros(int(duration*fs))
    signal[int((fs*duration)/10)] = 1

    filtered_signal = filterIIR(signal, a=a, b=b)
    plt.axvline(9, c=purples[3], label = "impulse")
    plt.plot(filtered_signal, c=purples[6], label="response")
    plt.title(title) #fontdict={'fontsize': 8}
    plt.yticks([])
    plt.legend()
    plt.show()
    return None

def plotfourimpulseresponses():
    titles = ['LPF', 'HPF','BPF 1', 'BPF 2']
    a = [[-0.9], [0.9], [-1.265, 0.810], [-1.702, 0.81]]
    b = [[0.1], [0.1], [0.135], [0.063]]

    # creating an impuls:
    signal = np.zeros(100)
    signal[10] = 1

    fig, axs = plt.subplots(4, sharex=True)
    fig.suptitle("Impulse Responses for Systems so Far")
    fig.set_size_inches(6, 10)
    for i in range(len(titles)):
        filtered_signal = filterIIR(signal, a=a[i], b=b[i])
        # axs[i].plot(signal, c=purples[6])
        axs[i].axvline(9, c=purples[3], label = "impulse")
        axs[i].plot(filtered_signal, c=purples[6], label="response")
        axs[i].set_title(titles[i]) #fontdict={'fontsize': 8}
        # axs[i].set_ylim(-0.20, 0.20)
        axs[i].get_yaxis().set_visible(False)
    plt.legend()
    plt.show()
    return None

def convolve(x, h = [1], h0=0, klimit=None): # since this is python; I am making h0 = 0 => an entirely causal filter
    if klimit == None:
        klimit = np.size(h)
    y = np.zeros(np.size(x)).astype(float)
    for n in range(np.size(x)):
        for k in range(n+1):
            if (n-k+h0) < klimit:
                y[n] += x[k]*h[n-k+h0]
    return y