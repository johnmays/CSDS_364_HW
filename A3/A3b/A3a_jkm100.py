# Some important imports
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wavfile

# importing A1b_jkm100.py purely for function default kwargs
import sys
sys.path.insert(1, '/Users/johnmays/Documents/GitHub/CSDS_364_HW/A1')
from A1b_jkm100 import *

purples = ["#0a0612", "#392249", "#482980", "#673ab7",
           "#7a52aa", "#9779bd", "#b59fd0", "#d3c5e3"]

def plot_sampled_function(g=sinewave, fs=1, f=1.0, tlim=(0, 2*math.pi), tscale=1, tunits="secs"):
    # all values are given in seconds; the time scale will alter the plot to show x-axis values on the order of the tscale
    t = np.linspace(tlim[0], tlim[1], 2001)
    signal = g(t=t, f=f)
    sampled_t = np.linspace(tlim[0], tlim[1], int(((tlim[1]-tlim[0])*fs)+1))
    sampled_signal = g(t=sampled_t, f=f)
    plt.figure(figsize=(8, 5), dpi=72)
    plt.plot(t*tscale, signal, c=purples[6], label="signal")
    markers, stems, base = plt.stem(sampled_t*tscale, sampled_signal, linefmt=purples[2], basefmt=" ", label="sampled signal")
    plt.setp(markers, 'color', purples[2])
    plt.setp(markers, 'marker', 'o')
    plt.xlim(tlim[0]*tscale, tlim[1]*tscale)
    plt.xlabel("Time ({units})".format(units=tunits))
    plt.ylabel("Amplitude")
    plt.title("Sampled Signal with Original Signal")
    plt.legend()
    plt.show()
    return None

def delta(t, fs=1):
    center_time_val = 0
    if type(t) == np.ndarray or type(t) == list:
        signal = []
        for time_val in t:
            if (time_val >= center_time_val - 1/(2*fs)) and (time_val < center_time_val + 1/(2*fs)):
                print("oh" + str(time_val))
                signal.append(1)
            else:
                signal.append(0)
        return signal
    else:  # t must be a number then
        if (t >= center_time_val - 1/(2*fs)) and (t <= center_time_val + 1/(2*fs)):
            return 1
        else:
            return 0

def u(t):
    if type(t) == np.ndarray or type(t) == list:
        signal = []
        for time_val in t:
            if time_val >= 0:
                signal.append(1)
            else:
                signal.append(0)
        return signal
    else:  # t must be a number then
        if t >= 0:
            return 1
        else:
            return 0

def gensignal(t, g, tau=0.025, T=0.1):
    if type(t) == list or type(t) == np.ndarray:
        signal = g(t-tau)
        for i in range(len(signal)):
            if t[i] < tau or t[i] >= tau + T:
                signal[i] = 0
        return np.array(signal)
    else:
        if t < tau or t >= tau + T:
            return 0
        else:
            return g(t-tau)


def energy(x):
    return np.linalg.norm(x, ord=2)**2

def power(x):
    return energy(x)/len(x)

def snr(Ps, Pn):
    return 10*np.log10(Ps/Pn)

def noisysignal(t, g, tau=0.025, T=0.1, sigma=1):
    if type(t) == np.ndarray or type(t) == list: # for t as an array
        signal = []
        for t_val in t:
            noise_val = np.random.normal(loc=0, scale=sigma)
            if t_val < tau or t_val >= tau + T:
                signal.append(noise_val)
            else:
                signal.append(noise_val + g(t_val-tau))
        return np.array(signal)
    else: # for single time values t
        noise_val = np.random.normal(loc=0, scale=sigma)
        if t < tau or t >= tau + T:
            return noise_val
        else:
            # CHECK
            if t-tau == 0:
                print("t-tau is = 0 at t=" + str(t))
            return noise_val + g(t-tau)

def plotsignal(g=lambda t: noisysignal(1, g = lambda t: u(t)), tlim=(0, 2*math.pi), tscale=1, tunits="secs", title="Signal", show_sampled=False, plotfs=1000, fs=10):
    # g should be an anonymized function
    t = np.linspace(tlim[0], tlim[1], int((tlim[1]-tlim[0])*plotfs+1))
    signal = []
    for t_val in t:
        signal.append(g(t_val))
    plt.figure(figsize=(8, 5), dpi=72)
    plt.plot(t*tscale, signal, c=purples[6], label="Signal with Noise")
    plt.xlim(tlim[0]*tscale, tlim[1]*tscale)
    plt.xlabel("Time ({units})".format(units=tunits))
    plt.ylabel("Amplitude")
    plt.title(title)
    if show_sampled:
        sampled_t = np.linspace(tlim[0], tlim[1], int(((tlim[1]-tlim[0])*fs)+1))
        sampled_signal = []
        for sampled_t_val in sampled_t:
            sampled_signal.append(g(sampled_t_val))
        markers, stems, base = plt.stem(sampled_t*tscale, sampled_signal, linefmt=purples[2], basefmt=" ", label="sampled signal")
        plt.setp(markers, 'color', purples[2])
        plt.setp(markers, 'marker', 'o')
        plt.legend()
    plt.show()
    return None

def snr2sigma(x, xrange=None, snr=10):
    if xrange == None:
        xrange = (0, len(x))
    Ps = power(x[xrange[0]:xrange[1]])
    sigma_s = math.sqrt(Ps)
    sigma_n = sigma_s/(10**(snr/20))
    return sigma_n

def extent(y, theta=0.01):
    max_abs_value = 0
    for y_val in y:
        if np.abs(y_val) > max_abs_value:
            max_abs_value = np.abs(y_val)
    threshold_value = theta*max_abs_value
    first_thresh_index = None
    last_thresh_index = None
    for i in range(len(y)):
        if np.abs(y[i]) > threshold_value:
            if first_thresh_index == None:
                first_thresh_index = i
            last_thresh_index = i
    return first_thresh_index, last_thresh_index    

def compose_sound(filename = "testfile.wav", fmin=10, fmax=10000, sigma=0.1, numgammatones=0):
    fs = 44100 # Hz
    T = 3 # seconds
    t = np.linspace(0, 3, 3*fs+1)

    delays = np.zeros(numgammatones)
    for i in range(numgammatones):
        delays[i] = np.random.uniform(0, T)
    
    frequencies = np.zeros(numgammatones)
    for i in range(numgammatones):
        frequencies[i] = np.random.uniform(fmin, fmax)

    signal = np.zeros(np.size(t))
    for i in range(np.size(t)):
        for j in range(numgammatones):
            signal[i] += gensignal(t[i], g = lambda x: gammatone(x, f=frequencies[j], normalize=False), tau=delays[j], T=T)

    # Normalize the gammatones alone
    signal = signal/np.max(signal)

    # Add Noise
    for i in range(np.size(t)):
        signal[i] += np.random.normal(loc=0, scale=sigma)

    # Normalize again
    signal = signal/np.max(signal)

    # Write to wave file
    wavfile.write(filename, fs, signal)
    return None