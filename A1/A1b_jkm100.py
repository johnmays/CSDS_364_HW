import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import wave  # for importing .wav files as iterables

# 1a


def sinewave(t, f=1.0, d=0.0):
    if(type(t) == list):
        t = np.array(t)
    phi = 2 * np.pi * d * f
    return np.sin(2 * np.pi * f * t + phi)

def cosinewave(t, f=1.0, d=0.0):
    if(type(t) == list):
        t = np.array(t)
    phi = 2 * np.pi * d * f
    return np.cos(2 * np.pi * f * t + phi)


def plot_sinewave(t, f=4.5, d=1.0):
    figure(figsize=(8, 6), dpi=80)
    sine_vals = sinewave(t, f=4.5, d=1.0)
    plt.plot(t, sine_vals)
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.show()


def plot_delayed_sinewave(t, f=1.0, d=0.0):
    figure(figsize=(8, 6), dpi=80)
    sine_vals = sinewave(t, f=5, d=0.0)
    sine_vals_delayed = sinewave(t, f=5, d=d)
    plt.plot(t, sine_vals)
    plt.plot(t, sine_vals_delayed, '--', c="cornflowerblue")
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.show()

# 1d


def gabor(t, a=1.0, sigma=1, f=1.0, phi=0.0):
    return a * np.exp((-t**2)/(2.0 * (sigma**2))) * np.cos(2*np.pi*f*t+phi)


def gabore(t, a=1.0, sigma=1, f=1.0):
    return gabor(t, a, sigma, f=f, phi=0.0)


def gaboro(t, a=1.0, sigma=1, f=1.0):
    return gabor(t, a, sigma, f=f, phi=(np.pi/2))


def gabor_norm(fs, sigma=1, f=1.0, phi=0.0):
    vanish_point = math.sqrt(-math.log(0.01)*2*(sigma)**2)
    time_val = -vanish_point
    t = []
    while time_val < vanish_point:
        t.append(time_val)
        time_val += 1/fs

    gabor_vals = []
    for t_val in t:
        gabor_vals.append(gabor(t_val, a=1.0, sigma=sigma, f=f, phi=phi))
    return np.linalg.norm(gabor_vals, ord=2)


def gabore_norm(fs, sigma=1, f=1.0):
    vanish_point = math.sqrt(-math.log(0.01)*2*(sigma)**2)
    time_val = -vanish_point
    t = []
    while time_val < vanish_point:
        t.append(time_val)
        time_val += 1/fs

    gabor_vals = []
    for t_val in t:
        gabor_vals.append(gabore(t_val, a=1.0, sigma=sigma, f=f))
    return np.linalg.norm(gabor_vals, ord=2)


def gaboro_norm(fs, sigma=1, f=1.0):
    vanish_point = math.sqrt(-math.log(0.01)*2*(sigma)**2)
    time_val = -vanish_point
    t = []
    while time_val < vanish_point:
        t.append(time_val)
        time_val += 1/fs

    gabor_vals = []
    for t_val in t:
        gabor_vals.append(gaboro(t_val, a=1.0, sigma=sigma, f=f))
    return np.linalg.norm(gabor_vals, ord=2)


def plot_gabor(t, sigma=4.0, f=1.0, a=1.0):
    gabore_vals = []
    for t_val in t:
        gabore_vals.append(gabore(t_val, sigma=sigma, f=f, a=a))

    figure(figsize=(8, 6), dpi=80)
    plt.plot(t, gabore_vals)
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.show()

# 1c


def gammatone(t, n: int = 4, f=1.0, phi=0.0, normalize = False):
    if(type(t) == list):
        t = np.array(t)
    b = 1.019*(24.7*(((4.37*f)/1000) + 1))
    gamma_value = (t**(n-1))*np.exp(-2*np.pi*b*t)*np.cos(2*np.pi*f*t + phi)
    if normalize:
        return gamma_value / gammatone_norm(n=n, f=f)
    else: 
        return gamma_value

def gammatone_norm(n: int = 4, f=1.0):
    return np.max(np.abs(gammatone(np.linspace(0, 3, 44100*3+1), n=n, f=f, normalize=False)))


def plot_gammatone(t, f=1.0, xlim=(0, 0.1)):
    gammatone_vals = gammatone(t, f=f)
    figure(figsize=(8, 6), dpi=80)
    plt.plot(t, gammatone_vals)
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.xlim(xlim)
    plt.show()

# 1a


def localmaxima(data):
    local_maxima_indices = []
    for i in range(len(data)):
        if i != 0 and i != len(data)-1:
            if data[i-1] < data[i] and data[i] > data[i+1]:
                local_maxima_indices.append(i)
    return local_maxima_indices


def plot_local_maxima(t, fn_vals, lm_t, lm_vals, xlim=(0.00, 0.1), ylabel="Amplitude"):
    figure(figsize=(8, 6), dpi=80)
    plt.plot(t, fn_vals)
    plt.scatter(lm_t, lm_vals, c='#ff0000')
    plt.ylabel(ylabel)
    plt.xlabel("Time (s)")
    plt.xlim((0, 0.05))
    plt.show()

# 1b


def crossings(data, threshold, dir: str = "both"):
    crossings_indices = []
    for i in range(len(data)):
        if i != 0:
            if data[i] >= threshold and data[i-1] < threshold:
                if dir == "both" or dir == "negpos":
                    crossings_indices.append(i)
            elif data[i] <= threshold and data[i-1] > threshold:
                if dir == "both" or dir == "posneg":
                    crossings_indices.append(i)
    return crossings_indices


def plot_crossings(t, fn_vals, crossing_t, crossing_vals, threshold=None):
    figure(figsize=(8, 6), dpi=80)
    plt.plot(t, fn_vals)
    if threshold != None:
        plt.hlines(threshold, xmin=0.0, xmax=0.5, colors="#ff8080")
    plt.scatter(crossing_t, crossing_vals, c="#ff0000")
    plt.ylabel("Amplitude")
    plt.xlabel("Time (s)")
    plt.show()

# 2c


def envelope(y, nblocks):
    num_samples_per_block = math.ceil(len(y)/nblocks)
    block_indices = [0]
    y_lower = []
    y_upper = []
    max_per_block = y[0]
    min_per_block = y[0]
    for i in range(len(y)):
        if i != 0 and i % num_samples_per_block == 0:
            block_indices.append(i)
            y_upper.append(max_per_block)
            y_lower.append(min_per_block)
            max_per_block = y[i]
            min_per_block = y[i]
        elif i == len(y) - 1:
            y_upper.append(max_per_block)
            y_lower.append(min_per_block)
        else:
            if y[i] > max_per_block:
                max_per_block = y[i]
            elif y[i] < min_per_block:
                min_per_block = y[i]
    return y_lower, y_upper, block_indices


def plot_raw_audio(wav_data):
    figure(figsize=(8, 6), dpi=80)
    plt.plot(range(len(wav_data)), wav_data, c="000000")
    plt.ylabel("Amplitude")
    plt.xlabel("Frame Index")
    plt.show()


def plot_envelope(y_upper, y_lower, block_indices):
    figure(figsize=(8, 6), dpi=80)
    plt.plot(block_indices, y_upper, c="000000")
    plt.plot(block_indices, y_lower, c="000000")
    plt.ylabel("Amplitude")
    plt.xlabel("Block Index")
    plt.show()
