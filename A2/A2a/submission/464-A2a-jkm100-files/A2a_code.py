import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


def genwaveform(N=100, alpha=0.4, A=1, sigma=1.0, noisetype="Gaussian"):
    signal = []
    event_locations = []
    for i in range(N):
        if noisetype == "Gaussian":
            noise_val = np.random.normal(scale=sigma)
        elif noisetype == "uniform":
            noise_val = np.random.uniform(low=-sigma/2, high=sigma/2)
        if (np.random.uniform(low=0.0, high=1.0) < alpha):
            signal_val = A
            event_locations.append(i)
        else:
            signal_val = 0
        signal.append(signal_val + noise_val)
    return signal, event_locations


def plotwaveform(signal, event_locations=[], title="Signal with Noise"):
    plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(range(len(signal)), signal, label="waveform")
    plt.vlines(event_locations, ymin=np.min(signal),
               ymax=np.max(signal), colors="#a3a3ff", label="event")
    plt.title(title)
    plt.ylabel("Amplitude")
    plt.xlabel("Index")
    plt.legend()
    plt.show()


def dectectioncounts(si=[], y=[], thresh=1.0):
    tp, fn, fp, tn = 0, 0, 0, 0
    for i in range(len(y)):
        if i in si:
            signal = True
        else:
            signal = False
        if y[i] >= thresh:
            detection = True
        else:
            detection = False

        if signal and detection:
            tp += 1
        elif signal and not detection:
            fn += 1
        elif not signal and detection:
            fp += 1
        else:
            tn += 1
    return tp, fn, fp, tn


def dectectionindices(si=[], y=[], thresh=1.0):
    tp, fn, fp, tn = [], [], [], []
    for i in range(len(y)):
        if i in si:
            signal = True
        else:
            signal = False
        if y[i] >= thresh:
            detection = True
        else:
            detection = False

        if signal and detection:
            tp.append(i)
        elif signal and not detection:
            fn.append(i)
        elif not signal and detection:
            fp.append(i)
        else:
            tn.append(i)
    return tp, fn, fp, tn


def plotdetectioncounts(signal=[], thresh=1.0, tp=[], fn=[], fp=[], tn=[], title="Threshold Detection Successes"):
    plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(range(len(signal)), signal, label="waveform")
    plt.hlines(thresh, xmin=0, xmax=len(signal), colors="#a3a3ff")
    plt.vlines(tp, ymin=np.min(signal), ymax=np.max(
        signal), colors="lime", label="True Positive")
    # plt.vlines(tn, ymin=np.min(signal), ymax=np.max(signal), colors="lime", label="True Negative", linestyles={'dashed'})
    plt.vlines(fp, ymin=np.min(signal), ymax=np.max(
        signal), colors="red", label="False Positive")
    plt.vlines(fn, ymin=np.min(signal), ymax=np.max(signal),
               colors="red", label="False Negative", linestyles={'dashed'})
    plt.title(title)
    plt.ylabel("Amplitude")
    plt.xlabel("Index")
    plt.legend()
    plt.show()


def falsepos(thresh=1.0, sigma=1.0):
    return 1 - 0.5*(1+math.erf(thresh/(sigma*math.sqrt(2))))


def falseneg(thresh=1.0, A=1.0, sigma=1.0):
    return 0.5*(1+math.erf((thresh-A)/(sigma*math.sqrt(2))))


def falseposestimate(N=100, thresh=1.0, sigma=1.0):
    false_pos_rates = []
    for i in range(N):
        signal, event_locations = genwaveform(N=100, sigma=sigma)
        tp, fn, fp, tn = dectectioncounts(
            si=event_locations, y=signal, thresh=thresh)
        false_pos_rates.append(fp/(fp+tn))
    return np.average(false_pos_rates)


def falsenegestimate(N=100, thresh=1.0, A=1.0, sigma=1.0):
    false_neg_rates = []
    for i in range(N):
        signal, event_locations = genwaveform(N=100, A=A, sigma=sigma)
        tp, fn, fp, tn = dectectioncounts(
            si=event_locations, y=signal, thresh=thresh)
        false_neg_rates.append(fn/(fn+tp))
    return np.average(false_neg_rates)


# N=1000, alpha=0.4, A=1, sigma=1.0, noisetype="Gaussian"
def plotROC(thresh_min: int = -2, thresh_max: int = 5, resolution: int = 10, sigma=1.0, A=1.0, title="ROC Curve"):
    thresholds = np.linspace(
        start=thresh_min, stop=thresh_max, num=resolution*(thresh_max-thresh_min)+1)
    false_positives = []
    true_positives = []
    for thresh in thresholds:
        false_positives.append(falsepos(thresh, sigma))
        true_positives.append(1-falseneg(thresh, A, sigma))
    plt.figure(figsize=(6, 6), dpi=80)
    plt.plot(false_positives, true_positives)
    plt.title(title)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.show()
