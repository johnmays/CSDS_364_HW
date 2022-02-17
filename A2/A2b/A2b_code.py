from collections import OrderedDict
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from scipy import stats
mpl.rcdefaults()
purples = ["#0a0612", "#291749", "#482980", "#673ab7", "#7a52aa", "#9779bd", "#b59fd0", "#d3c5e3"]

def randtimes(N, t1, t2):
    t = []
    for i in range(N):
        t.append((t2-t1) * np.random.uniform() + t1)
    return t


def plotflash(t, t1, t2=None, title="Poisson Distribution of Photons"):
    if t2 == None:
        t2 = np.max(t)
    plt.figure(figsize=(8, 5), dpi=80)
    heads = np.ones(shape=len(t))
    plt.stem(t, heads, markerfmt=" ")
    plt.title(title)
    plt.ylabel("Amplitude")
    plt.xlabel("Time")
    plt.xlim((t1, t2))
    plt.ylim((0, 2))
    plt.yticks([0, 1, 2])
    plt.show()


def plotpdfexp(lam=10):
    t = np.linspace(0.0, 5.0, 501)
    p = lam * np.exp(-lam*t)
    p = p / np.linalg.norm(p, ord=2)
    plt.figure(figsize=(8, 6), dpi=80)
    plt.plot(t, p)
    plt.title("P.D.F. of Exponential Distribution")
    plt.ylabel("$p(\Delta t | \lambda)$")
    plt.xlabel("$\Delta t$")
    plt.xlim((0, 2))
    plt.show()


def randintervals(N, lam, t1):
    t = np.zeros(N)
    for i in range(N):
        if i == 0:
            t[i] = t1 + np.random.exponential(scale=1/lam)
        else:
            t[i] = t[i-1] + np.random.exponential(scale=1/lam)
    return t


def broadcastfactorial(lis):
    # accepts a single value, or a list, or a 1D numpy array
    if type(lis) == list or type(lis) == np.ndarray:
        for i in range(len(lis)):
            lis[i] = math.factorial(lis[i])
        return np.array(lis)
    else:
        lis = math.factorial(lis)
        return lis


def pdfphotons(K=1, lam=10, T=.100):
    if type(K) == list or type(K) == np.ndarray:
        K_copy = np.copy(K)
        return (((lam*T)**K_copy)/(broadcastfactorial(K_copy)))*np.exp(-lam*T)
    else:
        return (((lam*T)**K)/(broadcastfactorial(K)))*np.exp(-lam*T)


def plotbarpdfphotons(K=np.arange(0, 10), lam=10, T=.100, xlimit=None, title="Probability of observing $n$ events"):
    p = pdfphotons(K=K, lam=lam, T=T)
    plt.figure(figsize=(8, 6), dpi=80)
    plt.bar(K, p)
    plt.title(title)
    plt.ylabel("$p(n | \lambda, T)$")
    plt.xlabel("$n$")
    if xlimit == None:
        plt.xlim((-1, len(K)))
    else:
        plt.xlim(xlimit)
    plt.show()


def detectionprob(K, lam=40, T=0.1):
    return 1-stats.poisson.cdf(k=(K-1), mu=lam*T)


def detectionprobsum(K, lam=40, T=0.1):
    # this function performs the series calculation of the cdf, more or less to confirm if I am doing it correctly or not
    K_index = K
    pdf = pdfphotons(K_index, lam=lam, T=T)
    sum = pdf
    while pdf > 0.0001:
        K_index += 1
        pdf = pdfphotons(K=K_index, lam=lam, T=T)
        sum += pdf
    return sum


def plotbarcdfphotons(K: list = np.arange(0, 10), lam=40, T=0.1, xlimit=None, title="Probability of observing $n\geq K$ events"):
    # the argument K should be a list or 1D numpy array
    p = []
    for k in K:
        p.append(detectionprob(K=k, lam=lam, T=T))
    plt.figure(figsize=(8, 6), dpi=80)
    plt.bar(K, p)
    plt.title(title)
    plt.ylabel("$p(n \geq K| \lambda, T)$")
    plt.xlabel("$K$")
    if xlimit == None:
        plt.xlim((-1, len(K)))
    else:
        plt.xlim(xlimit)
    plt.show()


def lightflash(lam=100, t1=0.8, t2=2.2):
    # submit lam parameter in photons/millisecond;  t1 and t2 parameters in milliseconds
    t = t1 + np.random.exponential(scale=1/lam)
    photon_times = []
    while t < t2:
        photon_times.append(t)
        t += np.random.exponential(scale=1/lam)
    return photon_times


def get_detected_photon_times(photon_times, alpha):
    # given an array of photon times and an alpha to describe probability of detection, this should return a subset of those times reflecting alpha
    detected_photon_times = []
    for pt in photon_times:
        if np.random.uniform() <= alpha:
            detected_photon_times.append(pt)
    return detected_photon_times


def cut_list(lis, min, max):
    new_lis = []
    for i in range(len(lis)):
        if lis[i] >= min and lis[i] <= max:
            new_lis.append(lis[i])
    return new_lis


def plotHSPsimulation(lam=100, t1=0.8, t2=2.2, f1=0.8, f2=2.2, s1=1.0, s2=2.0, alpha=0.06):
    photon_times = lightflash(lam=lam, t1=t1, t2=t2)
    photon_times_abbreviated = cut_list(photon_times, f1, f2)
    plt.figure(figsize=(8, 4), dpi=80)
    heads = np.ones(shape=len(photon_times_abbreviated))
    plt.stem(photon_times_abbreviated, heads, linefmt="#ccccff",
             markerfmt=" ", basefmt="#ccccff")
    plt.title(
        "Photons in Interval $[f_1, f_2]=[{f_1},{f_2}]$".format(f_1=f1, f_2=f2))
    plt.ylabel("Amplitude")
    plt.xlabel("Time (milliseconds)")
    plt.xlim((t1, t2))
    plt.ylim((0, 2))
    plt.show()

    photon_times_shutter = cut_list(photon_times_abbreviated, s1, s2)
    plt.figure(figsize=(8, 4), dpi=80)
    heads = np.ones(shape=len(photon_times_shutter))
    plt.stem(photon_times_shutter, heads, linefmt="#9999ff",
             markerfmt=" ", basefmt="#9999ff")
    plt.title("Photons that Pass through Shutter open during $[s_1, s_2]=[{s_1},{s_2}]$".format(
        s_1=s1, s_2=s2))
    plt.ylabel("Amplitude")
    plt.xlabel("Time (milliseconds)")
    # plt.xlim((s1, s2))
    plt.xlim((t1, t2))
    plt.ylim((0, 2))
    plt.show()

    photon_times_absorbed = get_detected_photon_times(
        photon_times_shutter, alpha=alpha)
    plt.figure(figsize=(8, 4), dpi=80)
    heads = np.ones(shape=len(photon_times_absorbed))
    plt.stem(photon_times_absorbed, heads, linefmt="#5555ff",
             markerfmt=" ", basefmt="#5555ff")
    plt.title("Photons that are Absorbed While Shudder is Open")
    plt.ylabel("Amplitude")
    plt.xlabel("Time (milliseconds)")
    # plt.xlim((s1, s2))
    plt.xlim((t1, t2))
    plt.ylim((0, 2))
    plt.show()


def probseeing(I, alpha=0.06, K=6):
    return 1-stats.poisson.cdf(k=(K-1), mu=alpha*I)


def plotdetectioncurve(alpha=0.5, K=6, seperatecurves = True):
    if type(alpha) != list and type(alpha) != np.ndarray:
        alpha = [alpha]
        K = [K]
    # if type(alpha) == list or type(alpha) == np.ndarray:
    # assert type(K) == list or type(K) == np.ndarray
    assert len(K) == len(alpha)
    plt.figure(figsize=(8, 5), dpi=80)
    for i in range(len(alpha)):
        p = []
        for I in np.linspace(0.01, 100, 10000):
            p.append(probseeing(I, alpha=alpha[i], K=K[i]))
        if seperatecurves:
            plt.plot(np.linspace(0.01, 100, 10000), p, label="alpha={a}, K={k}".format(a=alpha[i], k=K[i]), c=purples[i])
        else:
            plt.plot(np.linspace(0.01, 100, 10000), p, c="#ccccff")
    plt.title("Probability of Detection of a Flash w.r.t. Intensity")
    plt.ylabel("$p$(Detection|Flash)")
    plt.xlabel("Intensity")
    plt.xlim(0.01, 100)
    plt.xscale('log')
    if seperatecurves:
        plt.legend()
    plt.show()


def plotfit(alpha=3, K=3):
    plt.figure(figsize=(8, 5), dpi=80)

    # first plotting our expiremental results
    e_alpha = [24.1, 37.6, 58.6, 91.0, 141.9, 221.3]
    e_K = [0.0, 4.0, 18.0, 54.0, 94.0, 100.0]

    for i in range(len(e_alpha)):
        e_p = []
        for I in np.linspace(0.01, 100, 10000):
            e_p.append(probseeing(I, alpha=e_alpha[i], K=e_K[i]))
        plt.plot(np.linspace(0.01, 100, 10000), e_p,
                 c="#ccccff", label="experimental data")

    # then plotting the fitted K values
    if type(alpha) != list and type(alpha) != np.ndarray:
        alpha = [alpha]
        K = [K]
    assert len(K) == len(alpha)
    for i in range(len(alpha)):
        p = []
        for I in np.linspace(0.01, 100, 10000):
            p.append(probseeing(I, alpha=alpha[i], K=K[i]))
        plt.plot(np.linspace(0.01, 100, 10000), p, c="#8888ff", label="fits")

    plt.title("Some Fits: Probability of Detection of a Flash w.r.t. Intensity")
    plt.ylabel("$p$(Detection|Flash)")
    plt.xlabel("Intensity")
    plt.xlim(0.01, 100)
    plt.xscale('log')
    # Subsequent 3 lines mostly taken from: https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.show()
