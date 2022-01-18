import math
import matplotlib

def g (x, mu = 0, sigma = 1.0):
    return (1 / math.sqrt(2*math.pi*(sigma**2)))*math.exp(-((x-mu)**2)/(2*(sigma**2)))

