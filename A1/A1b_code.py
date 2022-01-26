from turtle import forward
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def plot_envelope(y_upper, y_lower, block_indices):
    figure(figsize=(8, 6), dpi=80, forward=True)
    plt.plot(block_indices, y_upper, c="000000")
    plt.plot(block_indices, y_lower, c="000000")
    plt.ylabel("Amplitude")
    plt.xlabel("Block Index")
    plt.show()
