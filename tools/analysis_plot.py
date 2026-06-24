import numpy as np
import matplotlib.pyplot as plt
import sys

import plot_settings

if __name__ == "__main__":

    extension = "png"
    resolution = 600

    fname = sys.argv[1]

    dat = np.loadtxt(fname, delimiter=",", skiprows=1)
    x = dat[:, 0]
    T = dat[:, 1]
    Texact = dat[:, 2]

    plt.figure()
    plt.plot(x, T, label="Thermos")
    plt.plot(x, Texact, label="Exact")
    plt.legend()
    plt.xlabel("$x$ [cm]")
    plt.ylabel("$T(x)$ [K]")
    plt.title("Thermos Temperature")
    plt.tight_layout()
    plt.savefig(fname.replace(".csv", "." + extension), dpi=resolution)

    plt.figure()
    plt.plot(x, T - Texact)
    plt.xlabel("$x$ [cm]")
    plt.ylabel("$T(x) - T_{exact}(x)$ [K]")
    plt.title("Thermos Temperature Difference")
    plt.tight_layout()
    plt.savefig(fname.replace(".csv", "_diff." + extension), dpi=resolution)

    plt.show()
