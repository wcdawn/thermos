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
    q = dat[:, 1]

    plt.figure()
    plt.plot(x, T)
    plt.xlabel("$x$ [cm]")
    plt.ylabel("$q(x)$ [W/cm^3]")
    plt.title("Thermos Source")
    plt.tight_layout()
    plt.savefig(fname.replace(".csv", "." + extension), dpi=resolution)

    plt.show()
