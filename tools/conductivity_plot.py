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
    k = dat[:, 1]

    plt.figure()
    plt.plot(x, k)
    plt.xlabel("$x$ [cm]")
    plt.ylabel("$k(T(x))$ [W/cm/K]")
    plt.title("Thermos Thermal Conductivity")
    plt.tight_layout()
    plt.savefig(fname.replace(".csv", "." + extension), dpi=resolution)

    plt.show()
