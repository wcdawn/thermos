import numpy as np
import subprocess
import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams["lines.linewidth"] = 2
# matplotlib.rcParams["mathtext.fontset"] = "stix"
# matplotlib.rcParams["font.family"] = "STIXGeneral"
matplotlib.rcParams["font.size"] = 14


def set_input(txt, refine):
    out = []
    for line in txt:
        if "refine" in line:
            out.append("refine {:d}\n".format(refine))
        else:
            out.append(line)
    return out


def run(exe, inp):
    result = subprocess.run([exe, inp], stdout=subprocess.PIPE)
    return result.stdout.decode("utf-8")


def get_linf(lines):
    for line in lines:
        if "Linf" in line:
            line = line.split()
            return float(line[3])


def get_nx(lines):
    ready = False
    for line in lines:
        if "after refinement" in line:
            ready = True
        elif "Number of cells" in line:
            nx = int(line.split()[3])
            if ready:
                return nx
    # in case refinement not performed
    return nx


if __name__ == "__main__":

    executable = "../../src/thermos.x"
    fname_base = "sph_lin.inp"
    max_refine = 20

    fname_run = fname_base.replace(".inp", "_run.inp")

    runtxt = open(fname_base, "r").readlines()

    linferr = np.zeros(max_refine)
    nx = np.zeros(max_refine, dtype=int)

    for r in range(max_refine):
        inp = set_input(runtxt, r)
        open(fname_run, "w").writelines(inp)

        out = run(executable, fname_run)
        if "CONVERGENCE!" not in out:
            print("failed to converge r=", r)
            sys.exit(1)
        out = out.split("\n")

        linferr[r] = get_linf(out)
        nx[r] = get_nx(out)

        print("r=", r, "linf= {:.2e}".format(linferr[r]))

    with open("result.csv", "w") as f:
        f.write("refine , linferr [K]\n")
        for ridx in range(max_refine):
            f.write("r{:d} , {:.2e}\n".format(ridx, linferr[ridx]))

    plt.figure()
    plt.loglog(nx, np.abs(linferr), "-o")
    plt.loglog(nx, np.abs(linferr[0]) / nx[0] / nx**2, "-k", lw=1, label="_hide")
    plt.xlabel("NX")
    plt.ylabel("$L_\\infty$ error [K]")
    plt.title("Spatial Refinement")
    plt.tight_layout()

    plt.show()
