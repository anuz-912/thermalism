import sys
import os

import numpy as np
import matplotlib
matplotlib.use("Agg") 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def load_csv(path):
    # columns: x, y, T_transient, T_steady, T@t0, T@t1, ... (one row per grid node)
    with open(path) as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)

    col = int(data[:, 0].max()) + 1
    row = int(data[:, 1].max()) + 1

    T = data[:, 2].reshape(row, col)   # transient final state, T[i][j] with i = y, j = x
    S = data[:, 3].reshape(row, col)   # steady state

    snaps = []                         # (time, matrix) for every saved snapshot
    for k, name in enumerate(header):
        if name.startswith("T@"):
            snaps.append((float(name[2:]), data[:, k].reshape(row, col)))
    return T, S, snaps


def main():
    if len(sys.argv) < 2:
        print("usage: python visualizer.py <file.csv>")
        sys.exit(1)

    csv_file = sys.argv[1]
    image_file = os.path.splitext(csv_file)[0] + ".png"   # same name as the csv file

    T, S, snaps = load_csv(csv_file)
    row, col = T.shape
    diff = T - S

    n = len(snaps)
    ncols = 4
    snap_rows = int(np.ceil(n / ncols))
    rows = snap_rows + 1

    fig = plt.figure(figsize=(18, 4.6 * rows))
    gs = GridSpec(rows, ncols, figure=fig)
    fig.suptitle("ThermalSim: 2D heat conduction (%d x %d grid), plate starts at %.6g C" %
                 (row, col, snaps[0][1][row // 2, col // 2]), fontsize=15)

    # 1. evolution of the plate: initial condition -> intermediate times -> steady state
    for k, (t, M) in enumerate(snaps):
        a = fig.add_subplot(gs[k // ncols, k % ncols])
        im = a.imshow(M, origin="lower", cmap="inferno")
        if k == 0:
            title = "Initial condition (t = 0 s)"
        elif k == n - 1:
            title = "Steady state reached (t = %.4g s)" % t
        else:
            title = "t = %.4g s" % t
        a.set_title(title)
        a.set_xlabel("x (grid index)")
        a.set_ylabel("y (grid index)")
        fig.colorbar(im, ax=a, label="Temperature")

    # 2. steady-state solution (SOR)
    a = fig.add_subplot(gs[rows - 1, 0])
    im = a.imshow(S, origin="lower", cmap="inferno")
    a.set_title("Steady-state solver (SOR)")
    a.set_xlabel("x (grid index)")
    a.set_ylabel("y (grid index)")
    fig.colorbar(im, ax=a, label="Temperature")

    # 3. difference between the final transient state and the steady-state solver
    a = fig.add_subplot(gs[rows - 1, 1])
    lim = max(abs(diff).max(), 1e-12)
    im = a.imshow(diff, origin="lower", cmap="coolwarm", vmin=-lim, vmax=lim)
    a.set_title("Transient - Steady (max %.1e)" % abs(diff).max())
    a.set_xlabel("x (grid index)")
    a.set_ylabel("y (grid index)")
    fig.colorbar(im, ax=a, label="Temperature difference")

    # 4. vertical centre line at every saved time
    a = fig.add_subplot(gs[rows - 1, 2:])
    colours = plt.cm.viridis(np.linspace(0, 0.9, n))
    for k, (t, M) in enumerate(snaps):
        a.plot(np.arange(row), M[:, col // 2], color=colours[k], label="t = %.4g s" % t)
    a.plot(np.arange(row), S[:, col // 2], "k--", label="Steady state")
    if max(abs(M).max() for _, M in snaps) > 10 * max(abs(S).max(), 1.0):
        a.set_yscale("symlog", linthresh=max(abs(S).max(), 1.0))   # initial temperature dwarfs the rest
    a.set_title("Temperature along the vertical centre line over time")
    a.set_xlabel("y (grid index)")
    a.set_ylabel("Temperature")
    a.grid(True, alpha=0.3)
    a.legend(fontsize=8, ncol=2)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(image_file, dpi=130)
    print("Plot saved to " + image_file)


if __name__ == "__main__":
    main()
