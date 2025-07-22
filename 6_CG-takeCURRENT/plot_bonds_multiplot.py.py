#!/usr/bin/env python3
"""
plot_bonds_multiplot.py

Reads the same .xvg files used in your gnuplot tutorial and
lays them out in a 2×4 grid, saving the result as a PDF.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Use Helvetica 12pt to match your gnuplot terminal settings
plt.rc('font', family='Helvetica')
plt.rcParams['font.size'] = 12

def read_xvg(fn):
    """
    Read an .xvg, skipping any line starting with '#' or '@',
    and return an (N×M) NumPy array.
    """
    rows = []
    with open(fn) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if parts:
                rows.append([float(x) for x in parts])
    return np.array(rows)

def main():
    # Configuration taken from your gnuplot script
    base_dirs = [
        "5_target-distr/bonds_mapped",
        "6_CG-takeCURRENT/bonds_mapped",
    ]
    colors = ["#1E90FF", "#DC143C"]
    labels = ["COG-mapped", "Martini"]
    titles = [
        r"$b_1$, 1-2",
        r"$b_2$, 1-4",
        r"$b_3$, 2-3",
        r"$b_4$, 2-5",
        r"$b_5$, 2-6",
        r"$b_6$, 3-6",
        r"$b_7$, 5-6",
    ]
    xlims = [
        (0.2, 0.4),
        (0.3, 0.5),
        (0.2, 0.4),
        (0.4, 0.6),
        (0.4, 0.6),
        (0.4, 0.6),
        (0.2, 0.4),
    ]

    # Create figure with 2 rows × 4 cols, size 8.25×3.3 inches
    fig, axes = plt.subplots(2, 4, figsize=(8.25, 3.3), sharey=True)
    axes = axes.flatten()

    for i in range(7):
        ax = axes[i]
        for base_dir, color in zip(base_dirs, colors):
            xvg = os.path.join(base_dir, f"distr_bond_{i}.xvg")
            data = read_xvg(xvg)
            if data.size == 0:
                continue
            ax.plot(data[:, 0], data[:, 1],
                    color=color, lw=2)

        ax.set_title(titles[i])
        ax.set_xlim(*xlims[i])
        ax.set_ylim(0, 200)
        # xtics every 0.05 as in gnuplot
        ax.set_xticks(np.arange(xlims[i][0],
                                xlims[i][1] + 1e-8,
                                0.05))
        ax.grid(False)       # unset grid
        # no legend (gnuplot had `unset key`)

    # Hide the unused 8th subplot
    axes[7].axis('off')

    fig.tight_layout()
    fig.savefig("AAvsCG-bond-tutorial-4x2.pdf", format='pdf', bbox_inches='tight')
    plt.close(fig)
    print("Written: AAvsCG-bond-tutorial-4x2.pdf")

if __name__ == "__main__":
    main()