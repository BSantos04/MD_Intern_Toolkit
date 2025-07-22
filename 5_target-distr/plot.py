#!/usr/bin/env python3
"""
plot_xvg_to_png.py

Reads a GROMACS .xvg file (skipping lines beginning with '#'),
extracts the plot title and axis labels from the @-directives, 
plots each data column against the first column, and writes the result
to a PNG file.
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

def read_xvg(filename):
    """
    Read an .xvg file into a NumPy array, skipping data-header lines.
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('@', '#')):
                continue
            parts = line.strip().split()
            if parts:
                data.append([float(x) for x in parts])
    return np.array(data)

def read_xvg_metadata(filename):
    """
    Parse @-lines for title, xaxis label, and yaxis label.
    Returns a dict with keys 'title', 'xlabel', 'ylabel'.
    """
    meta = {'title': None, 'xlabel': None, 'ylabel': None}
    title_re   = re.compile(r'@.*title\s+"([^"]+)"')
    xlabel_re  = re.compile(r'@.*xaxis\s+label\s+"([^"]+)"')
    ylabel_re  = re.compile(r'@.*yaxis\s+label\s+"([^"]+)"')
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('@'):
                continue
            if meta['title'] is None:
                m = title_re.match(line)
                if m:
                    meta['title'] = m.group(1)
                    continue
            if meta['xlabel'] is None:
                m = xlabel_re.match(line)
                if m:
                    meta['xlabel'] = m.group(1)
                    continue
            if meta['ylabel'] is None:
                m = ylabel_re.match(line)
                if m:
                    meta['ylabel'] = m.group(1)
                    continue
            # If we've got all three, we can stop looking
            if all(meta.values()):
                break
    return meta

def plot_xvg_to_png(xvg_path, png_path):
    # load data
    data = read_xvg(xvg_path)
    if data.size == 0:
        raise ValueError(f"No data found in {xvg_path}")

    # load metadata
    meta = read_xvg_metadata(xvg_path)
    title   = meta['title']   or "Title"
    xlabel  = meta['xlabel']  or "xx"
    ylabel  = meta['ylabel']  or "yy"

    # first column = x, rest = y-series
    x  = data[:, 0]
    ys = data[:, 1:]

    # plot
    plt.figure()
    for col in range(ys.shape[1]):
        plt.plot(x, ys[:, col], label=f'Col {col+2}')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()

    # save
    plt.savefig(png_path, dpi=300)
    plt.close()

if __name__ == '__main__':
    # Process all .xvg files in these directories
    input_dirs = ['bonds_mapped', 'angles_mapped', 'dihedrals_mapped']
    for input_dir in input_dirs:
        for root, dirs, files in os.walk(input_dir):
            for fname in files:
                if fname.endswith('.xvg'):
                    xvg_path = os.path.join(root, fname)
                    base, _ = os.path.splitext(fname)
                    png_path = os.path.join(root, base + '.png')
                    try:
                        plot_xvg_to_png(xvg_path, png_path)
                        print(f"Converted {xvg_path} -> {png_path}")
                    except Exception as e:
                        print(f"Error converting {xvg_path}: {e}")