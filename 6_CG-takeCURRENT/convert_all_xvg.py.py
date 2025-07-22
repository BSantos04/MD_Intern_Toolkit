#!/usr/bin/env python3
"""
convert_all_xvg.py

Scan the three mapping directories in the current folder and
convert every .xvg file into a PNG, using its @-headers
for title/xlabel/ylabel (falling back to Title/xx/yy).
"""

import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt

# --- Helpers to read data & metadata ----------------------------------------

def read_xvg_data(path):
    """Read numerical columns from an .xvg (skip lines starting with @ or #)."""
    arr = []
    with open(path, 'r') as f:
        for L in f:
            if L.startswith(('@', '#')):
                continue
            parts = L.strip().split()
            if parts:
                arr.append([float(x) for x in parts])
    return np.array(arr)

def read_xvg_metadata(path):
    """Extract title, x-label, and y-label from @-lines; else None."""
    meta = {'title': None, 'xlabel': None, 'ylabel': None}
    title_re  = re.compile(r'@.*title\s+"([^"]+)"')
    xl_re     = re.compile(r'@.*xaxis\s+label\s+"([^"]+)"')
    yl_re     = re.compile(r'@.*yaxis\s+label\s+"([^"]+)"')
    with open(path, 'r') as f:
        for L in f:
            if not L.startswith('@'):
                continue
            if meta['title'] is None:
                m = title_re.match(L)
                if m:
                    meta['title'] = m.group(1)
                    continue
            if meta['xlabel'] is None:
                m = xl_re.match(L)
                if m:
                    meta['xlabel'] = m.group(1)
                    continue
            if meta['ylabel'] is None:
                m = yl_re.match(L)
                if m:
                    meta['ylabel'] = m.group(1)
                    continue
            if all(meta.values()):
                break
    return meta

# --- Plotting function ------------------------------------------------------

def plot_xvg_to_png(xvg_path, png_path):
    data = read_xvg_data(xvg_path)
    if data.size == 0:
        raise RuntimeError(f"No data in {xvg_path!r}")
    meta = read_xvg_metadata(xvg_path)

    title  = meta['title']   or "Title"
    xlabel = meta['xlabel']  or "xx"
    ylabel = meta['ylabel']  or "yy"

    x = data[:,0]
    ys = data[:,1:]

    plt.figure()
    for i in range(ys.shape[1]):
        plt.plot(x, ys[:,i], label=f"Col {i+2}")
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    # Only add legend if there's more than one series
    if ys.shape[1] > 1:
        plt.legend()
    plt.tight_layout()
    plt.savefig(png_path, dpi=300)
    plt.close()

# --- Main loop ---------------------------------------------------------------

if __name__ == "__main__":
    # these three directories must exist under cwd
    input_dirs = ["bonds_mapped", "angles_mapped", "dihedrals_mapped"]

    for d in input_dirs:
        if not os.path.isdir(d):
            print(f"Warning: directory {d!r} not found, skipping.", file=sys.stderr)
            continue

        for fname in sorted(os.listdir(d)):
            if not fname.endswith(".xvg"):
                continue
            xvg_path = os.path.join(d, fname)
            png_name = os.path.splitext(fname)[0] + ".png"
            png_path = os.path.join(d, png_name)

            try:
                plot_xvg_to_png(xvg_path, png_path)
                print(f"✓ {xvg_path} → {png_path}")
            except Exception as e:
                print(f"✗ failed {xvg_path}: {e}", file=sys.stderr)