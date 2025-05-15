#!/usr/bin/env python3
# Roger Volden

'''
Usage:
    python3 plot_knees.py \
            -t bcstats.tsv \
            -o plots/out[.knee.png]
'''

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--tsv', '-t', type=str, required=True,
        help='Output tsv file from bcstats (can be gzipped)'
    )
    parser.add_argument(
        '--output', '-o', type=str, required=True,
        help='Output png prefix (`.knee.png` gets added)'
    )
    parser.add_argument(
        '--max_cells', '-m', default=-1, type=int,
        help='Force an x axis maximum instead of the mpl default'
    )
    parser.add_argument(
        '--estimate_percentile', '-e', default=None, type=int,
        help='Calculates 99th through Nth (inclusive) percentiles for real cell cutoff [None]'
    )
    return parser.parse_args()

def read_tsv(tsv, max_cells):
    zipped, first = False, True
    if tsv.endswith('.gz'):
        import gzip
        fh = gzip.open(tsv, 'rb')
        zipped = True
    else:
        fh = open(tsv, 'r')

    counts, last_real = [], -1
    for line in fh:
        if first:
            first = False
            continue
        if zipped:
            line = line.decode()
        _, nreads, rank, numis, real, _ = line.rstrip().split('\t')
        nreads, rank, numis = int(nreads), int(rank), int(numis)
        if rank >= max_cells and max_cells > 0:
            break
        if real == 'cell':
            last_real = rank
        counts.append(numis)

    fh.close()
    counts.sort(reverse=True)
    return counts, last_real

def plot_from_tsv(counts, ncells, args):
    max_cells, output = args.max_cells, args.output

    plt.figure(figsize=(5, 3.5))
    c = plt.axes([0.125, 0.125, 0.8, 0.8])

    pink = (223/255, 25/255, 149/255)

    c.plot(range(len(counts)), counts, lw=1, color='grey', zorder=10)
    c.plot(range(ncells), counts[:ncells], lw=1.5, color=pink, zorder=11)

    if args.estimate_percentile is not None:
        lb = args.estimate_percentile
        if lb < 80:
            print("WARNING: lower bound for the percentile has a min of 80. Truncating to 80.", file=sys.stderr)
            lb = 80
        colors = [
            (249/255, 157/255, 65/255), (255/255, 102/255, 204/255), (141/255, 109/255, 176/255), (68/255, 168/255, 223/255),
            (0/255, 184/255, 196/255), (106/255, 191/255, 105/255), (155/255, 167/255, 186/255), (225/255, 106/255, 44/255),
            (223/255, 25/255, 149/255), (95/255, 36/255, 159/255), (19/255, 131/255, 198/255), (0/255, 156/255, 162/255),
            (0/255, 157/255, 78/255), (103/255, 111/255, 127/255), (195/255, 74/255, 33/255), (158/255, 13/255, 59/255),
            (67/255, 31/255, 103/255), (0/255, 84/255, 150/255), (13/255, 77/255, 101/255), (0/255, 91/255, 66/255)
        ]
        percs, cutoffs = [], [0]*(99-lb+1)
        for percentile in range(99, lb-1, -1):
            percs.append(np.percentile(counts, percentile)*10)
        print(percs)
        for count in counts:
            for j, p in enumerate(percs):
                if count > p:
                    cutoffs[j] += 1
        for i, cutoff in enumerate(cutoffs):
            c.plot(
                cutoff, counts[cutoff], marker='o', ms=4, c=colors[i],
                label=str(list(range(99, lb-1, -1))[i]), zorder=15
            )
        c.legend(loc=3, prop={'size': 6})

    c.set_xlabel(r'Cell # (log$_{10}$)')
    c.set_xscale('log')
    c.set_yscale('log')
    c.set_ylabel(r'log$_{10}$(# of UMIs)')
    c.set_title('UMIs per cell')

    output += '.knee.png'
    plt.savefig(output, dpi=600)

def main(args):
    counts, ncells = read_tsv(args.tsv, args.max_cells)
    plot_from_tsv(counts, ncells, args)

if __name__ == '__main__':
    args = parse_args()
    main(args)
