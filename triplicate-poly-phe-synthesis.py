#!/usr/bin/env python
import sys 
import numpy as np
import scipy.optimize
import scipy.stats
import matplotlib.pyplot as plt


def plot_bar(data, fig_prefix):
    """
    Plot the data as a grouped bar graph
    :param data: input data
    :param fig_prefix: output figure prefix
    :return:
    """
    n_row, n_col = np.shape(data)
    ind = np.arange(n_col/2)
    w = 0.35

    # Calculate the mean and stdev of the data
    mean_data = np.average(data, axis=0)
    stdev_data = np.std(data, axis=0)

    # plot
    fig = plt.figure(figsize=(4, 4))
    left, width = 0.2, 0.7
    bottom, height = 0.2, 0.7
    rect_scatter = [left, bottom, width, height]
    ax = plt.axes(rect_scatter)

    rects1 = ax.bar(ind, mean_data[0:6], w, color="lightgrey", yerr=stdev_data[0:6], ecolor='k')
    rects2 = ax.bar(ind+w, mean_data[6:12], w, color="darkgrey", yerr=stdev_data[6:12], ecolor='k')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Polymerized Phe (pmol)')
    ax.set_ylim(-0.1, 4)
    ax.set_xticks(ind + w / 2)
    ax.set_xticklabels(('R1', 'R2', 'R3', 'R4', 'R5', 'R6'))

    ax.legend((rects1[0], rects2[0]), ('5 min', '20 min'), loc='upper left')

    fig.savefig(fig_prefix + '.eps')
    fig.savefig(fig_prefix + '.svg')


def main():
    # Read variables from the command line
    try:
        in_file_name = sys.argv[1]
        out_file_prefix = sys.argv[2]
    except IOError:
        print("Usage: ", sys.argv[0], "input file, output file prefix")
        sys.exit(1)

    # read the data
    # The input file has the following format:
    # Run 1 5 min Run 1 20 min
    # Run 2 5 min Run 2 20 min
    # ...
    data = np.loadtxt(in_file_name)

    plot_bar(data, out_file_prefix)


if __name__ == "__main__":
    main()
