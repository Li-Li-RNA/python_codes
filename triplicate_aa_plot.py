#!/usr/bin/env python
import sys
import numpy as np
import scipy.optimize
import scipy.stats
import matplotlib.pyplot as plt


def fitfunc(x, p0, p1):
    return p0 * (1 - np.exp(-p1*x))


def fit(conc_array, time_array):
    """
    :param conc_array: array of the product concentration
    :param time_array: array of time
    :return: pars, the fitted rate and r2, the R2 of the fit
    """
    pars, pcov = scipy.optimize.curve_fit(fitfunc, time_array, conc_array)

    # Calculate r squared
    fitted_data = fitfunc(time_array, pars[0], pars[1])
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(fitted_data, conc_array)
    print("max = %.3E" % pars[0])
    print("k = %.3E" % pars[1])
    print("R2 = %.3E" % r_value**2)
    return pars, r_value**2


class Reaction:
    def __init__(self, time_array, cpm_array, sa, aa_conc):
        self.cpm_array = cpm_array  # an array stores the raw CPM
        self.time_array = time_array  # an array stores the time
        self.sa = sa  # specific activity
        self.pars = 0  # fitted parameters
        self.r2 = 0  # R squared
        self.normalized_cpm_array = []  # an array stores the raw CPM - the CPM of time 0
        self.aa_conc = aa_conc  # concentration of amino acid corresponding to the specific activity
        self.conc_array = []  # array stores the calculated concentration of charged tRNA

    def calc_conc(self):
        # Convert the CPM to concentrations of charged tRNA
        cpm0 = self.cpm_array[0]
        self.normalized_cpm_array = self.cpm_array - cpm0
        self.conc_array = self.normalized_cpm_array / self.sa * self.aa_conc

    def fit_plateau(self):
        self.calc_conc()
        # determine the maximum concentration of charged tRNA
        self.pars, self.r2 = fit(self.conc_array, self.time_array)
        return self.pars


def read_data(in_file_name):
    """
    :param in_file_name: input file name
    :return: a list of Reaction objects
    Parse the input file
    The input file has three types of input line:
    1. conc XX
    Following keyword "conc", XX is the concentration of the total amino acid.
    2. YY ZZ
    With no keyword, this specifies the CPM (ZZ) at a given time point (YY).
    3. sa KK
    Following keyword "sa", KK is the CPM of the specific activity pad.
    Entries of one experiment starts with 1st type of line, followed by n 2nd
    type of lines, and ended by a 3rd type of line.
    An input file may contain multiple replicates of the same experiment.  Their
    average and standard deviations will be calculated.
    """
    # Initialization
    rxns = []  # an array contains objects of Reaction class

    with open(in_file_name, 'r') as ifile:
        aa_conc = 0  # concentration of total amino acid
        for line in ifile.readlines():
            entries = line.split()
            time = entries[0]
            cpm = entries[1]
            if time == "conc":
                # when the first item is "conc", initiate a new time and cpm array
                # and set the conc to be the second item
                time_array = []
                cpm_array = []
                aa_conc = float(cpm)
            elif time == "sa":
                # when the first item is "sa", generate a new object of the Reaction class
                # and append it to the rxns array
                rxns.append(Reaction(np.array(time_array), np.array(cpm_array), float(cpm), aa_conc))
            else:
                # otherwise, read the first item as time, and the second item as cpm.
                """ Add an assertion """
                time_array.append(float(time))
                cpm_array.append(float(cpm))
    return rxns


def print_data(counts_array, time_array, out_file_prefix):
    """
    :param counts_array: array of raw counts
    :param time_array: array of time
    :param out_file_prefix: output file prefix
    :return:
    """
    # print the data
    num_rxns, num_tp, num_lanes = counts_array.shape
    o_data_file_name = out_file_prefix + "-frac.dat"
    with open(o_data_file_name, 'w') as ofile:
        for rxn in range(num_rxns):
            ofile.write("Time/Fraction\t1\t2\t3\t4\n")
            for i in range(num_tp):
                ofile.write("%.1f\t" % time_array[i])
                for l in range(num_lanes):
                    ofile.write("\t%.3f" % counts_array[rxn, i, l])
                ofile.write("\n")


def plot_multiple_fit(rxns, fig_prefix):
    time_array = rxns[0].time_array
    num_rxns = len(rxns)
    num_data_points = len(time_array)
    # Concatenate all concentrations into an M x N array, where
    #   M is the number of reactions
    #   N is the number of time points
    concat_conc_array = []
    for rxn in rxns:
        concat_conc_array.append(rxn.conc_array)
    conc_matrix = np.reshape(concat_conc_array, [num_rxns, num_data_points])

    # Calculate the average and standard deviation of each time point
    ave_conc_array = np.average(conc_matrix, axis=0)
    std_conc_array = np.std(conc_matrix, axis=0)

    # plot
    num_points = 1001
    tt = np.linspace(0.0, np.amax(time_array), num_points)
    ave_pars, r2 = fit(ave_conc_array, time_array)
    fitted_conc = fitfunc(tt, ave_pars[0], ave_pars[1])

    # Plot
    fig = plt.figure(figsize=(4, 4))
    left, width = 0.2, 0.7
    bottom, height = 0.2, 0.7
    rect_scatter = [left, bottom, width, height]
    ax = plt.axes(rect_scatter)

    # plot the raw data
    ax.errorbar(time_array, ave_conc_array, yerr=std_conc_array, marker='o', mfc='b', color='b',
                ecolor='b', mec='b', linestyle='none')
    # plot the fitted curve
    plt.plot(tt, fitted_conc, ls='-', color='k', linewidth=1)

    # x,y-axis limit (need to be automated)
    ax.set_xlim(0, np.floor(np.amax(time_array))+1)
    ax.set_ylim(0, 4)

    plt.xticks([10, 20, 30])
    ticklabels = ax.get_xticklabels()
    ticklabels.extend(ax.get_yticklabels())
    for label in ticklabels:
        label.set_fontsize(12)
    plt.xlabel('Time (min)', fontsize=16)
    plt.ylabel('Charged tRNA (uM)', fontsize=16)

    majorticks = ax.xaxis.get_majorticklines()
    for majortickline in majorticks:
        majortickline.set_linewidth(1)
        majortickline.set_markeredgecolor('k')
        majortickline.set_markeredgewidth(1)

    majorticks = ax.yaxis.get_majorticklines()
    for majortickline in majorticks:
        majortickline.set_linewidth(1)
        majortickline.set_markeredgecolor('k')
        majortickline.set_markeredgewidth(1)

    fig.savefig(fig_prefix+'.eps')
    fig.savefig(fig_prefix+'.svg')


def main():
    # Read variables from the command line
    try:
        in_file_name = sys.argv[1]
        out_file_prefix = sys.argv[2]
    except IOError:
        print("Usage: ", sys.argv[0], "input file, output file prefix")
        sys.exit(1)

    # read the data
    rxns = read_data(in_file_name)

    for rxn in rxns:
        # Calculate the concentration of charged tRNA from CPM
        rxn.calc_conc()

    plot_multiple_fit(rxns, out_file_prefix)


if __name__ == "__main__":
    main()
