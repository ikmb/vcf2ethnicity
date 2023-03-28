#!/usr/bin/env python

import colorsys
import getopt
import sys

import matplotlib.pyplot as plot
import numpy as np
import pandas

from Simulation.msprime import MSPrimeClass

"""
please cite 
Anil Raj, Matthew Stephens, and Jonathan K. Pritchard. fastSTRUCTURE: Variational Inference of Population Structure in 
Large SNP Data Sets, (Genetics) June 2014 197:573-589 [Genetics, Biorxiv]
downloaded from https://github.com/rajanil/fastStructure
"""


def plot_admixture(admixture, population_indices, population_labels, title, subset_k=None):
    N, K = admixture.shape
    if subset_k:
        colors = [colorsys.hsv_to_rgb(h, 0.9, 0.7) for h in np.linspace(0, 1, subset_k + 1)[:-1]]
    else:
        colors = [colorsys.hsv_to_rgb(h, 0.9, 0.7) for h in np.linspace(0, 1, K + 1)[:-1]]
    text_color = 'k'
    bg_color = 'w'
    fontsize = 12

    figure = plot.figure(figsize=(16, 10))

    xmin = 0.13
    ymin = 0.2
    height = 0.6
    width = 0.74
    indiv_width = width / N
    subplot = figure.add_axes([xmin, ymin, width, height])
    [spine.set_linewidth(0.001) for spine in list(subplot.spines.values())]

    for k in range(K):
        if k:
            bottoms = admixture[:, :k].sum(1)
        else:
            bottoms = np.zeros((N,), dtype=float)

        lefts = np.arange(N) * indiv_width
        subplot.bar(lefts, admixture[:, k], width=indiv_width, bottom=bottoms, facecolor=colors[k], edgecolor=colors[k],
                    linewidth=0.4)

        subplot.axis([0, N * indiv_width, 0, 1])
        subplot.tick_params(axis='both', top=False, right=False, left=False, bottom=False)
        xtick_labels = tuple(map(str, [''] * N))
        subplot.set_xticklabels(xtick_labels)
        ytick_labels = tuple(map(str, [''] * K))
        subplot.set_yticklabels(ytick_labels)

    position = subplot.get_position()
    title_position = (0.5, 0.9)
    figure.text(title_position[0], title_position[1], title, fontsize=fontsize, \
                color='k', horizontalalignment='center', verticalalignment='center')

    for p, popname in enumerate(population_labels):
        indices = np.where(population_indices == p)[0]
        if indices.size > 0:
            vline_pos = ((indices.max() + 1) * indiv_width) - (indiv_width / 2)
            # print(vline_pos)
            # vline_pos = (indices.max() ) * indiv_width
            subplot.axvline(vline_pos, linestyle='-', linewidth=0.2, c='#888888')
            label_position = (xmin + (2 * indices.min() + indices.size) * 0.5 * indiv_width, ymin - 0.01)
            figure.text(label_position[0], label_position[1], popname, fontsize=6, color='k', \
                        horizontalalignment='right', verticalalignment='top', rotation=70)

    return figure


def get_admixture_proportions(params):
    # load admixture proportions
    handle = open('%s.%d.Q' % (params['inputfile'], params['K']), 'r')
    admixture = np.array([line.strip().split() for line in handle]).astype('float')
    handle.close()
    N, K = admixture.shape
    admixture = admixture / admixture.sum(1).reshape(N, 1)

    # get population labels
    if 'popfile' in params:
        handle = open(params['popfile'], 'r')
        populations = [line.strip() for line in handle]
        handle.close()
        population_labels = list(set(populations))

        # group populations by cluster similarity
        population_cluster = [np.mean(admixture[[i for i, p in enumerate(populations) if p == label], :], 0).argmax() \
                              for label in population_labels]
        population_labels = [population_labels[j] for j in np.argsort(population_cluster)]
        population_indices = np.array([population_labels.index(pop) for pop in populations])

        # re-order samples in admixture matrix
        order = np.argsort(population_indices)
        population_indices = population_indices[order]
        admixture = admixture[order, :]

    else:

        print(
            "file with population labels is not provided or does not exist .... \ncreating population labels based on inferred admixture proportions")
        population_labels = ['population %d' % i for i in range(1, K + 1)]
        population_indices = np.argmax(admixture, 1)

        # re-order samples in admixture matrix
        order = np.argsort(population_indices)
        population_indices = population_indices[order]
        admixture = admixture[order, :]
        order = np.arange(N)
        for k in range(K):
            order[population_indices == k] = order[population_indices == k][
                np.argsort(admixture[population_indices == k, :][:, k])[::-1]]
        admixture = admixture[order, :]

    return admixture, population_indices, population_labels


def parseopts(opts):
    """
    parses the command-line flags and options passed to the script
    """

    params = {}
    params['subset'] = None
    params['subset_k'] = None

    for opt, arg in opts:

        if opt in ["-K"]:
            params['K'] = int(arg)

        elif opt in ["--input"]:
            params['inputfile'] = arg

        elif opt in ["--output"]:
            params['outputfile'] = arg

        elif opt in ["--popfile"]:
            params['popfile'] = arg

        elif opt in ["--title"]:
            params['title'] = arg

        elif opt in ["--subset"]:
            params['subset'] = arg

        elif opt in ["--subset_k"]:
            params['subset_k'] = arg
    return params


def usage():
    """
    brief description of various flags and options for this script
    """

    print("\nHere is how you can use this script\n")
    print("Usage: python %s" % sys.argv[0])
    print("\t -K <int>  (number of populations)")
    print("\t --input=<file>  (/path/to/input/file; same as output flag passed to structure.py)")
    print("\t --output=<file> (/path/to/output/file)")
    print("\t --popfile=<file> (file with known categorical labels; optional)")
    print("\t --title=<figure title> (a title for the figure; optional)")
    print("\t --subset= comma separated population name to subset; optional")
    print("\t --subset_k= K for subset pop; optional. if given subset, default is number of subset populations")

def string_param_2_array(paramsstr, type='float'):
    """
    Just to break the param string send by the shell script to numpy array
    :param paramsstr: the string froamt of the paramters. comma seprated '1122,112,.443..'
    :type to tell what will be the format for all the params. default is float
    :return: will send the numpy array [1122.0,112.0,0.443]
    """
    params = paramsstr.split(',')
    params = [param.strip() for param in params]
    params = list(map(eval(type), params))

    return params
def subsetting(admixture, population_indices, population_labels, subset, subset_k=None):
    """
    subsetting the whole admixture plot to better to look at important population
    Args:
        admixture: admixture coming from get_admixture_proportions
        population_indices: population indices of series coming from get_admixture_proportions
        population_labels: population names of those population indices
        subset: susbet population names divivded by comma
        subset_k: K values for subset populations

    Returns: will return sub-setted admixture, population_indices, population_labels and subset_k for required format

    """
    subset = MSPrimeClass.MsPrimeSFSCommandGenerator.string_param_2_array(subset, type='str')
    if subset_k:
        subset_k = int(subset_k)
    else:
        subset_k = len(subset)
    names = [population_labels[indice] for indice in population_indices]
    admixturedf = pandas.DataFrame(admixture, index=names)
    admixturedf = admixturedf.loc[subset, :]
    population_indices = pandas.Categorical(admixturedf.index).codes
    population_labels = dict(zip(pandas.Categorical(admixturedf.index).codes, admixturedf.index)).values()
    topcol = admixturedf.sum().sort_values(ascending=False).head(subset_k).index
    admixture = admixturedf.iloc[:, topcol].values
    return admixture, population_indices, population_labels, subset_k


if __name__ == "__main__":

    # parse command-line options
    argv = sys.argv[1:]
    smallflags = "K:"
    bigflags = ["input=", "output=", "popfile=", "title=", "subset=", "subset_k="]
    try:
        opts, args = getopt.getopt(argv, smallflags, bigflags)
        if not opts:
            usage()
            sys.exit(2)
    except getopt.GetoptError:
        print("Incorrect options passed")
        usage()
        sys.exit(2)
    print(opts)
    params = parseopts(opts)

    # get the data to be plotted
    admixture, population_indices, population_labels = get_admixture_proportions(params)
    if 'title' in params:
        title = params['title']
    else:
        title = params['inputfile']
    if params['subset']:
        admixture, population_indices, population_labels, params['subset_k'] = subsetting(admixture=admixture,
                                                                                          population_indices=population_indices,
                                                                                          population_labels=population_labels,
                                                                                          subset=params['subset'],
                                                                                          subset_k=params['subset_k'])
    # plot the data
    figure = plot_admixture(admixture, population_indices, population_labels, title, subset_k=params['subset_k'])
    figure.savefig(params['outputfile'], dpi=300)
