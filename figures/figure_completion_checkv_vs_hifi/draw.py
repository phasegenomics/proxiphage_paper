#!/usr/bin/env python3
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import numpy as np


# plotting constants
PERCENTILE_CATEGORIES = [20, 50, 70, 90]
COLOR_CATEGORIES = ["#000033", "#000099", "#0000ff", "#6666ff"]
PERCENTILE_TICKS = [0, 20, 40, 60, 80, 100]
GRID_ALPHA = 0.2
LINE_ALPHA = 0.1
MARKER_ALPHA = 0.6
TITLE_FONT = 16
AXIS_FONT = 14
MAIN_FONT = 12
LABEL_FONT = 10
MAIN_COLOR = 'k'
SECONDARY_COLOR = "r"
MAIN_LINEWIDTH = 1
MARKER_SIZE = 2
BAR_WIDTH = 0.85
HISTOGRAM_BINS = 50
PNG_DPI = 300


def initialize_plot_style():
    """
    Set matplotlib visual style
    Returns: None

    """
    matplotlib.style.use('ggplot')
    plt.rcParams['lines.linewidth'] = MAIN_LINEWIDTH
    plt.rcParams['axes.facecolor'] = 'w'
    plt.rcParams['xtick.color'] = MAIN_COLOR
    plt.rc('xtick', labelsize=LABEL_FONT)
    plt.rc('ytick', labelsize=LABEL_FONT)


def load_checkv_data(filename):
    print(f"Loading data from file {filename}")
    completions = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[9] != "completeness":
                    raise ImportError("not a checkv data file")
                continue
            name = cut[0].split("|")[0]
            completion = cut[9]
            if completion != "NA":
                completion = float(completion)
                completions[name] = completion
    return completions


def load_hifi_data(filename):
    print(f"Loading data from file {filename}")
    completions = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[4] != "completion" or cut[5] != "contamination":
                    raise ImportError("not a long-read validation data file")
                continue
            name = cut[0]
            completion = float(cut[4])
            if completion > 0:
                completions[name] = completion
    return completions


def pair_data(cv_comp, lr_comp):
    xs = list()
    ys = list()
    for name, x in lr_comp.items():
        if name in cv_comp:
            y = cv_comp[name]
            xs.append(x)
            ys.append(y)
    return xs, ys


def main():
    cv_contig_comp = load_checkv_data("contigs_CV.tsv")
    cv_cluster_comp = load_checkv_data("clusters_CV.tsv")
    lr_contig_comp = load_hifi_data("contigs_LR.tsv")
    lr_cluster_comp = load_hifi_data("clusters_LR.tsv")
    initialize_plot_style()
    contig_xs, contig_ys = pair_data(cv_contig_comp, lr_contig_comp)
    cluster_xs, cluster_ys = pair_data(cv_cluster_comp, lr_cluster_comp)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(contig_xs, contig_ys, c="k", alpha=0.3, label="Viral contigs")
    ax.scatter(cluster_xs, cluster_ys, c="r", alpha=0.3, label="Viral MAGs")
    ax.grid(axis="both", ls="-", alpha=0.1, c="k", which="major")
    ax.set_xlim(-1, 101)
    ax.set_ylim(-1, 101)
    ax.set_xlabel("Completion estimated with HiFi reference (%)", fontsize=16)
    ax.set_ylabel("Completion estimated with CheckV (%)", fontsize=16)
    ax.legend(fontsize=14)

    plt.tight_layout()
    plt.savefig("figure.png", dpi=300)


main()
