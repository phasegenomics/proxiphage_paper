#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import random
import pandas
import numpy
import seaborn
from PIL import Image
import os


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


def load_mags(filename):
    print(f"Loading data from {filename}")
    data = dict()
    with open(filename) as input:
        for line in input:
            cut = line.strip().split("\t")
            contig = cut[0]
            mag = cut[1]
            data[contig] = mag
    return data


def load_host_links(filename):
    print(f"Loading data from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[0] != "mobile_cluster_name" or \
                        cut[5] != "cluster_name" or \
                        cut[14] != "mobile_element_copies_per_cell":
                    raise ImportError("not a good master table")
                continue
            virus = cut[0]
            bacteria = cut[5]
            copy_count = float(cut[14])
            if virus not in data:
                data[virus] = dict()
            data[virus][bacteria] = copy_count
    return data


def generate_column_color_labels(contig_mags, df):
    if contig_mags is None:
        col_colors = None
    else:
        bins = dict()
        colors = dict()
        mobile_clusters = dict()
        for contig, mag in contig_mags.items():
            if mag not in mobile_clusters:
                mobile_clusters[mag] = set()
            mobile_clusters[mag].add(contig)
        for cluster, contigs in mobile_clusters.items():
            represented_contigs = 0
            for contig in contigs:
                bins[contig] = cluster
                if contig in df:
                    represented_contigs += 1
            if represented_contigs > 1:
                random_color = (random.random(), random.random(), random.random())
                colors[cluster] = random_color
        col_colors = list()
        boring = True
        for mobile_element in df:
            if mobile_element in bins:
                cluster = bins[mobile_element]
                if cluster in colors:
                    col_colors.append(colors[cluster])
                    boring = False
                else:
                    col_colors.append((1, 1, 1))
            else:
                col_colors.append((1, 1, 1))
        if boring:
            col_colors = None
    return col_colors


def make_heatmap(data_dict, mobile_clusters=None, output_file="heatmap.png", x_label="Mobile elements",
                 y_label="MAG clusters"):
    """
    Make a clustered heatmap plot of mobile element to cluster HiC connectivity data
    Args:

        data_dict (dict[str[dict[str:float]]]): mobile elements pointing to clusters pointing to connectivity values
        output_file (str): path to output png file
        x_label (str): x-axix label string
        y_label (str): y-axix label string
        title (str): heatmap title
    Returns: None

    """
    print(f"Plotting heatmap {output_file}")
    df = pandas.DataFrame.from_dict(data_dict)
    df = df.fillna(0)
    # remove columns and rows with all 0s
    df = df.loc[:, (df != 0).any(axis=0)]
    df = df.loc[(df != 0).any(axis=1)]
    # log normalize
    df += 0.01
    df = numpy.log(df)

    rows, columns = df.shape
    print(f"Matrix size: {rows} X {columns}")
    xticklabels = False
    yticklabels = False
    col_colors = generate_column_color_labels(mobile_clusters, df)
    if col_colors is not None:
        print(f"Generated {len(col_colors)} randomized color labels for mobile element clusters")
    seaborn.set(font_scale=0.5)
    g = seaborn.clustermap(df, figsize=(6, 6), cmap="Greys_r", col_colors=col_colors,
                           xticklabels=xticklabels, yticklabels=yticklabels,
                           dendrogram_ratio=[0.15, 0.15], cbar_pos=[0.02, 0.85, 0.03, 0.1])
    ax = g.ax_heatmap
    g.ax_row_dendrogram.text(g.ax_row_dendrogram.get_xlim()[0] * 1,
                             g.ax_row_dendrogram.get_ylim()[0] / 2,
                             y_label, rotation=90, verticalalignment="center", fontsize=16)
    g.ax_col_dendrogram.text(g.ax_col_dendrogram.get_xlim()[1] / 2,
                             g.ax_col_dendrogram.get_ylim()[1] * 0.75,
                             x_label, rotation=0, horizontalalignment="center", fontsize=16)
    ax_color = g.ax_cbar
    ax_color.set_ylabel(f"Log10(copy_count)", fontsize=8)
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()


def manual_png(figure, x, y, fig, resize=False):
    path = os.path.realpath(figure)
    print(f"Inserting {path}")
    im = Image.open(path)
    if resize:
        plot_h = fig.bbox.ymax * 1
        plot_w = fig.bbox.xmax * 0.8
        w, h = im.size
        resize_ratio = min(plot_h / h, plot_w / w)
        h *= resize_ratio
        im = im.resize((int(w), int(h)), Image.ANTIALIAS)
    im = numpy.array(im).astype(numpy.float) / 255
    plt.figimage(im, x, y)


def main():
    hosts = load_host_links("unfiltered_clustered_master_table.tsv")
    initialize_plot_style()
    make_heatmap(hosts, output_file="unfiltered_clustered_master_table.png",
                 x_label="Viral MAGs (unfiltered hits)", y_label="Prokaryotic MAGs")

    hosts = load_host_links("filtered_clustered_master_table.tsv")
    initialize_plot_style()
    make_heatmap(hosts, output_file="filtered_clustered_master_table.png",
                 x_label="Viral MAGs (filtered hits)", y_label="Prokaryotic MAGs")
    plt.close()

    print("Make final merged figure")
    initialize_plot_style()
    fig = plt.figure(figsize=(13, 6))
    ax = fig.add_axes([0, 0, 1, 1])
    x_max = fig.bbox.xmax
    y_max = fig.bbox.ymax
    print(f"x_max: {x_max}, y_max: {y_max}")
    manual_png("unfiltered_clustered_master_table.png", 100, 0, fig)
    manual_png("filtered_clustered_master_table.png", 2000, 0, fig)
    ax.text(0.005, 0.93, "A", fontsize=30, zorder=2)
    ax.text(0.49, 0.93, "B", fontsize=30, zorder=2)
    ax.axis("off")

    plt.savefig("figure.png", dpi=300)


main()
