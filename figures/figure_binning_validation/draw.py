#!/usr/bin/env python3
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
    completions = list()
    contaminations = list()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[9] != "completeness":
                    raise ImportError("not a checkv data file")
                continue
            name = cut[0]
            completion = cut[9]
            if completion == "NA":
                completion = 0
            else:
                completion = float(completion)
            completions.append(completion)
            contaminations.append(0)
    completions.sort(reverse=True)
    return completions, contaminations


def load_hifi_data(filename):
    print(f"Loading data from file {filename}")
    completions = list()
    contaminations = list()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[4] != "completion" or cut[5] != "contamination":
                    raise ImportError("not a long-read validation data file")
                continue
            name = cut[0]
            completion = float(cut[4])
            contamination = float(cut[5])
            completions.append(completion)
            contaminations.append(contamination)
    contaminations = [x for y, x in sorted(zip(completions, contaminations), reverse=True)]
    completions = completions.copy()
    completions.sort(reverse=True)
    return completions, contaminations


def get_max_x(comp1, comp2):
    max_1 = 0
    for i in comp1:
        if i >= 20:
            max_1 += 1
    max_2 = 0
    for i in comp2:
        if i >= 20:
            max_2 += 1
    max_x = max(max_1, max_2)
    return max_x


def plot_lines(comp1, comp2, cont1, cont2, ax1):
    max_x = get_max_x(comp1, comp2)
    ax1.set_xlim(-1, max_x)
    #ax2 = ax1.twinx()
    plot_completion(comp1, ax1, color="grey", label="Contig completion")
    plot_completion(comp2, ax1, color="k", label="vMAG completion")
    plot_contamination(cont1, ax1)
    plot_contamination(cont2, ax1, label="vMAG contamination")
    ax1.set_ylim(-1, 101)
    #ax2.set_ylim(-1, 101)
    ax1.grid(axis="both", ls="--", alpha=GRID_ALPHA, c=MAIN_COLOR)
    #ax2.grid(axis="both", ls="--", alpha=GRID_ALPHA, c=MAIN_COLOR)


def plot_completion(data, ax, color='k', label=None):
    """
    Plot the descending completion values of clusters produced by one binning method (as points connected with a line)
    Args:
        data (list): list of completion values to be plotted, sorted in descending numerical order
        ax (axis): plot axis to plot on
    Returns:
        num_data (int): number of data points plotted

    """
    ax.plot(data, linewidth=1, color=color, alpha=1)
    ax.plot(data, "o", color=color, alpha=MARKER_ALPHA, markersize=MARKER_SIZE, label=label)
    ax.set_yticks(PERCENTILE_TICKS)
    ax.set_yticklabels(PERCENTILE_TICKS, color=MAIN_COLOR)
    ax.tick_params(axis='y', colors=MAIN_COLOR)
    ax.grid(axis="both", ls="--", alpha=GRID_ALPHA, c=MAIN_COLOR)
    return len(data)


def plot_contamination(data, ax, label=None):
    """
    Plot the contamination values of clusters produced by one binning method (as points)

    Args:
        data (list): list of contamination values to be plotted, sorted by descending numerical order of the
                     corresponding cluster completion values
        ax (axis): plot axis to plot on
    Returns:
        num_data (int): number of data points plotted
    """

    if sum(data) > 0:
        ax.plot(data, "o", color=SECONDARY_COLOR, alpha=MARKER_ALPHA, markersize=MARKER_SIZE, label=label)
    #ax.set_yticks(PERCENTILE_TICKS)
    #ax.set_yticklabels(PERCENTILE_TICKS, color=SECONDARY_COLOR)
    ax.grid(axis="both", ls="--", alpha=GRID_ALPHA, c=MAIN_COLOR)
    return len(data)


def plot_barplot(completions, names_list, ax):
    data = generate_categories(completions)
    patch_handles, patch_handle_labels = draw_barplot(data, names_list, ax)
    add_labels(patch_handles, patch_handle_labels, ax)
    ax.set_xlabel("Number of viral genomes", fontsize=AXIS_FONT)

    # remove borders
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_linewidth(MAIN_LINEWIDTH / 2)
    ax.spines['bottom'].set_color('black')
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(MAIN_LINEWIDTH / 2)
    ax.spines['left'].set_color('black')
    plt.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=True, left=False, right=False,
                    labelleft=True)
    ax.grid(axis="x", alpha=GRID_ALPHA, ls="--", c=MAIN_COLOR)
    format_barplot(ax)


def draw_barplot(data, binners, ax):
    """
    Draw barplot of number of viral mags at each threshold in each set being compared
    Args:
        data (dict[dict]): tally for number of good clusters at each threshold cut off in each bin/sequence set
        binners (list): names of the viral sets being compared
        ax (plt axis): axis to draw barplot on

    Returns:
        patch_handles (list): handles pointing to the barplot element locations
        patch_handle_labels (list): handles pointing to the barplot labels

    """
    initialize_plot_style()
    width = BAR_WIDTH
    prev = None
    patch_handles = list()
    patch_handle_labels = list()
    N = len(binners)
    ind = np.arange(N)[::-1]
    categories = sorted(data.keys(), reverse=True)
    for i, category in enumerate(categories):
        bars = data[category][:]
        patch_handles.append(ax.barh(ind, bars, width, left=prev, label=category, color=COLOR_CATEGORIES[i],
                                     edgecolor="w", tick_label=binners))
        if prev is None:
            prev = bars[:]
            labels = bars[:]
        else:
            for j, val in enumerate(bars):
                labels[j] -= val
                val = prev[j] + val
                prev[j] = val
        patch_handle_labels.append(bars)
    return patch_handles, patch_handle_labels


def add_labels(patch_handles, patch_handle_labels, ax, units=""):
    """
    Add numerical labels inside the barplot elements to make it easier to interpret
    Args:
        patch_handles (list): handles pointing to the barplot element locations
        patch_handle_labels (list): handles pointing to the barplot labels
        ax (axis): plot axis to plot on
        units (str): and units to append to labels [default: ""]
    Returns:
        created_labels (list): string labels added to barplot

    """
    ax.legend(bbox_to_anchor=(0.99, -0.13), ncol=len(PERCENTILE_CATEGORIES), facecolor=None,
               prop={'size': MAIN_FONT}, frameon=True, fontsize=MAIN_FONT,
               title="Categories of viral genome completion", columnspacing=0.3)
    created_labels = list()
    # add numerical labels
    cumulative_width = list()
    for _ in patch_handles[0]:
        cumulative_width.append(0)
    for i, handle in enumerate(patch_handles):
        for j, patch in enumerate(handle.get_children()):
            bl = patch.get_xy()
            y = 0.43 * patch.get_height() + bl[1]
            x = 1.0 * patch.get_width() / 2 + cumulative_width[j]
            cumulative_width[j] += patch_handle_labels[i][j]
            if patch_handle_labels[i][j] >= 5:
                label = str(int(patch_handle_labels[i][j])) + units
                ax.text(x, y, label, ha='center', color='w', fontsize=MAIN_FONT)
                created_labels.append(label)
    return created_labels


def format_barplot(ax):
    """
    Make the barplot more aesthetically-pleasing by making it more minimalistic.
    Args:
        ax (axis): axis on which the barplot is plotted
    Returns: None

    """
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_linewidth(MAIN_LINEWIDTH / 2)
    ax.spines['bottom'].set_color('black')
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(MAIN_LINEWIDTH / 2)
    ax.spines['left'].set_color('black')
    plt.tick_params(axis="both", which="both", bottom=False, top=False, labelbottom=True, left=False, right=False,
                    labelleft=True)
    ax.grid(axis="x", alpha=GRID_ALPHA, ls="--", c=MAIN_COLOR)
    plt.tight_layout(w_pad=0.5)


def generate_categories(completion_data, thresholds=PERCENTILE_CATEGORIES):
    """
    Calculate the number of vMAGs/sequences at certain completion/contamination quality thresholds
    Args:
        completion_data (list[list]): completion values of sequences in each viral bin/sequence set

    Returns:
        data (dict[dict]): tally for number of good clusters at each thresholf cut off in each bin/sequence set

    """
    data = dict()
    for percentile in thresholds:
        name = f">{percentile}%"
        data[name] = list()

    for completions in completion_data:
        counts = dict()
        for percentile in thresholds:
            counts[percentile] = 0
        for comp in completions:
            for percentile in sorted(thresholds, reverse=True):
                if comp >= percentile:
                    counts[percentile] += 1
                    break
        for percentile, count in counts.items():
            name = f">{percentile}%"
            data[name].append(count)
    return data


def main():
    cv_contig_comp, cv_contig_cont = load_checkv_data("contigs_CV.tsv")
    cv_cluster_comp, cv_cluster_cont = load_checkv_data("clusters_CV.tsv")
    lr_contig_comp, lr_contig_cont = load_hifi_data("contigs_LR.tsv")
    lr_cluster_comp, lr_cluster_cont = load_hifi_data("clusters_LR.tsv")

    print("plotting subplots")
    fig = plt.figure(figsize=(12, 6))
    initialize_plot_style()
    gs = gridspec.GridSpec(nrows=2, ncols=2, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[:, 1])

    plot_lines(cv_contig_comp, cv_cluster_comp, cv_contig_cont, cv_cluster_cont, ax1)
    plot_lines(lr_contig_comp, lr_cluster_comp, lr_contig_cont, lr_cluster_cont, ax2)
    ax1.set_ylabel("CheckV metrics", fontsize=AXIS_FONT)
    ax2.set_ylabel("HiFi metrics", fontsize=AXIS_FONT)
    ax2.set_xlabel("Viral genome count", fontsize=AXIS_FONT)
    lgnd = ax2.legend(bbox_to_anchor=(1.05, -0.3), ncol=2, facecolor=None,
                      prop={'size': 12}, frameon=True, fontsize=MAIN_FONT, columnspacing=0.3)
    lgnd.legendHandles[0]._legmarker.set_markersize(8)
    lgnd.legendHandles[1]._legmarker.set_markersize(8)
    lgnd.legendHandles[2]._legmarker.set_markersize(8)
    plot_barplot([cv_contig_comp, cv_cluster_comp, lr_contig_comp, lr_cluster_comp],
                 ["Contigs\n(CheckV)", "vMAGs\n(CheckV)", "Contigs\n(HiFi)", "vMAGs\n(HiFi)", ], ax3)

    ax1.set_title("A", fontsize=20, x=-0.15, y=0.88)
    ax2.set_title("B", fontsize=20, x=-0.15, y=0.88)
    ax3.set_title("C", fontsize=20, x=-0.18, y=0.95)


    plt.tight_layout()
    plt.savefig("figure.png", dpi=300)


main()
