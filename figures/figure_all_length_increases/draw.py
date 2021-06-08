#!/usr/bin/env python3
# Plot histogram of viral contigs/vMAGs lengths
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.gridspec as gridspec


# plotting constants
PERCENTILE_CATEGORIES = [20, 50, 70, 90]
COLOR_CATEGORIES = ["#000033", "#000099", "#0000ff", "#6666ff"]
PERCENTILE_TICKS = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240]
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


def load_checkv(filename):
    print(f"Loading data from file {filename}")
    comps = list()
    lengths = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            if i > 0:
                cut = line.strip().split("\t")
                name = cut[0]
                length = int(cut[1])
                comp = cut[9]
                if comp == "NA":
                    comp = 0
                else:
                    comp = float(comp)
                lengths[name] = length
                comps.append(comp)
    comps.sort(reverse=True)
    return lengths, comps


def initialize_plot_style():
    """
    Set matplotlib visual style
    Returns: None

    """
    matplotlib.style.use('ggplot')
    plt.rcParams['lines.linewidth'] = 1
    plt.rcParams['axes.facecolor'] = 'w'
    plt.rcParams['xtick.color'] = 'k'
    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)


def plot_barplot(completions, names_list, ax):
    data = generate_categories(completions)
    patch_handles, patch_handle_labels = draw_barplot(data, names_list, ax)
    add_labels(patch_handles, patch_handle_labels, ax)
    ax.set_xlabel(None, fontsize=AXIS_FONT)

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
    print("draw barplot")
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
            if patch_handle_labels[i][j] != 0:
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


def plot_violin(data, sample_names, ax):
    print("plotting violin plots")
    df = pd.DataFrame(columns=['Sample', "ID", 'Length'])
    for i, lengths in enumerate(data):
        sample = sample_names[i]
        id = sample
        for l in lengths.values():
            subdf = pd.DataFrame({'Sample': [sample], "ID": [id], "Length": [l / 1000]})
            df = df.append(subdf, ignore_index=True)

    sns.violinplot(y="ID", x="Length", data=df, linewidth=0, ax=ax, width=1)
    sns.swarmplot(y="ID", x="Length", data=df, size=1.5, edgecolor='black', color="k", linewidth=0, ax=ax)
    ax.set(ylabel=None, xlabel=None)
    ax.grid(axis="x", alpha=GRID_ALPHA, ls="--", c=MAIN_COLOR)


def main():
    fig = plt.figure(figsize=(6, 10))
    initialize_plot_style()
    gs = gridspec.GridSpec(nrows=4, ncols=1, figure=fig)
    axes = list()
    for i in range(0, 4):
        axes.append(fig.add_subplot(gs[i]))

    samples = ["sheep_gut", "human_gut", "rumen", "wastewater"]
    sample_names = ["Sheep fecal", "Human fecal", "Cow rumen", "Wastewater"]
    for i, sample in enumerate(samples):
        vmags_file = f"{sample}.vmags.checkv.tsv"
        contigs_file = f"{sample}.contigs.checkv.tsv"
        vmags_lengths, vmags_comps = load_checkv(vmags_file)
        contigs_lengths, contigs_comps = load_checkv(contigs_file)
        ax1 = axes[i]
        names = [f"{sample_names[i]}\nViral contigs", f"{sample_names[i]}\nViral MAGs"]
        plot_violin([contigs_lengths, vmags_lengths], names, ax1)

    # EXTRA FORMATTING
    labels = ["A", "B", "C", "D"]
    for i, ax in enumerate(axes):
        x = -0.25
        ax.set_title(labels[i], fontsize=20, x=x)

    axes[3].set_xlabel('Length (Kbp)', fontsize=14)


    plt.tight_layout()
    plt.savefig("figure.png", bbox_inches='tight', dpi=300)


main()
