#!/usr/bin/env python3
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
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


def load_host_links(filename):
    print(f"Loading data from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[0] != "mobile_contig_name" or \
                        cut[5] != "cluster_name" or \
                        cut[14] != "mobile_element_copies_per_cell":
                    raise ImportError("not a good master table")
                continue
            virus = cut[0]
            bacteria = cut[5]
            copy_count = float(cut[14])
            adjusted_inter_vs_intra_ratio = float(cut[16])
            if adjusted_inter_vs_intra_ratio < 0.1:
                continue
            if virus not in data:
                data[virus] = list()
            data[virus].append(copy_count)
    return data


def load_validation(filename):
    print(f"Loading data from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                continue
            threshold = float(cut[0])
            support = float(cut[-1])
            data[threshold] = support
    return data


def categorize_data_for_roc_analysis(host_data, min_value, max_value, increment):
    """
    Go over and generate copy count threshold values and store the number of links kept at each cut off point
    Args:
        host_data (dict)
        min_value (float): lowest value to test as threshold
        max_value (float): highest value to test as threshold
        increment (float): multiplier to use for incrementing threshold
    Returns:
        roc_data (dict[float:list]): host hits remaining at each copy count threshold
        threshold_values (list[float]): copy count threshold generated

    """
    print("Generating ROC categories")
    roc_data = dict()
    threshold_values = list()
    threshold = min_value
    while threshold < max_value:
        roc_data[threshold] = [0, set()]
        threshold_values.append(threshold)
        threshold *= increment
    for mobile_contig, copy_counts in host_data.items():
        for copy_count in copy_counts:
            for threshold in threshold_values:
                if threshold > copy_count:
                    break
                roc_data[threshold][0] += 1
                roc_data[threshold][1].add(mobile_contig)
    for threshold in threshold_values:
        roc_data[threshold][1] = len(roc_data[threshold][1])
    return roc_data, threshold_values


def generate_roc_curve_values(roc_data, threshold_values, min_value):
    """
    Generate the x and y values of the copy count thresholding ROC curve
    Args:
        roc_data (dict[float:list]): host hits remaining at each copy count threshold
        threshold_values (list[float]): copy count threshold generated
        min_value (float): lowest value to test as threshold
    Returns:
        true_positives (list[float]): y-axis values of ROC curve
        false_positives: (list[float]): x-axis values of ROC curve

    """
    print("Generating ROC curve")
    true_positives = list()
    false_positives = list()
    max_possible_hits_accepted = roc_data[min_value][0]
    max_possible_mobile_contigs_with_hosts = roc_data[min_value][1]
    for threshold in threshold_values:
        hits_accepted = roc_data[threshold][0]
        false_positives.append(hits_accepted / max_possible_hits_accepted)
        mobile_contigs_with_hosts = roc_data[threshold][1]
        true_positives.append(mobile_contigs_with_hosts / max_possible_mobile_contigs_with_hosts)
    return true_positives, false_positives


def calculate_area_under_curve(false_positives, true_positives):
    """
    Calculate the area under a ROC curve
    Args:
        false_positives (list): x-values of a ROC curve
        true_positives (list): y-values of a ROC curve
    Returns:
        auc (float): area under the curve

    """
    print("Calculating AUC")
    auc = 0
    for i, false_positive in enumerate(false_positives):
        if i + 1 == len(false_positives):
            break
        x_delta = false_positives[i] - false_positives[i + 1]
        height = (true_positives[i] + true_positives[i + 1]) / 2
        auc += x_delta * height
    print(f"ROC area under curve = {auc}")
    return auc


def get_optimal_threshold(threshold_values, true_positives, false_positives, min_fraction_without_hosts=0.8):
    """
    Run down the ROC curve and stop when the optimal threshold has been reached
    Args:
        threshold_values (list[float]): copy count threshold generated
        true_positives (list[float]): y-axis values of ROC curve
        false_positives: (list[float]): x-axis values of ROC curve
        min_fraction_without_hosts (float): lowest acceptable fraction of mobile elements still with a host
    Returns:
        optimal_threshold (float): the optimal copy count value
        fp_rate (float): the x-value on the ROC curve of the optimal cut-off
        tp_rate (float): the y-value on the ROC curve of the optimal cut-off

    """
    print("Calculating optimal threshold")
    optimal_threshold = 0
    fp_rate = 0
    tp_rate = 0
    for i, threshold in enumerate(threshold_values):
        optimal_threshold = threshold
        fp_rate = false_positives[i]
        tp_rate = true_positives[i]
        if fp_rate + tp_rate < 1 or tp_rate < min_fraction_without_hosts:
            break
    print(f"Optimal value = {optimal_threshold}")
    return optimal_threshold, fp_rate, tp_rate


def plot_roc_curve(ax, true_positives, false_positives, fp_rate=1, tp_rate=1, optimal_threshold=None, auc=None):
    print(f"Drawing ROC  curve")
    ax.plot(false_positives, true_positives, c="k", alpha=0.3)
    ax.scatter(false_positives, true_positives, c="k", alpha=0.6)
    ax.scatter([fp_rate], [tp_rate], c="r")
    ax.plot([0, 1], [1, 0], "--", c="r", alpha=0.5)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlabel("Fraction of all hits kept")
    ax.set_ylabel("Fraction of viruses with one host")
    if auc is not None:
        auc = round(auc, 2)
    optimal_threshold = round(optimal_threshold, 2)
    ax.text(0.25, 0.78, f"Chosen threshold\n(0.14 copies per cell)", c="r")
    ax.grid(axis="both", ls="--", alpha=0.1, c="k", which="major")


def plot_host_validation(ax, prophage_validation, color="k", label=None):
    xs = list()
    ys = list()
    chosen_validation = 0
    for x, y in prophage_validation.items():
        xs.append(x)
        ys.append(y)
        if x == 0.14:
            chosen_validation = y
    ys = [x for y, x in sorted(zip(xs, ys), reverse=True)]
    xs.sort(reverse=True)
    ax.plot(xs, ys, c=color, alpha=0.3)
    ax.axvline(0.14, linestyle="--", c="r", alpha=0.4)
    ax.scatter(xs, ys, c=color, alpha=0.6, label=label)
    ax.scatter([0.14], [chosen_validation], c="r")

    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-1, 101)
    ax.set_xlabel("Minimum copy count threshold")
    ax.set_ylabel("Percent validated prophage hosts")
    ax.text(0.17, 70, f"Chosen threshold\n(0.14 copies per cell)", c="r")
    ax.grid(axis="both", ls="--", alpha=0.1, c="k", which="major")



def load_rarefaction_data(directory):
    master_data = dict()
    for filename in os.listdir(directory):
        path = os.path.join(directory, filename)
        read_ct = int(filename.split(".")[0].split("_")[-1])
        data = load_host_links_all(path)
        master_data[read_ct] = data
    return master_data


def load_host_links_all(filename):
    print(f"Loading data from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[0] != "mobile_contig_name" or \
                        cut[5] != "cluster_name" or \
                        cut[14] != "mobile_element_copies_per_cell":
                    raise ImportError("not a good master table")
                continue
            virus = cut[0]
            viral_depth = float(cut[3])
            bacteria = cut[5]
            intra_links = float(cut[9])
            inter_links = float(cut[11])
            if inter_links < 5 or intra_links < 10:
                continue
            links = int(cut[11])
            host_depth = float(cut[8])
            norm_links = (links / host_depth) / viral_depth
            copy_count = float(cut[14])
            adjusted_inter_vs_intra_ratio = float(cut[16])
            if virus not in data:
                data[virus] = dict()
            data[virus][bacteria] = (links, norm_links, copy_count, adjusted_inter_vs_intra_ratio)
    return data


def get_contig_host_pool(master_data, min_count=4, max_to_plot=10):
    contig_host_pool = dict()
    for data in master_data.values():
        for virus, subdata in data.items():
            for host in subdata:
                contig_host = f"{virus}:{host}"
                if contig_host not in contig_host_pool:
                    contig_host_pool[contig_host] = 0
                contig_host_pool[contig_host] += 1
    filtered_contig_host_pool = set()
    for contig_host, count in contig_host_pool.items():
        if len(filtered_contig_host_pool) > max_to_plot:
            break
        if count >= min_count:
            filtered_contig_host_pool.add(contig_host)
    print(f"Kept {len(filtered_contig_host_pool)} out of {len(contig_host_pool)} virus:host links for plotting")
    return filtered_contig_host_pool


def plot_data(contig_host_pool, master_data, n, ax):
    print(f"plotting value {n}")
    for i, contig_host in enumerate(contig_host_pool):
        if i % 1000 == 0 and i > 0:
            print(f"processed {i} links")
        virus = contig_host.split(":")[0]
        host = contig_host.split(":")[1]
        xs = list()
        ys = list()
        for read_ct, data in master_data.items():
            if virus not in data:
                continue
            if host not in data[virus]:
                continue
            value = data[virus][host][n]
            xs.append(read_ct)
            ys.append(value)
        ys = [x for y, x in sorted(zip(xs, ys), reverse=True)]
        xs.sort(reverse=True)
        ax.plot(xs, ys, "-", c="k", alpha=0.03)
        # ax.plot(xs, ys, "o", c="k", alpha=0.05, markersize=1)


def make_histogram_table(contig_host_pool, master_data, n):
    print(f"making histogram data for value {n}")
    bins = 50
    st = 0
    fi = 2
    bin_size = (fi - st) / bins

    xs = list()
    ys = list()
    x = st
    while x <= fi:
        read_ct = max(master_data.keys())
        data = master_data[read_ct]
        count = 0
        for i, contig_host in enumerate(contig_host_pool):
            if i % 1000 == 0 and i > 0:
                print(f"processed {i} links")
            virus = contig_host.split(":")[0]
            host = contig_host.split(":")[1]
            value = data[virus][host][n]
            if x <= value < x + bin_size:
                count += 1
        xs.append(x)
        ys.append(count)
        x += bin_size
    for i, x in enumerate(xs):
        y = ys[i]
        print(x, y, sep="\t")


def main():
    master_data = load_rarefaction_data("rarefaction_data")
    contig_host_pool = get_contig_host_pool(master_data, min_count=1, max_to_plot=10000)

    host_data = load_host_links("unfiltered_master_table.tsv")
    all_prophage_validation = load_validation("validation_values_all_prophages.tsv")
    single_host_prophage_validation = load_validation("validation_values_singlehost_prophages.tsv")

    min_value = 0.0001
    roc_data, threshold_values = categorize_data_for_roc_analysis(host_data, min_value, 1000, 1.1)
    true_positives, false_positives = generate_roc_curve_values(roc_data, threshold_values, min_value)
    auc = calculate_area_under_curve(false_positives, true_positives)
    optimal_threshold, fp_rate, tp_rate = get_optimal_threshold(threshold_values, true_positives, false_positives)

    print("plotting subplots")
    fig = plt.figure(figsize=(15, 5.4))
    initialize_plot_style()
    gs = gridspec.GridSpec(nrows=2, ncols=3, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    plt.gca().invert_xaxis()
    ax2 = fig.add_subplot(gs[1, 0])
    plt.gca().invert_xaxis()
    ax3 = fig.add_subplot(gs[:, 1])
    ax4 = fig.add_subplot(gs[:, 2])

    #########################
    for i, ax in enumerate([ax1, ax2]):
        ax.set_xscale("log")
        plot_data(contig_host_pool, master_data, i + 2, ax)
        ax.set_xlim(105000000, 95000)
        ax.grid(axis="both", ls="--", alpha=0.1, c="k", which="major")
        ax.set_ylim(-0.05, 3.05)
    ax2.set_xlabel("Hi-C library size (read count)")
    ax1.set_ylabel("Copies per cell")
    ax2.set_ylabel("Connectivity ratio (R')")

    #########################
    plot_roc_curve(ax3, true_positives, false_positives, fp_rate=fp_rate, tp_rate=tp_rate,
                   optimal_threshold=optimal_threshold, auc=auc)
    plot_host_validation(ax4, all_prophage_validation, label="All prophages")
    plot_host_validation(ax4, single_host_prophage_validation, color="b", label="Single-host prophages")

    ax1.set_title("A", fontsize=20, x=-0.1)
    ax2.set_title("B", fontsize=20, x=-0.1, y=1.05)
    ax3.set_title("C", fontsize=20, x=-0.1)
    ax4.set_title("D", fontsize=20, x=-0.1)
    plt.tight_layout()
    plt.savefig("figure.png", dpi=300)


main()
