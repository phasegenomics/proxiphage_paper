#!/usr/bin/env python3
import sys
import numpy
import seaborn
import scipy
import matplotlib
import matplotlib.pyplot as plt
import pandas
import random


def load_checkm(filename):
    print(f"loading from {filename}")
    completions = dict()
    contaminations = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                continue
            mag = cut[0]
            comp = float(cut[11])
            cont = float(cut[12])
            completions[mag] = comp
            contaminations[mag] = cont
    return completions, contaminations


def load_hifi_validation(filename):
    print(f"loading from {filename}")
    completions = dict()
    contaminations = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                continue
            name = cut[0]
            comp = float(cut[4])
            cont = float(cut[5])
            completions[name] = comp
            contaminations[name] = cont
    return completions, contaminations


def select_hosts(possible_hosts, host_completions, host_contaminations, min_comp=20, max_cont=10):
    selected = set()
    for host in possible_hosts:
        if host_completions[host] > min_comp and host_contaminations[host] < max_cont:
            selected.add(host)
    return selected


def load_host_links(filename):
    print(f"Loading data from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[5] != "cluster_name" or \
                        cut[14] != "mobile_element_copies_per_cell":
                    raise ImportError(f"not a good master table: {cut[0], cut[5], cut[14]}")
                continue
            virus = cut[0]
            bacteria = cut[5]
            copy_count = float(cut[14])
            if virus not in data:
                data[virus] = dict()
            data[virus][bacteria] = copy_count
    return data


def select_vmags(possible_vmags, vmag_completions, vmag_contaminations, min_comp=0, max_cont=10):
    selected = set()
    for virus in possible_vmags:
        if vmag_completions[virus] > min_comp and vmag_contaminations[virus] < max_cont:
            selected.add(virus)
    return selected


def make_matrix(possible_hosts, vmags_contigs, host_link_data):
    data = dict()
    for vmag, contigs in vmags_contigs.items():
        for virus in contigs:
            if virus in host_link_data:
                data[vmag+":"+virus] = dict()
                for host in host_link_data[virus]:
                    if host in possible_hosts:
                        data[vmag+":"+virus][host] = host_link_data[virus][host]
    return data


def make_heatmap(data_dict, mobile_clusters=None, output_file="heatmap.png", x_label="Viruses",
                 y_label="MAG clusters", title="Matrix values"):
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
    df = pandas.DataFrame.from_dict(data_dict)
    df = df.fillna(0)
    # remove columns and rows with all 0s
    #df = df.loc[:, (df != 0).any(axis=0)]
    #df = df.loc[(df != 0).any(axis=1)]
    # log normalize
    df += 0.01
    df = numpy.log(df)

    rows, columns = df.shape
    if rows <= 50 and columns <= 50:
        xticklabels = True
        yticklabels = True
    else:
        xticklabels = False
        yticklabels = False
    col_colors = generate_column_color_labels(mobile_clusters, df)
    print(f"Plotting the heatmap data: {rows} X {columns}")
    g = seaborn.clustermap(df, figsize=(6, 6), cmap="Greys_r", col_colors=col_colors,
                           xticklabels=xticklabels, yticklabels=yticklabels,
                           col_cluster=False)
    ax = g.ax_heatmap
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title, fontsize=16, pad=100)
    ax_color = g.ax_cbar
    # ax_color.set_ylabel(f"Log10({title.split()[-1]})", fontsize=6)
    ax_color.set_ylabel(f"Log10", fontsize=8)
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()


def generate_column_color_labels(mobile_clusters_orig, df):
    """
    Generate color labels for heatmap columns to represent the mobile element clusters that the contigs are in
    Args:
        df (Pandas data frame): matrix being plotted in heatmap
    Returns:
        col_colors (list(tuple(float))): colors in order of column names

    """
    mobile_clusters = dict()
    for mag, contigs in mobile_clusters_orig.items():
        mobile_clusters[mag] = set()
        for contig in contigs:
            name = mag + ":" + contig
            mobile_clusters[mag].add(name)
    bins = dict()
    colors = dict()
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


def load_cluster_contigs(vmags, filename):
    print(f"loading from {filename})")
    clusters = dict()
    for cluster in vmags:
        clusters[cluster] = set()
    for line in open(filename):
        cut = line.strip().split("\t")
        contig = cut[0]
        cluster = cut[1]
        if cluster in vmags:
            clusters[cluster].add(contig)
    return clusters


def select_relevant(hosts, vmags, host_link_data):
    selected_hosts = set()
    selected_vmags = dict()

    # remove irrelevant vmags
    for vmag, contigs in vmags.items():
        contigs_with_hosts = 0
        max_host_ct = 0
        for contig in contigs:
            host_ct = 0
            if contig in host_link_data:
                for host in host_link_data[contig]:
                    if host in hosts:
                        host_ct += 1
            if host_ct > 0:
                contigs_with_hosts += 1
            if host_ct > max_host_ct:
                max_host_ct = host_ct
        if contigs_with_hosts < len(contigs):
            continue
        if contigs_with_hosts < 3:
            continue
        if len(contigs) < 3:
            continue
        if max_host_ct > 5:
            continue
        selected_vmags[vmag] = contigs

    # remove irrelevant hosts
    for contigs in selected_vmags.values():
        for contig in contigs:
            if contig in host_link_data:
                for host in host_link_data[contig]:
                    if host in hosts:
                        selected_hosts.add(host)

    return selected_hosts, selected_vmags


def main():
    host_completions, host_contaminations = load_checkm("proxiwrap_clusters.checkm.tsv")
    hosts = set(host_completions.keys())
    hosts = select_hosts(hosts, host_completions, host_contaminations, min_comp=0, max_cont=10)

    vmag_completions, vmag_contaminations = load_hifi_validation("viral_mags.hifi.tsv")
    vmags = set(vmag_completions.keys())
    vmags = select_vmags(vmags, vmag_completions, vmag_contaminations, min_comp=-10, max_cont=1000)
    vmags_contigs = load_cluster_contigs(vmags, "viral_mags.contigs.tsv")

    host_link_data = load_host_links("filtered_master_table.tsv")
    hosts, vmags_contigs = select_relevant(hosts, vmags_contigs, host_link_data)

    host_matrix = make_matrix(hosts, vmags_contigs, host_link_data)
    make_heatmap(host_matrix, mobile_clusters=vmags_contigs)

    n = 0
    for contigs in vmags_contigs.values():
        n += len(contigs)
    print(f"Kept {len(vmags_contigs)} vMAGs containing {n} contigs, linked to {len(hosts)} host MAGs")

    print(f"Exporting selected MAGs and vMAGs")
    with open("chosen_hosts.list", "w+") as output:
        for host in hosts:
            print(host, file=output)
    with open("chosen_vmags.list", "w+") as output:
        for vmag in vmags_contigs:
            print(vmag, file=output)


main()
