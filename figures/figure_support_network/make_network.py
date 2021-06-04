#!/usr/bin/env python3
import sys
import random


def load_nodes(filename):
    print(f"loading from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                header = cut
            else:
                contig = cut[0]
                vmag = cut[4]
                if vmag not in data:
                    data[vmag] = dict()
                data[vmag][contig] = line.strip()
    return data


def load_edges(filename):
    print(f"loading from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                header = cut
            else:
                contig = cut[0]
                data[contig] = line.strip()
    return data


def choose_random_vmags(vmags, vmag_count):
    print(f"down sampling from {len(vmags)} to {vmag_count} vMAGs")
    chosen_vmags = dict()
    while len(chosen_vmags) < vmag_count and len(vmags) > 0:
        n = random.randint(0, len(vmags) - 1)
        vmag = list(vmags.keys())[n]
        if len(vmags[vmag]) < 3:
            continue
        chosen_vmags[vmag] = vmags[vmag]
        del vmags[vmag]
    print(f"{len(vmags)} vmags left unselected")
    return chosen_vmags


def filter_edges(all_edges, nodes):
    print(f"filtering {len(all_edges)} edges relevant to {len(nodes)} nodes")
    relevant_contigs = set()
    for contigs in nodes.values():
        for contig in contigs:
            relevant_contigs.add(contig)
    edges = dict()
    for contig1, line in all_edges.items():
        contig2 = line.split("\t")[1]
        if contig1 in relevant_contigs or contig2 in relevant_contigs:
            edges[contig1] = line
    print(f"kept {len(edges)} relevant edges")
    return edges


def export_nodes(nodes, filename):
    print(f"exporting to {filename}")
    with open(filename, "w+") as output:
        print("contig", "type", "length", "cluster", "cluster_n", sep="\t", file=output)
        for line in nodes.values():
            print(line, file=output)


def export_edges(edges, filename):
    print(f"exporting to {filename}")
    with open(filename, "w+") as output:
        print("viral_contig", "ont_contig", "aligned_length", "percent_aligned", sep="\t", file=output)
        for alignments in edges.values():
            for line in alignments.values():
                print(line, file=output)


def all_contigs_aligned(vmags, alignments):
    print(f"scanning {len(vmags)} vmags to make sure all their contigs aligned to something")
    fixed_vmags = dict()
    for vmag, contigs in vmags.items():
        aligned_ct = 0
        for contig in contigs:
            if contig in alignments:
                if max(alignments[contig].values()) >= 1000:
                    aligned_ct += 1
        if aligned_ct == len(contigs):
            fixed_vmags[vmag] = contigs
    print(f"kept {len(fixed_vmags)} vmags")
    return fixed_vmags


def load_contig_lengths(filename):
    print(f"loading from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            contig = cut[0]
            l = int(cut[1])
            data[contig] = l
    return data


def load_alignments(filename):
    print(f"loading from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            contig = cut[0]
            reference = cut[1]
            pident = float(cut[2])
            l = int(cut[3])
            if pident > 95 and l > 100:
                if contig not in data:
                    data[contig] = dict()
                if reference not in data[contig]:
                    data[contig][reference] = 0
                data[contig][reference] += l
    filtered_data = dict()
    for contig, alignments in data.items():
        for reference, l in alignments.items():
            if l > 1000:
                if contig not in filtered_data:
                    filtered_data[contig] = dict()
                filtered_data[contig][reference] = l
    return filtered_data


def load_vmags(filename):
    print(f"loading from {filename}")
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            contig = cut[0]
            cluster = cut[1]
            if cluster not in data:
                data[cluster] = set()
            data[cluster].add(contig)
    return data


def make_nodes(vmags, contig_lengths):
    print(f"converting contigs from {len(vmags)} vmags into nodes")
    data = dict()
    for vmag, contigs in vmags.items():
        for contig in contigs:
            l = contig_lengths[contig]
            info = [contig, "short_contig", str(l), vmag, vmag.split("_")[-1]]
            line = "\t".join(info)
            data[contig] = line
    return data


def make_edges(nodes, alignements, contig_lengths):
    print(f"making edges")
    data = dict()
    for contig in nodes:
        if contig in alignements:
            for hifi, l in alignements[contig].items():
                percent_aligned = 100 * l / contig_lengths[contig]
                info = [contig, hifi, str(l), str(percent_aligned)]
                line = "\t".join(info)
                if contig not in data:
                    data[contig] = dict()
                data[contig][hifi] = line
    return data


def add_hifi_to_nodes(nodes, edges, reference_lengths):
    references = set()
    for alignments in edges.values():
        for reference in alignments:
            references.add(reference)
    print(f"adding {len(references)} reference sequences to nodes")
    for reference in references:
        l = reference_lengths[reference]
        info = [reference, "reference_phage", str(l), "", "0"]
        line = "\t".join(info)
        nodes[reference] = line
    return nodes


def main():
    #load data
    contig_lengths = load_contig_lengths("sh_assembly_1kb.lengths.tsv")
    reference_lengths = load_contig_lengths("long_read_excised_phages.lengths.tsv")
    alignements = load_alignments("long_read_excised_phages.blast.tsv")
    vmags = load_vmags("viral_mags.contigs.tsv")

    # choose what to plot
    vmags = all_contigs_aligned(vmags, alignements)
    vmags = choose_random_vmags(vmags, 50)

    # make network
    nodes = make_nodes(vmags, contig_lengths)
    edges = make_edges(nodes, alignements, contig_lengths)
    nodes = add_hifi_to_nodes(nodes, edges, reference_lengths)
    export_nodes(nodes, "nodes.tsv")
    export_edges(edges, "edges.tsv")
main()
