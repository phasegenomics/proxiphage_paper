#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Gherman Uritskiy
Computational Biologist, Phase Genomics

This script checks validity of prophage host assignment using a long-read assembly
Copyright 2021, Phase Genomics Inc. All rights reserved.

The contents of this file are proprietary and private and are not intended for
distribution or use by any person or entity except Phase Genomics. You may not
use, modify, or distribute it in any fashion. You may not copy this file. You
may not describe the contents of this file to any other party.
"""

import proxiphage.proxiphage_utilities as utils
import pgproxiwrap.pgproxiwrap_utilities as wutils
import proxiphage.plasmid_utilities as putils
import proxiphage.proxiphage_plotting as plotting
from matplotlib import pyplot as plt
import argparse

MIN_PIDENT = 95
MIN_ALIGNMENT_LEN = 100
MIN_MAG_COMPLETION = 50
MAX_MAG_CONTAMINATION = 10
MIN_VIRUS_LENGTH = 5000
MIN_BACT_LENGTH_NEXT_TO_PROPHAGE = 200000
MIN_OVERLAP_FOR_VALIDATION = 10000


def select_good_mags(mags, mag_stats, min_comp, max_cont):
    """
    Keep only MAGs with a significant completion
    Args:
        mags (dict[str:dict[str:str]]): all bins
        mag_stats (dict[str:[dict[float]]]): checkm metrics for all bins
        min_comp (float): min completion
        max_cont (float): max contamination
    Returns:
        good_mags (dict[str:dict[str:str]]): only good bins

    """
    good_mags = dict()
    for mag, contigs in mags.items():
        comp = mag_stats["Completeness"][mag]
        cont = mag_stats["Contamination"][mag]
        if comp > min_comp and cont < max_cont:
            good_mags[mag] = contigs
    print(f"selected {len(good_mags)} good MAGs from {len(mags)} MAGs")
    return good_mags


def select_good_viruses(viruses, min_len):
    """
    Keep only long viral contigs
    Args:
        viruses (dict[str:str]): viral contigs
        min_len (float): minimum viral lengths
    Returns:
        good_viruses (dict[str:str]): only long viral contigs

    """
    good_viruses = dict()
    for virus, seq in viruses.items():
        if len(seq) > min_len:
            good_viruses[virus] = seq
    print(f"selected {len(good_viruses)} good viruses from {len(viruses)} viruses")
    return good_viruses


def get_good_pairs(host_links, viruses, mags):
    """
    Keep only links between good viruses and good bins
    Args:
        host_links (dict[str:dict[str:dict]]): host linkage data between viruses and bins
        viruses (dict[str:str]): only long viral contigs
        mags (dict[str:dict[str:str]]): only good bins
    Returns:
        links (dict[str:dict[str:dict]]): host linkage data between good viruses and good bins

    """
    links = list()
    total_links = 0
    for virus, hosts in host_links.items():
        total_links += len(hosts)
        if virus in viruses:
            for mag, subdata in hosts.items():
                if mag in mags:
                    links.append((virus, mag))
    print(f"selected {len(links)} verifiable virus-host links from {total_links} links")
    return links


def load_blast(filename, min_ident, min_len):
    """
    Load alignment events from BLAST alignment file
    Args:
        filename (str): path to blast output (in outformat6 format)
        min_ident (float): minimum percent identity for the alignment to be considered [default: 97]
        min_len (int): minimum length of alignment to be considered [default: 300]
    Returns:
        hits (dict[str:dict[str:tuple]]): alignment events identified
    """
    print(f"Loading blast results from {filename}")
    hits = dict()
    with open(filename) as file_handle:
        for line in file_handle:
            info = putils.BlastOutfmt6Line(line)
            if info.pident >= min_ident and info.length >= min_len:
                if info.qseqid not in hits:
                    hits[info.qseqid] = dict()
                if info.sseqid not in hits[info.qseqid]:
                    hits[info.qseqid][info.sseqid] = list()
                hits[info.qseqid][info.sseqid].append(info)
    return hits


def load_lengths(filename):
    """
    Load contig lengths either from a tsv of fasta file
    Args:
        filename (str): path to file
    Returns:
        lengths (dict[str:int]): length of each contig

    """
    lengths = dict()
    if filename.endswith("fasta") or filename.endswith("fa"):
        long_seqs = utils.load_fasta(filename)
        for name, seq in long_seqs.items():
            lengths[name] = len(seq)
    else:
        with open(filename) as handle:
            for line in handle:
                cut = line.strip().split("\t")
                lengths[cut[0]] = int(cut[1])
    return lengths


def calculate_lengths(contigs):
    """
    Calculate contig lengths
    Args:
        contigs (dict[str:str]): contig sequences
    Returns:
        lengths (dict[str:int]): contig lengths

    """
    lengths = dict()
    for contig, seq in contigs.items():
        lengths[contig] = len(seq)
    return lengths


def get_best_virus_references(viruses, blast_hits):
    """
    Find the best long-read reference for each viral contig
    Args:
        viruses (dict[str:str]): viral contigs
        blast_hits (dict[str:dict[str:tuple]]): alignment events identified
    Returns:
        virus_best_references (dict[str:str]): best long-read reference for each viral contig

    """
    virus_best_references = dict()
    for virus in viruses:
        if virus in blast_hits:
            best_long_hit = None
            best_long_hit_aligned = 0
            for long_contig, alignments in blast_hits[virus].items():
                intervals = list()
                for info in alignments:
                    intervals = putils.merge_alignments(intervals, info.qstart, info.qend)
                aligned_total = putils.get_total_aligned_length(intervals)
                if aligned_total > best_long_hit_aligned:
                    best_long_hit_aligned = aligned_total
                    best_long_hit = long_contig
            virus_best_references[virus] = (best_long_hit, best_long_hit_aligned)
    print(f"found long-read reference for {len(virus_best_references)} viruses")
    return virus_best_references


def find_prophages(viruses, virus_best_references, long_read_contig_lengths, short_read_contig_lengths, min_extra_len):
    """
    Find likely prophages
    Args:
        viruses (dict[str:str]): all viral contigs
        virus_best_references (dict[str:str]): best long-read reference contig for each virus
        long_read_contig_lengths (dict[str:int]): long-read contig lengths
        short_read_contig_lengths (dict[str:int]): short-read contig lengths
        min_extra_len (int): minimum extra unaliged sequence for prophage identification
    Returns:
        prophages (dict[str:str]): only viral prophages
    """
    prophages = dict()
    for virus, seq in viruses.items():
        if virus in virus_best_references:
            best_long_hit, best_long_hit_aligned = virus_best_references[virus]
            percent_virus_aligned = 100 * best_long_hit_aligned / short_read_contig_lengths[virus]
            if percent_virus_aligned > 50:
                unaligned_on_long_read_contig = long_read_contig_lengths[best_long_hit] - best_long_hit_aligned
                if unaligned_on_long_read_contig > min_extra_len:
                    prophages[virus] = seq
    print(f"found {len(prophages)} prophages from a pool of {len(viruses)} viruses")
    return prophages


def count_hosts(host_virus_pairs):
    """
    Count the number of hosts found for each virus
    Args:
        host_virus_pairs (list[tuple[str]]): virus and host pairs
    Returns:
        host_counts (dict[str:int]): number of hosts found for each virus

    """
    host_counts = dict()
    for virus, mag in host_virus_pairs:
        if virus not in host_counts:
            host_counts[virus] = 1
        else:
            host_counts[virus] += 1
    return host_counts


def validate_host_links(host_virus_pairs, mags, virus_best_references, blast_hits, min_alignment_len,
                        remove_promiscuous_prophages=True, draw_histo=False):
    """
    Check if prophage hosts are supported by long-reads and generate a summary report
    Args:
        host_virus_pairs (list[tuple[str]]): prophage virus and host pairs
        mags (dict[str:dict[str:str]]): good genomic bins
        virus_best_references (dict[str:str]): best long-read contig reference for each virus
        blast_hits (dict[str:dict[str:tuple]]): alignment events identified
        min_alignment_len (int): minimum alignment length between long-read reference and bin to be a hit
        remove_promiscuous_prophages (bool): only keep prophages with 1 host
        draw_histo (bool): draw a histogram of number of hosts in confirmed and unconfirmed prophages
    Returns: None

    """
    validated = 0
    unvalidated = 0
    host_counts = count_hosts(host_virus_pairs)
    host_counts_validated = list()
    host_counts_unvalidated = list()
    for virus, mag in host_virus_pairs:
        #if remove_promiscuous_prophages and host_counts[virus] > 1:
        #    continue
        best_long_hit, _ = virus_best_references[virus]
        aligned_between_long_read_contig_and_mag = 0
        for contig in mags[mag]:
            if contig != virus and contig in blast_hits:
                if best_long_hit in blast_hits[contig]:
                    intervals = list()
                    for info in blast_hits[contig][best_long_hit]:
                        putils.merge_alignments(intervals, info.qstart, info.qend)
                    aligned_between_long_read_contig_and_mag += putils.get_total_aligned_length(intervals)
        if aligned_between_long_read_contig_and_mag > min_alignment_len:
            validated += 1
            host_counts_validated.append(host_counts[virus])
        else:
            unvalidated += 1
            host_counts_unvalidated.append(host_counts[virus])
    print(f"{validated} out of {unvalidated + validated} prophage host links were validated with the "
          f"long-read assembly ({int(100 * validated / (validated + unvalidated))}%)")
    if draw_histo:
        draw_histograms(host_counts_validated, host_counts_unvalidated)


def draw_histograms(data1, data2):
    """
    Make simple overlapping histogram
    Args:
        data1 (list[int]): first value list
        data2 (list[int]): second value list
    Returns: None

    """
    plotting.initialize_plot_style()
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(data1, bins=7, color="k", alpha=0.5)
    ax.hist(data2, bins=7, color="r", alpha=0.5)
    ax.set_xlabel("number of hosts")
    ax.set_ylabel("number of phages")
    ax.set_title("Number of hosts in the validated and unvalidated virus-host links categories")
    ax.grid(axis="both", ls="-", alpha=0.1, c="k", which="major")
    ax.grid(axis="both", ls="--", alpha=0.1, c="k", which="minor")
    plt.savefig("histogram.png", bbox_inches='tight', dpi=300)
    plt.close()


def run(alignment_file, bins_dir, bins_stats_file, viral_file, long_lengths_file, host_file):
    """
    Run prophage host validation analysis
    Args:
        alignment_file (str): path to blast alignment file
        bins_dir (str): path to MAG directory
        bins_stats_file (str): path to MAG CheckM stats
        viral_file (str): path to viral fasta file
        long_lengths_file (str): path to long-read contig length information
        host_file (str): path to filtered host master table
    Returns: None

    """
    blast_hits = load_blast(alignment_file, min_ident=MIN_PIDENT, min_len=MIN_ALIGNMENT_LEN)
    mags = utils.load_bin_set(bins_dir)
    mag_stats = wutils.load_checkm_results(bins_stats_file)
    mags = select_good_mags(mags, mag_stats, MIN_MAG_COMPLETION, MAX_MAG_CONTAMINATION)
    viruses = utils.load_fasta(viral_file)
    viruses = select_good_viruses(viruses, MIN_VIRUS_LENGTH)
    virus_best_references = get_best_virus_references(viruses, blast_hits)
    long_read_contig_lengths = load_lengths(long_lengths_file)
    short_read_contig_lengths = calculate_lengths(viruses)
    prophages = find_prophages(viruses, virus_best_references, long_read_contig_lengths,
                               short_read_contig_lengths, MIN_BACT_LENGTH_NEXT_TO_PROPHAGE)
    host_links = utils.load_host_connection_data(host_file)
    host_virus_pairs = get_good_pairs(host_links, prophages, mags)
    validate_host_links(host_virus_pairs, mags, virus_best_references, blast_hits, MIN_OVERLAP_FOR_VALIDATION)


def parse_args(desc, cli_args=None, parents=None):
    """
    Simple command line argument parser
    Args:
        desc (str): short program description, often __file__
        cli_args (list): all CLI arguments passed (such as sys.argv[1:])
        parents ([argparse.ArgumentParser]): List of parent parsers to add to the parser
    Returns:
        vars(args) (dict): Dict of command line args with key=argument destination, val=argument value

    """
    print("Parsing arguments")
    if parents is None:
        parents = list()
    parser = argparse.ArgumentParser(description=desc, parents=parents)

    parser.add_argument("--alignment", required=True, type=str, dest="alignment",
                        help="Long- vs short-read contig blast alignment")
    parser.add_argument("--bins", required=True, type=str, dest="bins",
                        help="ProxiWrap bins")
    parser.add_argument("--stats", required=True, type=str, dest="stats",
                        help="ProxiWrap bin stats")
    parser.add_argument("--viruses", required=True, type=str, dest="viruses",
                        help="Viral contig fasta")
    parser.add_argument("--lengths", required=True, type=str, dest="lengths",
                        help="Long-read contig lengths tsv or fasta")
    parser.add_argument("--hosts", required=True, type=str, dest="hosts",
                        help="Filtered master hosts table")
    args = parser.parse_args(args=cli_args)
    return vars(args)


def main():
    """
    Run script
    Returns: None

    """
    c_args = parse_args(__file__)
    run(c_args["alignment"],
        c_args["bins"],
        c_args["stats"],
        c_args["viruses"],
        c_args["lengths"],
        c_args["hosts"])


if __name__ == '__main__':
    """ Launch the script
    """
    main()
