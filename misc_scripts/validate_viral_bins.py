#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Gherman Uritskiy
Computational Biologist, Phase Genomics

proxiphage/misc_scripts/validate_against_long_assembly.py

This script loads in viral MAGs from ProxiPhage and compares their accuracy by comparing them against a long-read
assembly (a pre-made blast alignment is required).

Copyright 2020, Phase Genomics Inc. All rights reserved.

The contents of this file are proprietary and private and are not intended for
distribution or use by any person or entity except Phase Genomics. You may not
use, modify, or distribute it in any fashion. You may not copy this file. You
may not describe the contents of this file to any other party.
"""

import os
import proxiphage.proxiphage_utilities as utils
import proxiphage.plasmid_utilities as putils
import proxiphage.proxiphage_plotting as plotting
import numpy as np
import argparse


MIN_COMPLETION = 50
MAX_CONTAMINATION = 10
MIN_PIDENT = 95
MIN_ALI_LEN = 100
MIN_PCOVER = 0
END_DISTANCE = 10
MIN_TOTAL_ALIGNED = 5000
MIN_TOTAL_PALI = 80

CONFLICTING = "conflicting"
UNVALIDATED = "unvalidated"
PARTIALLY_VALIDATED = "partially_validated"
PERFECTLY_VALIDATED = "perfectly_validated"
TOTAL_CONTIGS = "total_contigs"
TOTAL_LENGTH = "total_length"
TOTAL_MAGS = "total_mags"


def load_blast(filename, short_lengths, long_lengths, min_ident=MIN_PIDENT, min_len=MIN_ALI_LEN):
    """
    Load alignment events from BLAST alignment file
    Args:
        filename (str): path to blast output (in outformat6 format)
        short_lengths (dict): lengths of contig_pool in short-read assembly
        long_lengths (dict): lengths of contig_pool in long-read assembly
        min_ident (float): minimum percent identity for the alignment to be considered [default: 97]
        min_len (int): minimum length of alignment to be considered [default: 300]
    Returns:
        hits (dict[str:dict[str:BlastOutfmt6Line]]): alignment events identified

    """
    print(f"Loading blast results from {filename}")
    hits = dict()
    with open(filename) as file_handle:
        for line in file_handle:
            if len(line) == 0:
                continue
            info = putils.BlastOutfmt6Line(line)
            if info.qseqid in short_lengths:
                info.query_len = short_lengths[info.qseqid]
            if info.sseqid in long_lengths:
                info.subject_len = long_lengths[info.sseqid]
            if info.pident >= min_ident and info.length >= min_len:
                if info.qseqid not in hits:
                    hits[info.qseqid] = dict()
                if info.sseqid not in hits[info.qseqid]:
                    hits[info.qseqid][info.sseqid] = list()
                hits[info.qseqid][info.sseqid].append(info)
    return hits


def print_validation_report(mags, hits, outdir, stats=None, report_bad=True, hosts=None):
    """
    Print out support information for the viral MAGs based on the long-read assembly
    Args:
        mags (dict[dict]): the viral clusters
        hits (dict): alignment events between the long-read and short-read assemblies
        stats (dict): CheckV quality of viral clusters and sequences
        report_bad (bool): should bad MAGs be reported on in detail? [default: True]
        hosts (dict[str:dict[str:dict[str]]]): host information for each virus
        outdir (str): path to desired output to save info to
    Returns: None

    """
    print("Testing cluster support")
    mag_support_handle = open(os.path.join(outdir, "mag_support.tsv"), "w+")
    summary_handle = open(os.path.join(outdir, "support_summary.txt"), "w+")
    n_support_counts = {
        CONFLICTING: 0,
        UNVALIDATED: 0,
        PARTIALLY_VALIDATED: 0,
        PERFECTLY_VALIDATED: 0,
        TOTAL_MAGS: len(mags),
        TOTAL_CONTIGS: 0,
        TOTAL_LENGTH: 0
    }
    print("cluster", "n_contigs", "n_contigs_aligned_to_best_reference", "n_contigs_aligned", "verdict",
          "cluster_length", "reference_length", "completion", "contamination", sep="\t", file=mag_support_handle)
    for mag, contigs in list(mags.items()):
        mag_length = 0
        for seq in contigs.values():
            n_support_counts[TOTAL_CONTIGS] += 1
            mag_length += len(seq)
            n_support_counts[TOTAL_LENGTH] += mag_length
        ont_hits, aligned_contigs, max_support = tally_contig_hits(contigs, hits)
        completion, contamination, _ = calculate_completion_contamination(contigs, hits)
        hit_aligned_intervals, short_lengths, long_lengths, mag_length = summarize_alignments_of_mag(contigs, hits)
        best_long_read_contig, _ = get_best_reference_contig(hit_aligned_intervals, short_lengths)
        verdict = type_of_support(contigs, ont_hits, aligned_contigs, max_support)
        if best_long_read_contig is None:
            reference_length = 0
        else:
            reference_length = long_lengths[best_long_read_contig]
        n_support_counts[verdict] += 1
        print(mag, len(contigs), max_support, aligned_contigs, verdict, 
              mag_length, reference_length, str(completion)[:4], str(contamination)[:4],
              sep="\t", file=mag_support_handle)
        if verdict == CONFLICTING and report_bad:
            print_detailed_mistake(outdir, hits, mag, contigs, stats, hosts)
    print_summary(n_support_counts, out_handle=summary_handle)
    print_summary(n_support_counts)
    mag_support_handle.close()
    summary_handle.close()


def plot_completion_contamination(mags, hits, outdir, long_lengths):
    """
    Print out support information for the viral MAGs based on the long-read assembly
    Args:
        mags (dict[dict]): the viral clusters
        hits (dict): alignment events between the long-read and short-read assemblies
        outdir (str): path to desired output to save info to
        long_lengths (dict[str:int]): lengths of references
    Returns: None

    """
    print("Plotting completion and contamination metrics")
    individual_contig_mags = generate_single_contig_mags(mags)

    sample_names = list()
    completion_data = dict()
    contamination_data = dict()

    mag_names, reference_names, completions, contaminations = make_completion_contamination_lists(
        individual_contig_mags, hits)
    sample_name = "viral_contigs"
    export_completion_contamination_tsv(individual_contig_mags, mag_names, reference_names, completions, contaminations,
                                        os.path.join(outdir, f"{sample_name}.tsv"), long_lengths)
    sample_names.append(sample_name)
    completion_data[sample_name] = completions
    contamination_data[sample_name] = contaminations

    mag_names, reference_names, completions, contaminations = make_completion_contamination_lists(mags, hits)
    sample_name = "viral_clusters"
    export_completion_contamination_tsv(mags, mag_names, reference_names, completions, contaminations,
                                        os.path.join(outdir, f"{sample_name}.tsv"), long_lengths)
    sample_names.append(sample_name)
    completion_data[sample_name] = completions
    contamination_data[sample_name] = contaminations

    plotting.make_report(sample_names, completion_data, contamination_data,
                         os.path.join(outdir, "quality_report.png"), min_yval=20)
    plotting.make_barplot(sample_names, completion_data, contamination_data,
                          os.path.join(outdir, "quality_barplot.png"))


def export_completion_contamination_tsv(mags, mag_names, reference_names, completions, contaminations, file_name,
                                        reference_lengths):
    """
    Export tsv report with MAG completion and contamination
    Args:
        mags (dict[str:dict[str:str]]): the viral clusters
        reference_names (list[str]): ordered liks of reference long-read contigs
        mag_names (list[str]: ordered MAG names
        completions (list[int]): ordered completion values
        contaminations (list[int]): ordered contamination values
        file_name (str): path to output file
        reference_lengths (dict[str:int]): lengths of reference sequences
    Returns: None

    """
    print(f"Exporting completion/contamination tsv report {file_name}")
    good_mags = 0
    contaminated_mags = 0
    with open(file_name, "w+") as output_handle:
        print("name", "reference", "contig_count", "length", "completion", "contamination",
              sep="\t", file=output_handle)
        for i, mag_name in enumerate(mag_names):
            reference = reference_names[i]
            if reference is not None:
                reference_length = reference_lengths[reference]
            else:
                reference_length = 0
            completion = completions[i]
            contamination = contaminations[i]
            if contamination > MAX_CONTAMINATION:
                contaminated_mags += 1
            if completion > MIN_COMPLETION and contamination < MAX_CONTAMINATION:
                good_mags += 1
            n_contigs = len(mags[mag_name])
            l_contigs = 0
            for seq in mags[mag_name].values():
                l_contigs += len(seq)
            print(mag_name, f"{reference}: {reference_length}bp", n_contigs, l_contigs, completion, contamination,
                  sep="\t", file=output_handle)
    print(f"Found {good_mags} good viral MAGs and {contaminated_mags} over-contaminated MAGs")


def generate_single_contig_mags(mags):
    """
    Generate single-contig MAGs
    Args:
        mags (dict[dict[str:str]]): the viral clusters
    Returns:
        single_contig_mags (dict[dict[str:str]]): the single-contig clusters

    """
    print("Generate single-contig 'MAGs'")
    single_contig_mags = dict()
    for contigs in mags.values():
        for contig, seq in contigs.items():
            single_contig_mags[contig] = dict()
            single_contig_mags[contig][contig] = seq
    return single_contig_mags


def make_completion_contamination_lists(mags, hits):
    """
    Generate sorter lists of completion and contamination values for MAGs
    Args:
        mags (dict[dict[str:str]]): the viral clusters
        hits (dict): alignment events between the long-read and short-read assemblies
    Returns:
        mag_names (list[str]): names of mag in order with completion values
        reference_names (list[str]): names of mag reference long-read contigs in order with completion values
        completions (list[float]): sorted list of completion values
        contaminations (list[float]): list of contamination values sorted based on the completion values

    """
    mag_names = list()
    reference_names = list()
    completions = list()
    contaminations = list()
    for mag, contigs in list(mags.items()):
        mag_names.append(mag)
        completion, contamination, best_reference = calculate_completion_contamination(contigs, hits)
        if completion is None:
            completion = 0
            contamination = 0
        completions.append(completion)
        contaminations.append(contamination)
        reference_names.append(best_reference)
    mag_names = [x for _, x in sorted(zip(completions, mags), reverse=True)]
    reference_names = [x for _, x in sorted(zip(completions, reference_names), reverse=True)]
    contaminations = [x for _, x in sorted(zip(completions, contaminations), reverse=True)]
    completions.sort(reverse=True)
    return mag_names, reference_names, completions, contaminations


def print_summary(n_support_counts, out_handle=None):
    """
    Print summary of support counts
    Args:
        n_support_counts (dict[str:int]): counts of different support categiries
        out_handle (handle): file to write to
    Returns: None

    """
    print(f"\nNumber of contigs in vMAGs: {n_support_counts['total_contigs']}", file=out_handle)
    print(f"Length of contigs in vMAGs: {n_support_counts['total_length']}", file=out_handle)
    print(f"Number of total vMAGs: {n_support_counts['total_mags']}", file=out_handle)
    print(f"Number of conflicting vMAGs: {n_support_counts['conflicting']} "
          f"({int(100 * n_support_counts['conflicting'] / n_support_counts['total_mags'])}%)", file=out_handle)
    print(f"Number of unvalidated vMAGs: {n_support_counts['unvalidated']} "
          f"({int(100 * n_support_counts['unvalidated'] / n_support_counts['total_mags'])}%)", file=out_handle)
    print(f"Number of partially validated vMAGs: {n_support_counts['partially_validated']} "
          f"({int(100 * n_support_counts['partially_validated'] / n_support_counts['total_mags'])}%)", file=out_handle)
    print(f"Number of perfectly validated vMAGs: {n_support_counts['perfectly_validated']} "
          f"({int(100 * n_support_counts['perfectly_validated'] / n_support_counts['total_mags'])}%)", file=out_handle)


def print_detailed_mistake(outdir, hits, mag, contigs, stats, hosts):
    """
    Print detailed mistakes report of conflicting MAGs
    Args:
        mag (str): MAG name
        contigs (dict[str:str]): contig sequences in MAG
        stats (dict): CheckV quality of viral clusters and sequences
        hosts (dict[str:dict[str:dict[str]]]): host information for each virus
        outdir (str): path to desired output to save info to
        hits (dict[str:dict[str:BlastOutfmt6Line]]): blast alignments
    Returns: None

    """
    with open(os.path.join(outdir, "detailed_mistakes.txt"), "w+") as out_handle:
        print(mag, file=out_handle)
        info = [mag]
        if mag in stats:
            info += stats[mag]
        else:
            info.append("N/A")
        print(("\t".join(info)), file=out_handle)

        for contig, seq in list(contigs.items()):
            if contig in hits:
                contig_hits = hits[contig]
            else:
                contig_hits = ["N/A"]
            if contig in hosts:
                contig_hosts = list(hosts[contig])
            else:
                contig_hosts = ["N/A"]
            print(f"{contig}: {len(seq)}\t{';'.join(contig_hits)}\t{';'.join(contig_hosts)}", file=out_handle)

        for contig in contigs:
            contig_id = "_".join(contig.split("_")[:2])
            if contig_id in stats:
                info = [contig]
                info += stats[contig_id]
                print(("\t".join(info)), file=out_handle)
        print("", file=out_handle)


def tally_contig_hits(contigs, hits):
    """
    Summarize the blast hits of contigs in a cluster
    Args:
        contigs (dict[str:str]): contigs in a cluster
        hits (dict[str:dict[str:tuple]]): blast hits of short-read contigs to long-read contigs
    Returns:
        ont_hits (set): long-read contigs hit by the shot-read contig
        aligned_contigs (int): number of contigs that were aligned to the long-read contigs
        max_support (int): maximum number of short-read contigs aligned to a single long-read contig

    """
    ont_hits = set()
    supports = dict()
    aligned_contigs = 0
    for contig in contigs:
        if contig in hits:
            aligned_contigs += 1
            ont_hits = ont_hits.union(set(hits[contig].keys()))
            for subcontig in hits[contig]:
                if subcontig not in supports:
                    supports[subcontig] = 1
                else:
                    supports[subcontig] += 1
    if len(supports) >= 1:
        max_support = max(supports.values())
    else:
        max_support = 0
    return ont_hits, aligned_contigs, max_support


def calculate_completion_contamination(contigs, blast_hits):
    """
    Calculate the completion and contamination of a MAG based on long read assembly
    Args:
        contigs (dict[str:dict[str:str]]): contigs in MAG
        blast_hits (dict[str:dict[str:tuple]]): blast alignment data
    Returns:
        completion (float): estimated percent completion of MAG
        contamination (float): estimated percent contamination of MAG
        best_long_read_contig (str): name of best reference long-read contig

    """
    hit_aligned_intervals, short_lengths, long_lengths, mag_length = summarize_alignments_of_mag(contigs, blast_hits)
    best_long_read_contig, best_long_read_contig_aligned = get_best_reference_contig(
        hit_aligned_intervals, short_lengths)
    unalignable_length = get_unalignable_length(contigs, blast_hits, best_long_read_contig)
    if unalignable_length > mag_length - best_long_read_contig_aligned:
        unalignable_length = mag_length - best_long_read_contig_aligned

    if best_long_read_contig is None:
        completion = None
        contamination = None
    else:
        best_long_read_contig_coverage = calculate_coverage_of_best_contig(best_long_read_contig, contigs, blast_hits)
        completion = 100 * best_long_read_contig_coverage / long_lengths[best_long_read_contig]
        contamination = 100 * (mag_length - best_long_read_contig_aligned - unalignable_length) / \
                        (mag_length - unalignable_length)
    return completion, contamination, best_long_read_contig


def get_unalignable_length(contigs, blast_hits, intended_reference_contig):
    """
    Calculate length of sequence in MAG that did not align to anything
    Args:
        contigs (dict[str:dict[str:str]]): contigs in MAG
        blast_hits (dict[str:dict[str:tuple]]): blast alignment data
        intended_reference_contig (str): intended long-read contig reference
    Returns:
        unalignable_length (int): length of sequence in MAG that did not align to anything

    """
    unalignable_length = 0
    if intended_reference_contig is not None:
        for contig, seq in contigs.items():
            contig_length = len(seq)
            hit_aligned_intervals, short_lengths, _, mag_length = summarize_alignments_of_mag({contig: seq}, blast_hits)
            _, best_long_read_contig_aligned = get_best_reference_contig(hit_aligned_intervals, short_lengths)
            edge_alignment = False
            aligned_to_intended_reference = 0
            if contig in blast_hits:
                if intended_reference_contig in blast_hits[contig]:
                    alignments = blast_hits[contig][intended_reference_contig]
                    intervals = list()
                    for info in alignments:
                        intervals = putils.merge_alignments(intervals, info.qstart, info.qend)
                        if min(info.sstart, info.send) <= END_DISTANCE or \
                                max(info.sstart, info.send) >= info.subject_len - END_DISTANCE:
                            # this is a edge alignment and any missing sequence should not be held against bin
                            edge_alignment = True
                    aligned_to_intended_reference = putils.get_total_aligned_length(intervals)
            if edge_alignment:
                unalignable_contig_length = contig_length - aligned_to_intended_reference
            else:
                unalignable_contig_length = contig_length - best_long_read_contig_aligned
            unalignable_length += unalignable_contig_length

    return unalignable_length


def calculate_coverage_of_best_contig(best_long_read_contig, contigs, blast_hits):
    """
    Calculate the length of the best long-read contig that aligned to the MAG
    Args:
        best_long_read_contig (str): name of best long-read contig reference
        contigs (dict[str:dict[str:str]]): contigs in MAG
        blast_hits (dict[str:dict[str:tuple]]): blast alignment data
    Returns:
        best_long_read_contig_coverage

    """
    # find coverage of best long-read hit
    intervals = list()
    for short_read_contig in contigs:
        if short_read_contig in blast_hits:
            if best_long_read_contig in blast_hits[short_read_contig]:
                for info in blast_hits[short_read_contig][best_long_read_contig]:
                    intervals = putils.merge_alignments(intervals, info.sstart, info.send)
    best_long_read_contig_coverage = putils.get_total_aligned_length(intervals)
    return best_long_read_contig_coverage


def summarize_alignments_of_mag(contigs, blast_hits):
    """
    Summarize aligned intervals of a mag to all long-read contigs
    Args:
        contigs (dict[str:dict[str:str]]): contigs in MAG
        blast_hits (dict[str:dict[str:tuple]]): blast alignment data
    Returns:
        hit_aligned_intervals (dict[str:dict[str:tuple[int]]]): alignment intervals between long- and short-read contigs
        short_lengths (dict[str:int]): lengths of short-read contigs
        long_lengths (dict[str:int]): lengths of long-read contigs
        mag_length (int): total length of this mag

    """
    hit_aligned_intervals = dict()
    short_lengths = dict()
    long_lengths = dict()
    mag_length = 0
    for short_read_contig, seq in contigs.items():
        mag_length += len(seq)
        if short_read_contig in blast_hits:
            for long_read_contig, alignments in blast_hits[short_read_contig].items():
                for info in alignments:
                    short_lengths[short_read_contig] = info.query_len
                    long_lengths[long_read_contig] = info.subject_len
                    if long_read_contig not in hit_aligned_intervals:
                        hit_aligned_intervals[long_read_contig] = dict()
                    if short_read_contig not in hit_aligned_intervals[long_read_contig]:
                        hit_aligned_intervals[long_read_contig][short_read_contig] = [(info.qstart, info.qend)]
                    else:
                        hit_aligned_intervals[long_read_contig][short_read_contig] = putils.merge_alignments(
                            hit_aligned_intervals[long_read_contig][short_read_contig], info.qstart, info.qend)
    return hit_aligned_intervals, short_lengths, long_lengths, mag_length


def get_best_reference_contig(hit_aligned_intervals, short_lengths):
    """
    Find the best (the one that covers the most of the MAG) long-read contig to use as the reference for this MAG
    Args:
        hit_aligned_intervals (dict[str:dict[str:tuple[int]]]): alignment intervals between long- and short-read contigs
        short_lengths (dict[str:int]): length of short read contigs
    Returns:
        best_long_read_contig (str): name of best long-read contig
        best_long_read_contig_aligned (int): length of MAG covered by the best contig

    """
    best_long_read_contig = None
    best_long_read_contig_aligned = 0
    for long_read_contig, short_read_alignments in hit_aligned_intervals.items():
        total_aligned_length = 0
        # this is important - check that at least one contig is fully aligned to the reference, otherwise this is not it
        closest_contig_aligned_percentage = 0
        for short_contig, intervals in short_read_alignments.items():
            aligned_len = putils.get_total_aligned_length(intervals)
            total_aligned_length += aligned_len
            aligned_percentage = 100 * aligned_len / short_lengths[short_contig]
            if aligned_percentage > closest_contig_aligned_percentage and \
                    short_lengths[short_contig] > MIN_TOTAL_ALIGNED:
                closest_contig_aligned_percentage = aligned_percentage
        if total_aligned_length > best_long_read_contig_aligned and closest_contig_aligned_percentage > MIN_TOTAL_PALI:
            best_long_read_contig_aligned = total_aligned_length
            best_long_read_contig = long_read_contig
    return best_long_read_contig, best_long_read_contig_aligned


def type_of_support(contigs, ont_hits, aligned_contigs, max_support):
    """
    Classify the support for a given cluster
    Args:
        contigs (dict[str:str]): contigs in cluter
        ont_hits (set): long-read contigs hit by the shot-read contig
        aligned_contigs (int): number of contigs that were aligned to the long-read contigs
        max_support (int): maximum number of short-read contigs aligned to a single long-read contig
    Returns:
        verdict (str): type of support

    """
    if len(ont_hits) == 0:
        # the mag is not supported by any ont contig_pool - nothing to do
        verdict = UNVALIDATED
    elif aligned_contigs <= 1:
        # only one contig in the mag was aligned to the ont assembly - nothing to do
        verdict = UNVALIDATED
    elif len(contigs) > aligned_contigs == max_support:
        # not all contigs were aligned, but those that were are supported
        verdict = PARTIALLY_VALIDATED
    elif len(ont_hits) == 1:
        # all contigs in this mag are supported by a single ont contig - perfect!
        verdict = PERFECTLY_VALIDATED
    elif aligned_contigs == max_support:
        # some ambiguous hits, but all contig_pool are supported by a single contig - good
        verdict = PERFECTLY_VALIDATED
    else:
        # the contigs in this mag hit multiple ont contigs - not good
        verdict = CONFLICTING
    return verdict


def load_all_stats(bin_stats_file=None, proxiphage_dir=None):
    """
    Load CheckV quality information on viral sequencs and clusters from ProxiPhage output
    Args:
        bin_stats_file (str): path to tabular clusters stats file
        proxiphage_dir (str): path to ProxiPhage output dir
    Returns:
        seeds (set): viral seed sequences identified by ProxiPhage
        seed_stats (dict): CheckV quality information on both sequences

    """
    if proxiphage_dir is None:
        seeds = set()
        all_stats = dict()
    else:
        print("Loading all CheckV stats from proxiphage output")
        seed_stats = utils.load_checkv_results(os.path.join(proxiphage_dir, "viral_seed_sequences.tsv"))
        seeds = set(seed_stats)
        all_stats = seed_stats
        secondary_stats = utils.load_checkv_results(os.path.join(proxiphage_dir, "viral_secondary_sequences.tsv"))
        mag_stats = utils.load_checkv_results(os.path.join(proxiphage_dir, "viral_mags_final.tsv"))
        all_stats.update(secondary_stats)
        all_stats.update(mag_stats)
    if bin_stats_file is not None:
        print("Loading cluster stats from provided stats file")
        cluster_stats = utils.load_checkv_results(bin_stats_file)
        all_stats.update(cluster_stats)
    return seeds, all_stats


def curate_stats(stats):
    """
    Keep only CheckV or CheckM relevant stats
    Args:
        stats (dict[str:dict]): stats for a given set of bins of contigs
    Returns:
        curated_stats (dict[str:dict]): only expected stats for a given set of bins of contigs

    """
    checkv_fields = ['checkv_quality', 'miuvig_quality', 'completeness', 'contamination', 'gene_count', 'viral_genes',
                     'host_genes', 'termini', 'prophage', 'contig_length']
    plasmid_fields = ['cluster_name', 'contig_count', 'length', 'aligned_length', 'percent_aligned',
                      'percent_similarity', 'best_hit_id', 'reference_length', 'percent_completion',
                      'percent_unaligned']
    curated_stats = dict()
    for sequence, info in stats.items():
        curated_stats[sequence] = list()
        for field, value in info.items():
            if field in checkv_fields or field in plasmid_fields:
                curated_stats[sequence].append(value)
    return curated_stats


def calc_contig_lengths_in_mags(mags):
    """
    Calculate length of viral MAGs
    Args:
        mags (dict): viral mags
    Returns:
        lengths (dict[int]): lengths of viral MAGs

    """
    print("Calculating lengths of contigs in the bins")
    lengths = dict()
    for mag, contigs in list(mags.items()):
        for contig, seq in list(contigs.items()):
            lengths[contig] = len(seq)
    return lengths


def import_lengths(filename):
    """
    Import MAG lengths from tsv file
    Args:
        filename (str): path to tsv file with MAG lengths
    Returns:
        lengths (dict[int]): total lengths of MAGs

    """
    print(f"Loading contig lengths from {filename}")
    lengths = dict()
    with open(filename) as handle:
        for line in handle:
            cut = line.rstrip().split("\t")
            lengths[cut[0]] = int(cut[1])
    return lengths


def make_network_nodes(clusters, alignments, output_file, keep_partial_hits=True):
    """
    Generate network nodes for visualizing long-read support of clusters
    Args:
        clusters (dict): viral clusters to support
        alignments (dict): alignment events between long-read and short-read assemblies
        output_file (str): file to write nodes to
        keep_partial_hits (boo): show all clusters even if not every contig aligned to long-read contigs
    Returns:
        None

    """
    hit_long_contigs = dict()
    kept_clusters = dict()
    for cluster, sequences in clusters.items():
        all_found = True
        for contig in sequences.keys():
            if contig not in alignments:
                all_found = False
        if keep_partial_hits:
            all_found = True
        if all_found and len(sequences) > 1:
            kept_clusters[cluster] = sequences
            for contig in sequences:
                if contig in alignments:
                    for long_contig, infos in alignments[contig].items():
                        hit_long_contigs[long_contig] = infos[0].subject_len
    with open(output_file, "w+") as handle:
        n_nodes = 0
        long_contig_to_cluster_hits = dict()
        print("contig\ttype\tlength\tlog_length\tcluster\tcluster_n", file=handle)
        for cluster, sequences in kept_clusters.items():
            for contig, seq in sequences.items():
                if contig in alignments:
                    for long_contig in alignments[contig]:
                        if long_contig not in long_contig_to_cluster_hits:
                            long_contig_to_cluster_hits[long_contig] = dict()
                        if cluster not in long_contig_to_cluster_hits[long_contig]:
                            long_contig_to_cluster_hits[long_contig][cluster] = 0
                        long_contig_to_cluster_hits[long_contig][cluster] += 1
                n_nodes += 1
                length = len(seq)
                log_len = np.log10(length)
                line = "\t".join([contig, "short_contig", str(length), str(log_len), cluster, cluster.split("_")[1]])
                print(line, file=handle)
        for long_contig, length in hit_long_contigs.items():
            n_nodes += 1
            length = -10000
            log_len = 0
            affiliated_cluster = "None"
            affiliations = 0
            for cluster, links in long_contig_to_cluster_hits[long_contig].items():
                if links > affiliations:
                    affiliations = links
                    affiliated_cluster = cluster
            line = "\t".join([long_contig, "ont_contig", str(length), str(log_len), affiliated_cluster, "0"])
            print(line, file=handle)
    return kept_clusters


def make_network_edges(clusters, alignments, output_file):
    """
    Generate network edges for visualizing long-read support of clusters
    Args:
        clusters (dict): viral clusters
        alignments (dict): alignment events between long-read and short-read assemblies
        output_file (str): file to write edges to
    Returns: None

    """
    long_contig_aligned_counts = dict()
    binned_contigs = dict()
    for contigs in clusters.values():
        binned_contigs.update(contigs)
    for contig, info in alignments.items():
        if contig in binned_contigs:
            for long_contig in info:
                if long_contig not in long_contig_aligned_counts:
                    long_contig_aligned_counts[long_contig] = set()
                long_contig_aligned_counts[long_contig].add(contig)

    with open(output_file, "w+") as handle:
        n_edges = 0
        print("viral_contig\tont_contig\taligned_length\tpercent_aligned", file=handle)
        for cluster, sequences in clusters.items():
            for contig, seq in sequences.items():
                length = len(seq)
                if contig in alignments:
                    for long_contig, blast_alignments in alignments[contig].items():
                        if len(long_contig_aligned_counts[long_contig]) <= 1:
                            continue
                        intervals = list()
                        for info in blast_alignments:
                            intervals = putils.merge_alignments(intervals, info.qstart, info.qend)
                        aligned_length = putils.get_total_aligned_length(intervals)
                        percent_aligned = 100 * aligned_length // length
                        percent_aligned = min(percent_aligned, 100)
                        line = "\t".join([contig, long_contig, str(aligned_length), str(percent_aligned)])
                        print(line, file=handle)
                        n_edges += 1


def load_long_lengths(long_read_assembly, out_dir):
    """
    Get the contig lengths of the long-read assembly; check if a lenghts file is already generated for speed
    Args:
        long_read_assembly (str): path to long-read assembly fasta
        out_dir (str): path to desired output dir to save long-read contig lengths for later re-runs

    Returns:
        long_lengths (dict[str:int]): lengths of contigs

    """
    print("Loading long-read assembly information")
    lengths_file = os.path.join(out_dir, "long_read_lengths.tsv")
    if os.path.isfile(lengths_file):
        long_lengths = import_lengths(lengths_file)
    elif long_read_assembly.endswith(".tsv"):
        long_lengths = import_lengths(long_read_assembly)
    else:
        long_lengths = dict()
    if len(long_lengths) <= 1:
        long_contigs = utils.load_fasta(long_read_assembly)
        with open(lengths_file, "w+") as out_handle:
            for contig, seq in long_contigs.items():
                length = len(seq)
                long_lengths[contig] = length
                print(f"{contig}\t{length}", file=out_handle)
    return long_lengths


def load_contig_abund(filename):
    """
    Load contig abundances from metabat2 file
    Args:
        filename (str): path to file (i.e. contig_depth.txt)

    Returns:
        contig_depths (dict[str:float]): average abundance of each contig
    """
    contig_depths = dict()
    if filename is None:
        return contig_depths
    print(f"Loading contig read depths from {filename}")
    with open(filename) as input:
        for i, line in enumerate(input):
            if i > 0:
                cut = line.split("\t")
                contig_depths[cut[0]] = float(cut[2])
    return contig_depths


def keep_only_abundant_viruses(bins, abundances, min_depth=5):
    """

    Args:
        bins (dict[str:dict[str:int]]): contigs in each viral mag
        abundances (dict[str:float]): read coverage of each contig in assembly
        min_depth (float): minimum mean coverage of viral contigs in cluster

    Returns:
        curated_mags (dict[str:dict[str:int]]): contigs in each viral mag that passed thresholding
    """
    if len(abundances) == 0:
        return bins
    print("Removing low coverage viral contigs and bins")
    curated_mags = dict()
    for bin_name, contigs in bins.items():
        abund_total = 0
        length_total = 0
        for contig, seq in contigs.items():
            length = len(seq)
            abund_total += abundances[contig] * length
            length_total += length
        abundance = abund_total / length_total
        if abundance >= min_depth:
            curated_mags[bin_name] = contigs
    return curated_mags


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
    parser.add_argument("-o", "--output", required=True, type=str, dest="output",
                        help="Path to desired output directory")
    parser.add_argument("-b", "--bins", required=True, type=str, dest="bins",
                        help="Path to directory with clusters or bins to test")
    parser.add_argument("-s", "--bin-stats", required=False, type=str, dest="bin_stats",
                        help="Path to tabular file with completion/contamination information on the clusters")
    parser.add_argument("-l", "--long-assembly", required=True, type=str, dest="long_assembly",
                        help="Path to fasta file with trusted long-read assembly")
    parser.add_argument("-a", "--alignments", required=True, type=str, dest="alignments",
                        help="Path to blast alignment file in outfmt6 format (first column is the short-read contig)")
    parser.add_argument("--proxiphage_out", required=False, type=str, dest="proxiphage_out", default=None,
                        help="Path to proxiphage output directory with all the CheckV stats generated")
    parser.add_argument("--contig-abund", required=False, type=str, dest="contig_abund", default=None,
                        help="Path contig abundance/depth file from proxiwrap_out/contig_depth.txt")
    parser.add_argument("--min-depth", required=False, type=float, dest="min_depth", default=5,
                        help="Min average contig depth in mobile element cluster to be considered for validation")
    parser.add_argument("--hosts", required=False, type=str, dest="hosts_file", default=None,
                        help="Filtered master table with contig hosts")
    args = parser.parse_args(args=cli_args)
    return vars(args)


def main():
    """
    Launch analysis, print report, and make support network
    Returns: None
    """

    # load in arguments for run
    c_args = parse_args(__file__)
    bins_dir = c_args["bins"]
    bin_stats_file = c_args["bin_stats"]
    long_read_assembly = c_args["long_assembly"]
    blast_alignment = c_args["alignments"]
    proxiphage_out_dir = c_args["proxiphage_out"]
    out_dir = c_args["output"]
    contig_abund = c_args["contig_abund"]
    min_depth = c_args["min_depth"]
    hosts_file = c_args["hosts_file"]

    # load data for run
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    bins = dict()
    for name, contigs in utils.load_bin_set(bins_dir).items():
        if len(contigs) > 1:
            bins[name] = contigs
    abundances = load_contig_abund(contig_abund)
    bins = keep_only_abundant_viruses(bins, abundances, min_depth=min_depth)
    short_lengths = calc_contig_lengths_in_mags(bins)
    long_lengths = load_long_lengths(long_read_assembly, out_dir)
    seeds, stats = load_all_stats(bin_stats_file=bin_stats_file, proxiphage_dir=proxiphage_out_dir)
    curated_stats = curate_stats(stats)
    hosts = utils.load_host_connection_data(hosts_file)
    hits = load_blast(blast_alignment, short_lengths, long_lengths)

    # check if mags are supported
    plot_completion_contamination(mags=bins, hits=hits, outdir=out_dir, long_lengths=long_lengths)
    print_validation_report(mags=bins, hits=hits, stats=curated_stats, hosts=hosts, report_bad=True, outdir=out_dir)
    kept_bins = make_network_nodes(bins, hits, os.path.join(out_dir, "nodes.txt"))
    make_network_edges(kept_bins, hits, os.path.join(out_dir, "edges.txt"))


if __name__ == '__main__':
    """ Launch the script
    """
    main()
