#!/usr/bin/env python3
import os
import sys
import shutil


def load_bin_set(bin_dir):
    """
    Load in sequence information from a directory of bins or clusters
    Args:
        bin_dir (str): path to bin or cluster directory
    Returns:
         bin_contigs (dict[str:dict[str:int]]): contigs present in each bin, and their lengths

    """
    print(f"Loading binning information from {bin_dir}")
    bin_contigs = dict()
    for filename in os.listdir(bin_dir):
        if filename.startswith(".") or "unbinned" in filename:
            continue
        # stores the contigs contained in each bin, and their lengths
        bin_name = os.path.splitext(filename)[0]
        bin_contigs[bin_name] = load_fasta(os.path.join(bin_dir, filename), quiet=True)
    return bin_contigs


def load_list(filename):
    print(f"loading from {filename}")
    data = set()
    for line in open(filename):
        data.add(line.strip())
    return data


def load_fasta(fasta_file, sequences=None, quiet=False, logger=None):
    """
    Load the sequence names and sequences from any fasta file

    Args:
        fasta_file (str): path to the fasta file for loading
        sequences (dict[str:str]): previous fasta entries to add to
        quiet (bool): should the program run with no messages? [default: False]
    Returns:
        contig_sequences (dict[str:str]): contig names pointing to their full sequence

    """
    if sequences is None:
        sequences = dict()
    if not quiet:
        print(f"Loading fasta file {fasta_file} into memory")
    with open(fasta_file, "r") as fasta_handle:
        for i, line in enumerate(fasta_handle):
            if i == 0 and not line.startswith(">"):
                raise ImportError(f"{fasta_file} does not appear to be a fasta file")
            if line.startswith(">"):
                contig = universal_contig_naming(line.rstrip())
                sequences[contig] = list()
            else:
                sequences[contig].append(line.rstrip())
    for contig, seqs in sequences.items():
        sequences[contig] = "".join(seqs)
    return sequences


def universal_contig_naming(contig_name):
    """
    Fix common errors and modifications in contig names introduced by various software
    Args:
        contig_name (str): name of contig
    Returns (str): fixed contig name

    """
    if contig_name.startswith(">"):
        contig_name = contig_name[1:]
    if "contig:" in contig_name:
        return contig_name.split()[1].split(":")[1]
    elif "=" in contig_name:
        return "_".join(contig_name.split()[0].split("="))
    else:
        return contig_name.split()[0]


def load_depths(filename):
    print(f"loading data from {filename}")
    contig_lengths = dict()
    contig_depths = dict()
    for i, line in enumerate(open(filename)):
        if i > 0:
            cut = line.strip().split("\t")
            contig = cut[0]
            l = int(cut[1])
            depth = float(cut[2])
            contig_lengths[contig] = l
            contig_depths[contig] = depth
    return contig_lengths, contig_depths


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


def load_checkv(filename):
    print(f"loading from {filename}")
    completions = dict()
    gene_count = dict()
    viral_gene_count = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.strip().split("\t")
            if i == 0:
                if cut[9] != "completeness":
                    raise ImportError(f"not a checkmv file: {cut[9]}")
                continue
            name = cut[0].split("|")[0]
            genes = int(cut[4])
            viral_genes = int(cut[5])
            comp = cut[9]
            if comp == "NA":
                comp = 0
            else:
                comp = float(comp)
            completions[name] = comp
            gene_count[name] = genes
            viral_gene_count[name] = viral_genes
    return completions, gene_count, viral_gene_count


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


def reverse_map_custers(clusters):
    contig_bins = dict()
    for name, contigs in clusters.items():
        for contig in contigs:
            contig_bins[contig] = name
    return contig_bins


def load_kraken(filename):
    print(f"loading {filename}")
    taxonomy = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            cut = line.split("\t")
            contig = cut[0]
            taxa = cut[1].strip().split(";")
            if len(taxa) > 2:
                taxonomy[contig] = taxa
    return taxonomy


def calculate_phyla(chosen_hosts, contig_taxonomy, lengths, n=4):
    phyla = dict()
    for host, contigs in chosen_hosts.items():
        tallies = dict()
        for contig in contigs:
            if contig in contig_taxonomy:
                taxa = contig_taxonomy[contig]
                if len(taxa) > 4:
                    phylum = taxa[3] + ":" + taxa[4]
                    if phylum not in tallies:
                        tallies[phylum] = 0
                    tallies[phylum] += lengths[host]
        best_phylum = None
        best_phylum_len = 0
        for phylum, l in tallies.items():
            if l > best_phylum_len:
                best_phylum_len = l
                best_phylum = phylum
        phyla[host] = best_phylum
    return phyla


def calculate_contig_gc(clusters):
    gcs = dict()
    for contigs in clusters.values():
        for contig, seq in contigs.items():
            all = 0
            gc = 0
            for c in seq:
                if c=="G" or c=="C":
                    gc +=1
                all += 1
            gc_percent = 100 * gc / all
            gcs[contig] = gc_percent
    return gcs


def calculate_cluster_gc(clusters):
    gcs = dict()
    for cluster, contigs in clusters.items():
        all = 0
        gc = 0
        for contig, seq in contigs.items():
            for c in seq:
                if c=="G" or c=="C":
                    gc +=1
                all += 1
        gc_percent = 100 * gc / all
        gcs[cluster] = gc_percent
    return gcs


def load_newick(filename):
    print(f"loading from {filename}")
    for line in open(filename):
        tree = line.strip()
        break
    return tree


def load_contigs_for_clusters(mag_names, filename):
    print(f"loading from {filename}")
    data = dict()
    for line in open(filename):
        cut = line.strip().split("\t")
        contig = cut[0]
        cluster = cut[1]
        if cluster in mag_names:
            if cluster not in data:
                data[cluster] = set()
            data[cluster].add(contig)
    return data


def load_gc(filename):
    print(f"loading from {filename}")
    data = dict()
    for line in open(filename):
        cut = line.strip().split("\t")
        contig = cut[0]
        gc = float(cut[2])
        data[contig] = gc
    return data


def subset_alignment(inputfile, outputfile, chosen_hosts):
    print(f"scanning {inputfile}")
    p = True
    with open(inputfile) as input:
        with open(outputfile, "w+") as output:
            for line in input:
                if line[0] == ">":
                    mag = line[1:].split()[0]
                    if mag in chosen_hosts:
                        p = True
                    else:
                        p = False
                if p:
                    print(line.rstrip(), file=output)


def run_anvio_fastree(alignement, tree):
    print("running anvio tree making")
    cmd = f"source activate anvio-env && anvi-gen-phylogenomic-tree -f {alignement} -o {tree}"
    print(f"running command {cmd}")
    code = os.system(cmd)
    print(f"finished with exit code {code}")


def main():
    # load all metadata for viral contigs
    chosen_vmags = load_list("chosen_vmags.list")
    chosen_vmags = load_contigs_for_clusters(chosen_vmags, "viral_mags.contigs.tsv")
    contig_gc = load_gc("contig_gc.tsv")
    contig_vmags = reverse_map_custers(chosen_vmags)
    contig_lengths, contig_depths = load_depths("contig_depth.txt")
    host_link_data = load_host_links("filtered_master_table.tsv")
    viral_contig_cv_completions, viral_contig_genes, viral_contig_viral_genes = load_checkv("viral_contigs.checkv.tsv")
    viral_mag_cv_completions, viral_mag_genes, viral_mag_viral_genes = load_checkv("viral_mags.checkv.tsv")
    viral_contig_completions, viral_contig_contaminations = load_hifi_validation("viral_contigs.hifi.tsv")
    viral_mag_completions, viral_mag_contaminations = load_hifi_validation("viral_mags.hifi.tsv")

    # load all metadata for bacterial MAGs
    chosen_hosts = load_list("chosen_hosts.list")
    chosen_hosts = load_contigs_for_clusters(chosen_hosts, "proxiwrap_clusters.contigs.tsv")
    host_gc = load_gc("proxiwrap_clusters.gc.tsv")
    contig_taxonomy = load_kraken("contig_kraken2_taxonomy.tsv")
    mag_lengths, mag_depths = load_depths("bin_depth.txt")
    mag_completions, mag_contaminations = load_checkm("proxiwrap_clusters.checkm.tsv")
    mag_phyla = calculate_phyla(chosen_hosts, contig_taxonomy, mag_lengths)

    print("making viral info table")
    with open("viral_info.tsv", "w+") as out:
        print("item_name", "vMAG", "Length", "Gene count", "Viral genes", "Read depth", "GC content",
              "Contig completion", "vMAG completion", sep="\t", file=out)
        for vmag, contigs in chosen_vmags.items():
            for contig in contigs:
                print(contig, vmag, contig_lengths[contig], viral_contig_genes[contig],
                      viral_contig_viral_genes[contig], contig_depths[contig], contig_gc[contig],
                      viral_contig_cv_completions[contig], viral_mag_cv_completions[vmag],
                      sep="\t", file=out)

    print("making host info table")
    with open("host_info.tsv", "w+") as out:
        print("Host", "Length",	"GC content", "Read depth",	"MAG completion", "MAG contamination", "Phylum",
              sep="\t", file=out)
        for host in chosen_hosts:
            print(host, mag_lengths[host], host_gc[host], mag_depths[host], mag_completions[host],
                  mag_contaminations[host], mag_phyla[host], sep="\t", file=out)

    print("making virus-host matrix")
    with open("host_connection.tsv", "w+") as out:
        host_list = list(chosen_hosts.keys())
        header = "item_name"
        for host in host_list:
            header += "\t" + host
        print(header, file=out)
        for contig in contig_vmags:
            line = contig
            for host in host_list:
                value = 0
                if contig in host_link_data:
                    if host in host_link_data[contig]:
                        value = host_link_data[contig][host]
                line += "\t" + str(value)
            print(line, file=out)

    print("making bin collection file")
    with open("bins_collection.txt", "w+") as out:
        for contig, mag in contig_vmags.items():
            print(contig, mag, sep="\t", file=out)

    print("making layer ordering file")
    subset_alignment("concatenated-proteins.fa", "chosen_hosts-proteins.fa", chosen_hosts)
    run_anvio_fastree("chosen_hosts-proteins.fa", "phylogenomic-tree.txt")
    host_newick_tree = load_newick("phylogenomic-tree.txt")

    custom_order = "names,vMAG,Length,Gene count,Viral genes,Read depth,GC content,Contig completion," \
                   "vMAG completion,"
    for vmag in chosen_vmags:
        custom_order += vmag + ","
    with open("ordering_info.tsv", "w+") as out:
        print("item_name", "data_type", "data_value", sep="\t", file=out)
        print("bact_tree", "newick", host_newick_tree, sep="\t", file=out)
        print("order", "basic", custom_order, sep="\t", file=out)



main()
