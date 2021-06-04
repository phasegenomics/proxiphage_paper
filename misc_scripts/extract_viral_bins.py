#!/usr/bin/env python

import sys,os


viral_fasta = sys.argv[1]
bins_dir = sys.argv[2]
out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

viruses=set()
for line in open(viral_fasta):
    if line[0]==">":
        viruses.add(line[1:-1])

for filename in os.listdir(bins_dir):
    filepath = bins_dir + "/" + filename
    outpath = out_dir + "/" + filename
    out = None
    for line in open(filepath):
        if line[0]==">":
            contig = line[1:-1]
            if contig.startswith("CL"):
                contig = contig.split()[1].split(":")[-1]
            contig_line = ">"+contig+"\n"
        else:
            if contig in viruses:
                if out is None:
                    out = open(outpath, "w+")
                out.write(contig_line)
                out.write(line)
    if out is not None:
        out.close()



