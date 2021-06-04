#!/usr/bin/env python3
import sys
for line in open(sys.argv[1]):
    if line[0]==">":
        contig=line[1:-1].split()[0]
    else:
        l = len(line.strip())
        print(contig, l, sep="\t")
