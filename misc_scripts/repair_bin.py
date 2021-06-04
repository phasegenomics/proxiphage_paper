#!/usr/bin/env python3
import sys

prev_contig = ""
for line in open(sys.argv[1]):
    if line[0]==">":
        contig = line
        if contig != prev_contig:
            prev_contig = contig
            print(line.strip())
    else:
        print(line.strip())

