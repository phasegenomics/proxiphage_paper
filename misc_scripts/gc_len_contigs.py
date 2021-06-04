#!/usr/bin/env python3
import sys
for line in open(sys.argv[1]):
	if line[0]==">":
		contig=line[1:-1]
	if line[0]!=">":
		seq = line.strip()
		gc = 0 
		tot = 0
		for c in seq:
			if c=="G" or c=="C":
				gc+=1
			tot += 1
		l = len(seq)
		gc=100*gc/tot
		print(contig, l, gc, sep="\t")
