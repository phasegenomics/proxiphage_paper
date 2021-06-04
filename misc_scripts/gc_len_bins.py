#!/usr/bin/env python3
import sys
seq=""
for line in open(sys.argv[1]):
	if line[0]!=">":
		seq += line.strip()
gc=0
t=0
for c in seq:
	if c=="G" or c=="C":
		gc+=1
	t+=1

print(sys.argv[1].split("/")[-1].split(".")[0], len(seq), 100*gc/t, sep="\t")
