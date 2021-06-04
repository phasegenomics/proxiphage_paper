#!/usr/bin/env python3
import sys
import os

bins = {}
for line in open("viral_mags.contigs.tsv"):
	cut = line.strip().split("\t")
	contig = cut[0]
	vmag = cut[1]
	if vmag not in bins:
		bins[vmag] = set()
	bins[vmag].add(contig)

data = ""
for line in open("viral_phylogeny_tree.raw.txt"):
	data = line.strip()
	break

with open("viral_phylogeny_tree.txt", "w+") as output:
	s1 = data.split(",")
	for i, s in enumerate(s1):
		if "vMAG" in s:
			s2 = s.split(":")
			mag_start = s2[0].split("@")[0].split("(")
			mag = mag_start[-1]
			contigs = bins[mag]
			fixed = f"({','.join(contigs)})"
			mag_start[-1] = fixed
			s2[0] = "(".join(mag_start)
			s = ":".join(s2)
		s1[i] = s
	data = ",".join(s1)
	print(data, file=output)
