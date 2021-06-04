#!/usr/bin/env python3
import sys
threshold = float(".".join(sys.argv[1].split(".")[:-1]).split("_")[-1])
for line in open(sys.argv[1]):
	if "prophage host links were validated with" in line:
		cut = line.split()
		validated = int(cut[0])
		all = int(cut[3])
		percent = validated * 100.0 / all
		print(threshold, validated, all, percent, sep="\t")

