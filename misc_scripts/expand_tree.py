#!/usr/bin/env python3
import sys,os

bin_dir = sys.argv[1]
bins={}
for file in os.listdir(bin_dir):
	path=bin_dir+"/"+file
	bin = file.split(".")[0]
	bins[bin]=[]
	for line in open(path):
		if line[0]==">":
			name = line[1:-1]
			bins[bin].append(name)

for line in open(sys.argv[2]):
	data = line.strip()
s1 = data.split(",")
for i,s in enumerate(s1):
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
print(data)
