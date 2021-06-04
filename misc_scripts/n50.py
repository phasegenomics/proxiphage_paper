#!/usr/bin/env python3
import sys
lengths = list()
for line in open(sys.argv[1]):
    lengths.append(int(line.split("\t")[1]))
half = sum(lengths) / 2
lengths = sorted(lengths, reverse=False)
print(len(lengths))
tally = 0
for l in lengths:
    tally += l
    if tally >= half:
        print(l)
        break
print(sum(lengths))
