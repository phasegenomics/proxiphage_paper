#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
import proxiphage.proxiphage_utilities as utils
import proxiphage.proxiphage_plotting as plotting
import numpy

def load_lengths(filename):
    data = dict()
    with open(filename) as input:
        for i, line in enumerate(input):
            if i > 0:
                cut = line.strip().split("\t")
                name = cut[0]
                length = int(cut[3])
                data[name] = length
    return data


lengths1 = load_lengths(sys.argv[1])
lengths2 = load_lengths(sys.argv[2])

data1 = list(lengths1.values())
data2 = list(lengths2.values())

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_title("Length distribution of viral contigs (red) and MAGs (blue)")
ax.set_xlabel("Length (bp)")
ax.set_ylabel("Frequency (counts)")
ax.grid(axis="both", ls="-", alpha=0.1, c="k", which="major")
ax.grid(axis="both", ls="--", alpha=0.1, c="k", which="minor")
ax.set_yscale('log')
ax.set_xlim(0, 200000)
bins = numpy.linspace(0, 200000, 21)
ax.hist(data2, bins=bins, color="b", alpha=0.7, label="Viral MAGs")
ax.hist(data1, bins=bins, color="r", alpha=0.7, label="Viral Contigs")
plt.legend(loc='upper right')
plt.xticks(numpy.linspace(0, 200000, 11))

plt.savefig("lengths_histogram.png", bbox_inches='tight', dpi=300)
plt.close()
