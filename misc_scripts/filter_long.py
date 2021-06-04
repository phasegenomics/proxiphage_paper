#!/usr/bin/env python3
import sys

for line in open(sys.argv[1]):
    if line[0]==">":
        name=line.strip()
    else:
        l = len(line.strip())
        if l >= int(sys.argv[2]):
            print(name)
            print(line.strip())
