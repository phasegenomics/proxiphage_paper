#!/usr/bin/env python3
import sys
import pandas as pd

viruses = set()
for line in open(sys.argv[1]):
	viruses.add(line.strip())

hosts = set()
for line in open(sys.argv[2]):
	hosts.add(line.strip())

data = dict()
for i, line in enumerate(open(sys.argv[3])):
	cut = line.strip().split("\t")
	if i==0: continue
	virus = cut[0]
	host = cut[5]
	copies = cut[14]
	if virus in viruses and host in hosts:
		if virus not in data:
			data[virus] = dict()
		data[virus][host] = copies

hosts_list = list(hosts)
for i, (virus, subdata) in enumerate(data.items()):
	if i==0:
		line = "item_name"
		for host in hosts_list:
			line += "\t" + host
	else:
		line = virus
		for host in hosts_list:
			if host in subdata:
				value = subdata[host]
			else:
				value = "0"
			line += "\t" + value
	print(line)
