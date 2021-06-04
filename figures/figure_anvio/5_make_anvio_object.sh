#!/usr/bin/env bash
source activate anvio-env
# make initial object, then quit
anvi-interactive -t viral_phylogeny_tree.txt -p anvio-object.db --title 'ProxiPhage_viral_bins_and_hosts' --manual

