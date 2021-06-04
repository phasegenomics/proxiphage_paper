#!/usr/bin/env bash
source activate anvio-env
# re-launch the session
anvi-interactive -t viral_phylogeny_tree.txt -p anvio-object.db --title 'ProxiPhage_viral_bins_and_hosts' --manual
