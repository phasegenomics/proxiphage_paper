#!/usr/bin/env bash

# annotate genes
for i in `ls *fasta | awk 'BEGIN{FS=".fasta"}{print $1}'`; do anvi-gen-contigs-database -T 8 -f $i.fasta -o $i.db; anvi-run-hmms -c $i.db; done

# make mapping file
echo -e "name\tcontigs_db_path" > genomes.txt
for i in *db; do echo -e "${i%.*}\t${i}" >> genomes.txt; done

# find common marker genes
anvi-get-sequences-for-hmm-hits --external-genomes genomes.txt --hmm-source Bacteria_71 --list-available-gene-names

# extract and contatenate marker genes
anvi-get-sequences-for-hmm-hits --external-genomes genomes.txt -o concatenated-proteins.fa --hmm-source Bacteria_71 --gene-names Ribosom_S12_S23,Ribosomal_L1,Ribosomal_L13,Ribosomal_L14,Ribosomal_L16,Ribosomal_L17,Ribosomal_L18p,Ribosomal_L19,Ribosomal_L2,Ribosomal_L20,Ribosomal_L21p,Ribosomal_L22,Ribosomal_L23,Ribosomal_L27,Ribosomal_L27A,Ribosomal_L28,Ribosomal_L29,Ribosomal_L3,Ribosomal_L32p,Ribosomal_L35p,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6,Ribosomal_L9_C,Ribosomal_S10,Ribosomal_S11,Ribosomal_S13,Ribosomal_S15,Ribosomal_S16,Ribosomal_S17,Ribosomal_S19,Ribosomal_S2,Ribosomal_S20p,Ribosomal_S3_C,Ribosomal_S6,Ribosomal_S7,Ribosomal_S8,Ribosomal_S9 --return-best-hit --get-aa-sequences --concatenate

# make newick tree
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa -o phylogenomic-tree.txt

