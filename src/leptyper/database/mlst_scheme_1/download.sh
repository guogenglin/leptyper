#!/bin/bash

# profiles
curl https://rest.pubmlst.org/db/pubmlst_leptospira_seqdef/schemes/1/profiles_csv \
 -o profiles.tsv
 
# loci
for locus in glmU_1 pntA_1 sucA_1 tpiA_1 pfkB_1 mreA_1 caiB_1
do
    curl "https://rest.pubmlst.org/db/pubmlst_leptospira_seqdef/loci/${locus}/alleles_fasta" \
     -o ${locus}.fasta
done

# convert fasta to one line
for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"
done