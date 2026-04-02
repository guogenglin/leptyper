#!/bin/bash

# profiles
curl https://rest.pubmlst.org/db/pubmlst_leptospira_seqdef/schemes/2/profiles_csv \
 -o profiles.tsv
 
# loci
for locus in adk_2 glmU_2 icdA_2 lipL32_2 lipL41_2 mreA_2 pntA_2
do
    curl "https://rest.pubmlst.org/db/pubmlst_leptospira_seqdef/loci/${locus}/alleles_fasta" \
     -o ${locus}.fasta
done

# convert fasta to one line
for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"
done