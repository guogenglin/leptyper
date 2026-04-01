#!/bin/bash

# profiles
curl https://rest.pubmlst.org/db/pubmlst_leptospira_seqdef/schemes/1/profiles_csv \
 -o profiles.tsv
 
# loci
for locus in adk_3 icdA_3 lipL41_3 lipL32_3 rrs2_3 secY_3
do
    curl "https://rest.pubmlst.org/db/pubmlst_leptospira_seqdef/loci/${locus}/alleles_fasta" \
     -o ${locus}.fasta
done

# convert fasta to one line
for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"
done