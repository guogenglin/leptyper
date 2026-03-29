Leptyper is a command-line tool inspired by Kleborate (Mash-based species identification) and Kaptive (locus-based typing), designed for:

Species identification: Compare input genome assemblies (FASTA) against a built-in Mash sketch database.
Serovar (LPS locus) typing: Align reference LPS locus gene sequences to the input assemblies using minimap2, reporting the best-matching serovar along with coverage and identity.

Installation (from source)
python -m pip install -e .

External dependencies (must be in PATH):
mash
minimap2

Usage
leptyper -i *.fasta -o report.tsv
Output (TSV)

Columns include:
sample
predicted_species
mash_distance
mash_p_value
mash_shared_hashes
predicted_serovar
lps_percent_coverage
lps_percent_identity

Database files
By default, the tool looks for the following in the repository root:
leptospira_ref.msh (Mash species identification database; an error is raised if missing)
Leptospira_lps_locus_reference.gbk (LPS locus reference database; provided)

Leptyper (prototype)
Minimal command-line prototype for:
Species identification: Mash against leptospira_ref.msh
Serovar typing: Best-hit locus comparison against Leptospira_lps_locus_reference.gbk
Output: Kleborate-style TSV (one row per input assembly)

Requirements
Python 3.9+

External tools on PATH:
mash
minimap2

Install (editable)
python -m pip install -e .
Run
leptyper run \
  --assemblies path/to/a.fasta path/to/b.fasta \
  --mash-db leptospira_ref.msh \
  --lps-db Leptospira_lps_locus_reference.gbk \
  --out results.tsv

If the console script is not installed, you can also run:

python -m leptyper run --assemblies ... --mash-db ... --lps-db ... --out results.tsv
Output columns (TSV)

sample
predicted_species
mash_distance
mash_p_value
mash_shared_hashes
predicted_serovar
locus_hit
locus_coverage
locus_identity
locus_bitscore


**Citation**

