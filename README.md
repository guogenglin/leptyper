# Leptyper
Leptyper is a command-line tool inspired by Kleborate (Mash-based species identification) and Kaptive (locus-based typing). It is designed for species identification, MLST, and in silico serotyping based on rfb locus sequences.

# Preface
I have not done much work on the algorithm itself. Instead, I combined the species identification and MLST logic from Kleborate with the serotyping logic from Kaptive. I learned a great deal from reading and understanding both tools, and I am deeply grateful to their development teams.

I will give a brief explanation of how this analysis works. For Leptospira, things are more complex, as this genus shares serovars and contains dozens of species and serovars. I will explain how this process works specifically for Leptospira.

# Dependencies
`mash`, `minimap2`, `biopython`, `python` >= 3.8

# Usage
```
leptyper [-h] -i INPUT [INPUT ...] [-o OUTPUT] [-t THREADS]
         [--min-cov MIN_COV] [--n-best N_BEST]
         [--percent-expected PERCENT_EXPECTED]
         [--no_cps_sequence] [-f] [-V] [--version]

Input and Output:
  -i, --input             Input FASTA files
  -r, --reference         GenBank format reference database (default: Leptospira_lps_locus_reference.gbk)
  -o, --output            Output file (default: leptyper_output.txt)

Parameters:
  -t, --threads           Number of alignment threads or 0 for all available CPUs (available cpus or 32(maximum))
  --min-cov               Minimum gene %coverage (blen/q_len*100) to be used for scoring (default: 50.0)
  --n-best                Number of best loci from the 1st round of scoring to be
                          fully aligned to the assembly (default: 2)
  --percent-expected      Typeable if >= % expected genes (default: 50.0)
  --no_singleton          Do not run alignment against the extra singleton database.
  --no_rfb_sequence       Suppress output of rfb sequence file
  -f, --figure            Export the gene structure map of rfb locus
  -V, --verbose           Print debug messages to stderr
  --version               Show version number and exit
  -h, --help              show this help message and exit
```

## Quick start

You can install this tool via pip, or run it directly as follows:

``` Python
PYTHONPATH=src(the root folder) python -m leptyper.cli -i *.fasta
```
I will also try to submit this tool to Bioconda so that it can be installed via conda in the future. I will update the documentation once it is available.

# Species identification
To date, there are 74 recognized species in the genus Leptospira, along with additional uncharacterized bacteria. I collected all available genome sequences from the GenBank public database (updated on 2026-03-06; 1,069 genomes in total). Considering that taxonomic annotations may be incorrect, I used NCBI’s best ANI match to determine species assignments. A threshold of 96% ANI was applied (no isolates fall between 95% and 96%, so lowering the threshold to 95% would not change the results). If the ANI is below 96%, the isolate is assigned as Leptospira sp..

Unfortunately, only 34 distinct species could be confidently identified, as some species currently have no available genome sequences.

After assigning species labels to these genomes, I built a Mash database. Species identification is then performed using Mash, with Mash distance used as the criterion for assigning species.

# MLST sequence type
The MLST database were obtained from [pubMLST](https://pubmlst.org/), update: 2026-04-01
  
  There are 3 different MLST schemes for Leptospira, all of them will be included and screened by Leptyper.
  
  *glmU_1*, *pntA_1*, *sucA_1*, *tpiA_1*, *pfkB_1*, *mreA_1*, *caiB_1* for MLST scheme 1.
  
  *adk_2*, *glmU_2*, *icdA_2*, *lipL32_2*, *lipL41_2*, *mreA_2*, *pntA_2* for MLST scheme 2.
  
  *adk_3*, *icdA_3*, *lipL41_3*, *lipL32_3*, *rrs2_3*, *secY_3* for MLST scheme 3.

I have included a Bash script in the database that can be used to update the MLST database. 

The result will returns either the exact allele number or the closest matching allele number (marked with an asterisk). These results are then used to determine the Sequence Type (ST). If some genes do not have an exact match, the resulting ST will include a suffix indicating the number of variant loci. If too many genes lack exact matches, no ST will be assigned.

# Serovar prediction
I directly adopted the logic from the latest version of Kaptive to predict Leptospira serovars. The main effort in this part is the construction of a reference database. However, although many Leptospira isolates are annotated with serovars, I do not consider all of these annotations reliable. I have observed cases where isolates labeled with the same serovar show substantial differences in their rfb loci. In such cases, these serovars are excluded from the database. In addition, many serovars are represented by only a single isolate, so their assignments still require further validation.

For some serovars, multiple isolates are available but show minor variation. In these cases, I prioritize predominant patterns—for example, when more than half to two-thirds of isolates share the same rfb locus, I treat this as the correct representation. In most such cases, the total number of isolates is small (often 3–4), so a single divergent isolate is likely to be misannotated. I also assign higher confidence to isolates described as type strains in published studies. 

In total, 49 distinct rfb loci from different serovars were collected.

The original version of Kaptive used BLAST to search serovar-determinant loci in genome sequences, identifying full-length loci in query genomes and assigning serovars based on the best match. In the latest version, BLAST has been replaced by minimap2, and a gene-first approach has been introduced. This workflow first scans reference genes for each serovar, selects the top two (user-adjustable) best matches, and then compares full locus alignment to determine the closest match.

I have made only minimal modifications to this component, so the results are largely produced using the same underlying logic.

Moving forward, my goal is to go beyond the current rfb reference database by constructing an isolate-specific rfb dataset. For isolates that cannot be serotyped, I aim to identify whether they share identical or highly similar rfb loci with known isolates.

Finally, it should be noted that for some serovars—such as Icterohaemorrhagiae and Copenhageni—the rfb loci are highly similar, which may introduce ambiguity and bias in serovar assignments.

Here is the explanation of the problems mark:
Problems: characters indicating issues with the locus match. An absence of any such characters indicates a very good match.
```
? = the match was not in a single piece, possible due to a poor match or discontiguous assembly.The number of pieces will be indicated with an integer.. 
- = one or more genes expected in the locus were not found.
+ = one or more extra genes were found in the locus.
* = one or more expected genes was found but with identity below the minimum threshold (default threshold is 80%)
! = one or more locus genes is truncated
```
However, based on my experience, these markers are often over-interpreted or overemphasized. In many cases, I think we do not need to be taken so seriously, instead of focusing on coverage and identity (That said, because serovar-determinant loci within the same serovars can vary substantially, while different serotypes may sometimes appear highly similar, determining an appropriate threshold will raise another complex issue.)

In this tool, a sequence is considered typable if it has both identity and coverage ≥ 95%. Note that this criterion reflects sequence similarity only and does not necessarily correspond to the exact serovar.

## Adjustment of the Pending Logic
Kaptive provides two parameters—weight_metric and score_metric—to evaluate gene alignment results. However, its default approach relies solely on the total alignment score (AS) to determine the best-matching locus.

During testing, I identified a limitation in this strategy. In some cases, irrelevant matches may be detected within the input genome, resulting in multiple small alignment fragments. The combined AS of these fragments can exceed that of the true match, leading to incorrect prioritization.

For example, an incorrect match may produce: `AS = 1815 + 2745 + 1952 + 3098 + 2844 + 3315`

while the true match has: `AS = 10861`

In such cases, the true match may be incorrectly filtered out.

To address this issue, I revised the scoring logic by incorporating identity and coverage as weighting factors:

$$
score = \sum_{i=1}^{n} AS_i \cdot \left(\frac{mlen_i}{blen_i}\right) \cdot \left(\frac{q_{end,i} - q_{start,i}}{q_{len,i}}\right)
$$

## Singleton Database
An additional singleton database has been constructed. For genomes that do not yield a typable result in the reference database, their rfb sequences are collected and stored.

Users can optionally query this supplementary database (enabled by default) to check for potential matches. In total, 136 rfb loci have been included, combined with the reference database, this enables coverage all rfb loci from Leptospira genomes available in the NCBI RefSeq database (2026-03-06).

# Output
A tab-separated output file will be generated, so as the rfb locus recognized in every query genomes, there may be also a figure generated by dnafeatures, I'm still thinking to do it or not.
example:
```
Isolate	Species	Hit_strength	Distance	P-val	Shared_hashes	mlst_scheme_1	caiB_1	glmU_1	mreA_1	pfkB_1	pntA_1	sucA_1	tpiA_1	Extra_Info	mlst_scheme_2	adk_2	glmU_2	icdA_2	lipL32_2	lipL41_2	mreA_2	pntA_2	Extra_Info	mlst_scheme_3	adk_3	icdA_3	lipL32_3	lipL41_3	rrs2_3	secY_3	Extra_Info	Serovar	Confidence	Problems	Identity	Coverage	Length discrepancy	Expected genes in locus	Expected genes in locus, details	Missing expected genes	Other genes in locus	Other genes in locus, details	Expected genes outside locus	Expected genes outside locus, details	Other genes outside locus	Other genes outside locus, details	Truncated genes, details	Extra genes, details
GCF_000007685.1	Leptospira_interrogans	strong	0.0	0.0	1000/1000	ST17	8	1	4	10	1	2	2	-	ST47	5	1	3	2	5	4	2	0	ST2	1	1	2	2	1	1	-	Copenhageni	Typeable	+!	100.00	100.00	0 bp	90 / 90 (100.00%)	Copenhageni_0001,100.00%,100.00%;Copenhageni_0002_fcl,100.00%,100.00%;Copenhageni_0003,100.00%,100.00%;Copenhageni_0004,100.00%,100.00%;Copenhageni_0005,100.00%,100.00%;Copenhageni_0006,100.00%,100.00%;Copenhageni_0007,100.00%,100.00%;Copenhageni_0008,100.00%,100.00%;Copenhageni_0009,100.00%,100.00%;Copenhageni_0010,100.00%,100.00%;Copenhageni_0011,100.00%,100.00%;Copenhageni_0012,100.00%,100.00%;Copenhageni_0013_tktC,100.00%,100.00%;Copenhageni_0014,100.00%,100.00%;Copenhageni_0015,100.00%,100.00%;Copenhageni_0016_mutT,100.00%,100.00%;Copenhageni_0017,100.00%,100.00%;Copenhageni_0018,100.00%,100.00%;Copenhageni_0019,100.00%,100.00%;Copenhageni_0020_per,100.00%,100.00%;Copenhageni_0021,100.00%,100.00%;Copenhageni_0022_fcI,100.00%,100.00%;Copenhageni_0023,100.00%,100.00%;Copenhageni_0024,100.00%,100.00%;Copenhageni_0025,100.00%,100.00%;Copenhageni_0026,100.00%,100.00%;Copenhageni_0027,100.00%,100.00%;Copenhageni_0028,100.00%,100.00%;Copenhageni_0029,100.00%,100.00%;Copenhageni_0030_neuA,100.00%,100.00%;Copenhageni_0031,100.00%,100.00%;Copenhageni_0032,100.00%,100.00%;Copenhageni_0033,100.00%,100.00%;Copenhageni_0034_neuB,100.00%,100.00%;Copenhageni_0035,100.00%,100.00%;Copenhageni_0036,100.00%,100.00%;Copenhageni_0037,100.00%,100.00%;Copenhageni_0038_sas,100.00%,100.00%;Copenhageni_0039,100.00%,100.00%;Copenhageni_0040,100.00%,100.00%;Copenhageni_0041,100.00%,100.00%;Copenhageni_0042_moaA,100.00%,100.00%;Copenhageni_0043,100.00%,100.00%;Copenhageni_0044,100.00%,100.00%;Copenhageni_0045,100.00%,100.00%;Copenhageni_0046_nagB,100.00%,100.00%;Copenhageni_0047,100.00%,100.00%;Copenhageni_0048,100.00%,100.00%;Copenhageni_0049,100.00%,100.00%;Copenhageni_0050,100.00%,100.00%;Copenhageni_0051,100.00%,100.00%;Copenhageni_0052,100.00%,100.00%;Copenhageni_0053,100.00%,100.00%;Copenhageni_0054,100.00%,100.00%;Copenhageni_0055,100.00%,100.00%;Copenhageni_0056,100.00%,100.00%;Copenhageni_0057_rfbG,100.00%,100.00%;Copenhageni_0058,100.00%,100.00%;Copenhageni_0059_galE,100.00%,100.00%;Copenhageni_0060,100.00%,100.00%;Copenhageni_0061,100.00%,100.00%;Copenhageni_0062,100.00%,100.00%;Copenhageni_0063,100.00%,100.00%;Copenhageni_0064,100.00%,100.00%;Copenhageni_0065,100.00%,100.00%;Copenhageni_0066,100.00%,100.00%;Copenhageni_0067,100.00%,100.00%;Copenhageni_0068,100.00%,100.00%;Copenhageni_0069_rffE,100.00%,100.00%;Copenhageni_0070_wcaI,100.00%,100.00%;Copenhageni_0071,100.00%,100.00%;Copenhageni_0072,100.00%,100.00%;Copenhageni_0073,100.00%,100.00%;Copenhageni_0074_hlpA,100.00%,100.00%;Copenhageni_0075,100.00%,100.00%;Copenhageni_0076,100.00%,100.00%;Copenhageni_0077,100.00%,100.00%;Copenhageni_0078,100.00%,100.00%;Copenhageni_0079,100.00%,100.00%;Copenhageni_0080,100.00%,100.00%;Copenhageni_0081,100.00%,100.00%;Copenhageni_0082_rfbC,100.00%,100.00%;Copenhageni_0083_rfbD,100.00%,100.00%;Copenhageni_0084_rfbB,100.00%,100.00%;Copenhageni_0085_rfbA,100.00%,100.00%;Copenhageni_0086,100.00%,100.00%;Copenhageni_0087_rfbF,100.00%,100.00%;Copenhageni_0088,100.00%,100.00%;Copenhageni_0089,100.00%,100.00%;Copenhageni_0090,100.00%,100.00%		3	Yeoncheon_0067,100.00%,50.00%,truncated;Lai_0075,80.00%,100.00%;Pyrogenes_0033,100.00%,4.54%,truncated	0 / 90 (0.00%)		18	Pomona_0035,97.70%,25.92%;Pomona_0035,98.85%,25.92%;Pomona_0035,98.85%,25.92%;Pomona_0035,97.70%,25.92%;Pomona_0035,98.85%,25.92%;Pomona_0035,97.70%,25.92%;Pingchang_0054,92.35%,55.42%;Panama_0012,83.51%,35.47%;Pingchang_0055,72.73%,21.53%;Bratislava_0062,96.08%,137.50%;Semaranga_0047,72.32%,46.17%;Bananal_0028,92.59%,96.43%;Cynopteri_0028,73.17%,128.12%;Cynopteri_0028,73.17%,128.12%;Panama_0043,78.79%,8.91%;Hardjo-bovis_0064,93.55%,53.38%;Bataviae_0059,79.63%,33.33%;Pomona_0030,100.00%,100.00%	Yeoncheon_0067,100.00%,50.00%,truncated;Pyrogenes_0033,100.00%,4.54%,truncated	
GCF_000013945.1	Leptospira_borgpetersenii	strong	0.0	0.0	1000/1000	ST152	29	26	29	39	30	28	35	-	ST175	29	54	66	49	67	31	55	12	ST145-1LV	57	54	35	39*	20	47	-	Hardjo-bovis	Typeable	+!	99.59	100.00	0 bp	82 / 82 (100.00%)	Hardjo-bovis_0001,100.00%,100.00%;Hardjo-bovis_0002,100.00%,100.00%;Hardjo-bovis_0003,100.00%,100.00%;Hardjo-bovis_0004,100.00%,100.00%;Hardjo-bovis_0005,100.00%,100.00%;Hardjo-bovis_0006,100.00%,100.00%;Hardjo-bovis_0007,100.00%,100.00%;Hardjo-bovis_0008,100.00%,100.00%;Hardjo-bovis_0009,100.00%,100.00%;Hardjo-bovis_0010,100.00%,100.00%;Hardjo-bovis_0011,100.00%,100.00%;Hardjo-bovis_0012,100.00%,100.00%;Hardjo-bovis_0013,100.00%,100.00%;Hardjo-bovis_0014,100.00%,100.00%;Hardjo-bovis_0015,100.00%,100.00%;Hardjo-bovis_0016,100.00%,100.00%;Hardjo-bovis_0017,100.00%,100.00%;Hardjo-bovis_0018,100.00%,100.00%;Hardjo-bovis_0019,100.00%,100.00%;Hardjo-bovis_0020,100.00%,100.00%;Hardjo-bovis_0021,100.00%,100.00%;Hardjo-bovis_0022,100.00%,100.00%;Hardjo-bovis_0023,100.00%,100.00%;Hardjo-bovis_0024,100.00%,100.00%;Hardjo-bovis_0025,100.00%,100.00%;Hardjo-bovis_0026,100.00%,100.00%;Hardjo-bovis_0027,100.00%,100.00%;Hardjo-bovis_0028,100.00%,100.00%;Hardjo-bovis_0029,100.00%,100.00%;Hardjo-bovis_0030,100.00%,100.00%;Hardjo-bovis_0031_pgl,100.00%,100.00%;Hardjo-bovis_0032,100.00%,100.00%;Hardjo-bovis_0033,100.00%,100.00%;Hardjo-bovis_0034,100.00%,100.00%;Hardjo-bovis_0035_kdsB,100.00%,100.00%;Hardjo-bovis_0036,100.00%,100.00%;Hardjo-bovis_0037,100.00%,100.00%;Hardjo-bovis_0038,100.00%,100.00%;Hardjo-bovis_0039,66.67%,230.00%;Hardjo-bovis_0040,100.00%,100.00%;Hardjo-bovis_0041,100.00%,100.00%;Hardjo-bovis_0042,100.00%,100.00%;Hardjo-bovis_0043,100.00%,100.00%;Hardjo-bovis_0044,100.00%,100.00%;Hardjo-bovis_0045,100.00%,100.00%;Hardjo-bovis_0046,100.00%,100.00%;Hardjo-bovis_0047,100.00%,100.00%;Hardjo-bovis_0048_hisH,100.00%,100.00%;Hardjo-bovis_0049,100.00%,100.00%;Hardjo-bovis_0050,100.00%,100.00%;Hardjo-bovis_0051,100.00%,100.00%;Hardjo-bovis_0052,100.00%,100.00%;Hardjo-bovis_0053,100.00%,100.00%;Hardjo-bovis_0054,100.00%,100.00%;Hardjo-bovis_0055,100.00%,100.00%;Hardjo-bovis_0056_wecB,100.00%,100.00%;Hardjo-bovis_0057,100.00%,100.00%;Hardjo-bovis_0058,100.00%,100.00%;Hardjo-bovis_0059,100.00%,100.00%;Hardjo-bovis_0060,100.00%,100.00%;Hardjo-bovis_0061,100.00%,100.00%;Hardjo-bovis_0062,100.00%,100.00%;Hardjo-bovis_0063,100.00%,100.00%;Hardjo-bovis_0064,100.00%,100.00%;Hardjo-bovis_0065,100.00%,100.00%;Hardjo-bovis_0066,100.00%,100.00%;Hardjo-bovis_0067,100.00%,100.00%;Hardjo-bovis_0068,100.00%,100.00%;Hardjo-bovis_0069,100.00%,100.00%;Hardjo-bovis_0070,100.00%,100.00%;Hardjo-bovis_0071,100.00%,100.00%;Hardjo-bovis_0072,100.00%,100.00%;Hardjo-bovis_0073,100.00%,100.00%;Hardjo-bovis_0074,100.00%,100.00%;Hardjo-bovis_0075_rfbC,100.00%,100.00%;Hardjo-bovis_0076_rfbD,100.00%,100.00%;Hardjo-bovis_0077_rfbB,100.00%,100.00%;Hardjo-bovis_0078_rfbA,100.00%,100.00%;Hardjo-bovis_0079,100.00%,100.00%;Hardjo-bovis_0080,100.00%,100.00%;Hardjo-bovis_0081,100.00%,100.00%;Hardjo-bovis_0082,100.00%,100.00%		1	Ricardi_0039,68.18%,6.87%,truncated	3 / 82 (3.66%)	Hardjo-bovis_0038,99.73%,100.00%;Hardjo-bovis_0038,100.00%,28.57%,truncated;Hardjo-bovis_0038,100.00%,100.00%;Hardjo-bovis_0038,99.73%,100.00%;Hardjo-bovis_0038,100.00%,100.00%;Hardjo-bovis_0064,97.74%,100.00%;Hardjo-bovis_0064,93.98%,100.00%;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,97.22%,56.69%,truncated	34	Bataviae_0059,92.05%,132.69%;Bataviae_0059,92.05%,132.69%;Tarassov_0032,98.08%,100.00%;Tarassov_0032,98.63%,100.00%;Tarassov_0032,98.35%,100.00%;Tarassov_0032,98.35%,100.00%;Tarassov_0032,98.90%,100.00%;Tarassov_0032,96.15%,28.57%;Tarassov_0034,75.00%,2.52%;Pingchang_0054,85.63%,109.34%;Pingchang_0054,85.63%,109.34%;Pingchang_0054,85.02%,109.34%;Pingchang_0054,85.02%,109.34%;Pingchang_0054,85.02%,109.04%;Pingchang_0054,85.32%,109.04%;Pomona_0038,90.79%,226.03%;Pomona_0038,89.13%,63.01%;Sejroe_0053,99.07%,100.00%;Bim_0018,85.96%,43.40%;Guaricura_0021,40.43%,100.70%;Tarassov_0016,92.73%,80.00%;Sejroe_0068,99.07%,99.07%;Balcanica_0020,98.68%,99.34%;Pingchang_0082,66.67%,388.89%;Sokoine_0078,84.38%,32.33%;Pomona_0031,28.00%,102.08%;Pomona_0031,28.00%,102.08%;Pomona_0031,28.00%,104.17%;Pomona_0028,79.49%,263.41%;Pomona_0028,82.05%,263.41%;Tarassov_0033,32.00%,57.25%;Pomona_0035,80.95%,16.49%;Tarassov_0015,91.43%,23.49%;Balcanica_0019,100.00%,100.00%	Ricardi_0039,68.18%,6.87%,truncated;Hardjo-bovis_0038,100.00%,28.57%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,87.83%,90.55%,truncated;Hardjo-bovis_0065,97.22%,56.69%,truncated	
```
  
