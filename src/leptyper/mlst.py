# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:09:21 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

from pathlib import Path
from typing import Dict, List, Any

from .log import log
from .alignment import Assembly
from .utils import LeptyperError

def load_st_profiles(database_path: Path, gene_names: List[str], verbose: bool = False) -> List[tuple]:
    '''
    Load MLST profiles from a TSV file. Each line should have the format:
    ST\tgene1_allele\tgene2_allele\t...\tgeneN_allele\textra_info
    '''
    profiles = []
    log(f'Loading ST profiles from {database_path}', verbose = verbose)
    header_line = []
    with open(database_path, 'rt') as file:
        header_line = next(file).rstrip('\n').split('\t')[1:-1]  # read header line
        for line in file:
            parts = line.rstrip('\n').split('\t')
            st = int(parts[0])
            alleles = {gene: int(parts[i + 1]) for i, gene in enumerate(header_line)}
            
            extra_info = parts[len(gene_names) + 1] if len(parts) > len(gene_names) + 1 else None   # Optional extra info column
            
            profiles.append((st, alleles, extra_info))
    
    log(f'Loaded {len(profiles)} ST profiles', verbose = verbose)

    return profiles

def _get_allele_num(hit: Any) -> int:
    '''
    Extract allele number from name (e.g. 'glmU_1' -> 1)
    '''
    if not hit:
        return 0
    try:
        return int(hit.q.split('_')[-1])
    except (ValueError, IndexError):
        return 0

def solve_mlst(profiles: List[tuple], hits_per_gene: Dict[str, List[Any]], gene_names: List[str], 
               required_exact_matches: int = 3, verbose: bool = False) -> tuple:
    '''
    Logic to determine the ST based on hits.
    '''

    # 1. Get best hits for each gene (highest identity, then highest score)
    best_hits_per_gene = {}
    for gene in gene_names:
        hits = hits_per_gene[gene]  # There could be multiple hits for this gene
        if not hits:
            best_hits_per_gene[gene] = []
            log(f'{gene}: no hits', verbose = verbose)
            continue
        max_id = max(h.identity for h in hits)
        top_hits = [h for h in hits if h.identity == max_id]   # Multiple hits could have the same max identity
        max_score = max(h.score for h in top_hits)   # Filter further by score
        best_hits_per_gene[gene] = [h for h in top_hits if h.score == max_score]

        log(f'{gene}: retained {len(best_hits_per_gene[gene])} best hits ' f'(identity={max_id:.2f})', 
            verbose = verbose)

    # 2. Find matching profile
    best_st, best_extra, best_match_count = 0, '-', 0
    best_alleles = {gene: 0 for gene in gene_names}

    for st, alleles, extra in profiles:
        matches = 0
        for gene_name in gene_names:
            for h in best_hits_per_gene[gene_name]:
                if _get_allele_num(h) == alleles[gene_name]:
                    matches += 1
                    break
        
        if matches > best_match_count:
            best_st, best_alleles, best_extra, best_match_count = st, alleles, extra, matches
    
    log(f'Best ST candidate: ST{best_st} with {best_match_count} matches', verbose = verbose)

    # 3. Finalize allele strings and ST name
    exact_count, mismatch_count, final_alleles = 0, 0, {}
    
    for gene in gene_names:
        hits = best_hits_per_gene[gene]
        # Find if any hit matches the profile allele
        match_hit = next((h for h in hits if _get_allele_num(h) == best_alleles[gene]), None)
        
        if not hits:
            final_alleles[gene] = '-'
            mismatch_count += 1
            log(f'{gene}: no hit', verbose = verbose)
        else:
            hit = match_hit if match_hit else sorted(hits, key=lambda x: _get_allele_num(x))[0]
            a_num = _get_allele_num(hit)
            
            if hit.is_exact() and a_num == best_alleles[gene]:
                final_alleles[gene] = str(a_num)
                exact_count += 1
            else:
                final_alleles[gene] = f"{a_num}*"
                mismatch_count += 1

    # ST naming
    if exact_count < required_exact_matches:
        st_name = 'NA'
        best_extra = '-'
        log(f'Not enough exact matches ({exact_count}), setting ST = NA', verbose = verbose)

    elif mismatch_count == 0:
        st_name = f'ST{best_st}'
    else:
        st_name = f'ST{best_st}-{mismatch_count}LV'

    log(f'Final result: {st_name} ' f'(exact={exact_count}, mismatch={mismatch_count})', verbose = verbose)

    return st_name, best_extra, final_alleles

def lepto_mlst(assembly_obj: Assembly, verbose: bool = False) -> Dict[str, List[Any]]:
    '''
    Main function to perform MLST typing. Returns a dictionary with scheme names as keys and results as values.
    '''
    log('Starting MLST typing pipeline', verbose = verbose)

    results = {}
    base_db_path = Path(__file__).parents[0] / 'database'
    
    # Define the MLST schemes to run.
    schemes = ['mlst_scheme_1', 'mlst_scheme_2', 'mlst_scheme_3']
    log(f'MLST schemes to run: {schemes}', verbose = verbose)
    
    for scheme_name in schemes:
        log(f'Processing scheme: {scheme_name}', verbose = verbose)
        db_dir = base_db_path / scheme_name
        if not db_dir.exists():
            raise LeptyperError(f'MLST database directory {db_dir} not found for scheme {scheme_name}')
            
        # Load allele sequences and profiles 
        profiles_path = db_dir / 'profiles.tsv'
        # Load allele sequences  eg. glmU_1.fasta
        allele_files = list(db_dir.glob('*.fasta'))
        gene_names = [f.stem for f in allele_files]
        
        log(f'{scheme_name}: found {len(gene_names)} loci', verbose = verbose)
        # Load ST profiles: a list of tuples (ST, [allele1, allele2, ...], extra_info)
        profiles = load_st_profiles(profiles_path, gene_names, verbose = verbose)
        
        # Execute mapping: use -x asm5 as it is often best for intra-species gene mapping
        hits_per_gene = {gene: [] for gene in gene_names}
        # Map each gene's fasta file separately
        for gene in gene_names:
            log(f'{scheme_name}: mapping gene {gene}', verbose = verbose)

            gene_fasta = db_dir / f'{gene}.fasta'
            for aln in assembly_obj.map(gene_fasta, extra_args = '--end-bonus=10 --eqx -x asm5', 
                                        verbose = verbose, use_stdin = False):
                # Filters
                if aln.identity >= 90.0 and aln.coverage >= 80.0:
                    hits_per_gene[gene].append(aln)

        # Run MLST logic
        st, extra_info, allele_results = solve_mlst(profiles, hits_per_gene, gene_names, verbose = verbose)

        # 5. Format output: eg. {'ST': 'ST1', 'glmU': '1', 'pntA': '2*', 'extra_info': 'CC1'}
        scheme_result_list = {'ST': st}
        for gene in gene_names:
            scheme_result_list[gene] = allele_results.get(gene, '-')
        scheme_result_list['Extra_Info'] = extra_info if extra_info else '-'
        
        results[scheme_name] = scheme_result_list
        log(f'Finished scheme {scheme_name}', verbose = verbose)

    log('MLST typing completed', verbose = verbose)

    return results


