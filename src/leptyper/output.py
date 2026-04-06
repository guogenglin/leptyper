# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:14:05 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

from __future__ import annotations

from pathlib import Path

from .log import log, warning

def fmt_pct(value: float) -> str:
    """Format a percentage, capping at 100.00."""
    return f"{min(value, 100.0):.2f}"


def is_low_quality(serovar) -> bool:
    return serovar.percent_coverage < 50 and serovar.percent_identity < 50


def serovar_score(serovar) -> float:
    return len(serovar.best_match) * serovar.percent_identity * serovar.percent_coverage


def write_rfb_sequence(serovar, root_dir: Path) -> None:
    fasta_dir = root_dir / 'rfb_locus'
    fasta_dir.mkdir(exist_ok=True)
    fasta_path = fasta_dir / f'{serovar.sample_name}_rfb_locus_sequence.fasta'
    with open(fasta_path, 'wt') as handle:
        handle.write(''.join(p.format('fna') for p in serovar.pieces))


def select_rfb_serovar(sv0, sv1):
    """
    Return the serovar whose sequence should be written, or None if both are too low quality.
    Priority: Typable > partial hit > score comparison > give up.
    """
    if sv0 and sv0.confidence == 'Typable':
        return sv0
    if sv1 and sv1.confidence == 'Typable':
        return sv1

    low0 = sv0 and is_low_quality(sv0)
    low1 = sv1 and is_low_quality(sv1)

    # Only one serovar present
    if sv0 and not sv1:
        return None if low0 else sv0
    if sv1 and not sv0:
        return None if low1 else sv1

    # Both present: prefer the non-low-quality one, or pick by score
    if low0 and not low1:
        return sv1
    if low1 and not low0:
        return sv0
    if not low0 and not low1:
        return sv0 if serovar_score(sv0) >= serovar_score(sv1) else sv1

    return None  # Both low quality


def build_headers(result_file: dict) -> list[str]:
    headers = ['Isolate', 'Species', 'Hit_strength', 'Distance', 'P-val', 'Shared_hashes',
               'mlst_scheme_1']
    for scheme in ['mlst_scheme_1', 'mlst_scheme_2', 'mlst_scheme_3']:
        if scheme != 'mlst_scheme_1':
            headers.append(scheme)
        headers += list(result_file['MLST'][scheme].keys())[1:]
    headers += [
        'Serovar', 'Confidence', 'Problems', 'Identity', 'Coverage',
        'Singleton', 'Sin_Problems', 'Sin_Identity', 'Sin_Coverage',
        'Length discrepancy', 'Expected genes in locus', 'Expected genes in locus, details',
        'Missing expected genes', 'Other genes in locus', 'Other genes in locus, details',
        'Expected genes outside locus', 'Expected genes outside locus, details',
        'Other genes outside locus', 'Other genes outside locus, details',
        'Truncated genes, details', 'Extra genes, details',
    ]
    return headers


def build_serovar_fields(sv0, sv1, min_locus_id: float = 95.0, min_locus_cov: float = 95.0) -> list[str]:
    """Build the 21 serovar-related output fields."""
    na = 'n/a'

    if sv0 is None and (sv1 is None or sv1.confidence == 'Untypeable'):
        return [na] * 21

    # --- Primary serovar (sv0) ---
    sv0_fields = (
        [sv0.output_name(min_locus_id, min_locus_cov), sv0.confidence, sv0.problems,
         fmt_pct(sv0.percent_identity), fmt_pct(sv0.percent_coverage)]
        if sv0 else [na] * 5
    )

    # --- Singleton serovar (sv1) ---
    sv1_fields = (
        [sv1.output_name(min_locus_id, min_locus_cov), sv1.problems,
         fmt_pct(sv1.percent_identity), fmt_pct(sv1.percent_coverage)]
        if (sv1 and sv1.confidence == 'Typeable') else [na] * 4
    )

    # --- Length discrepancy ---
    active = next((sv for sv in (sv0, sv1) if sv and sv.confidence == 'Typable'), None)
    length_discrepancy = (
        f"{active.__len__() - len(active.best_match)} bp"
        if active and len(active.pieces) == 1
        else na
    )

    # --- Gene stats (from sv0 if available) ---
    if sv0 is not None:
        n_expected = len(sv0.best_match.genes)
        n_inside = len({g.gene.name for g in sv0.expected_genes_inside_locus})
        n_outside = len({g.gene.name for g in sv0.expected_genes_outside_locus})
        gene_fields = [
            f"{n_inside} / {n_expected} ({100 * n_inside / n_expected:.2f}%)",
            ';'.join(map(str, sv0.expected_genes_inside_locus)),
            ';'.join(sv0.missing_genes),
            str(len(sv0.unexpected_genes_inside_locus)),
            ';'.join(map(str, sv0.unexpected_genes_inside_locus)),
            f"{n_outside} / {n_expected} ({100 * n_outside / n_expected:.2f}%)",
            ';'.join(map(str, sv0.expected_genes_outside_locus)),
            str(len(sv0.unexpected_genes_outside_locus)),
            ';'.join(map(str, sv0.unexpected_genes_outside_locus)),
            ';'.join(str(g) for g in sv0 if g.phenotype == 'truncated'),
            ';'.join(map(str, sv0.extra_genes)),
        ]
    else:
        gene_fields = [na] * 11

    return sv0_fields + sv1_fields + [length_discrepancy] + gene_fields


def generate_output(result_file: dict, output_path: str, min_locus_cov: float = 95.0, min_locus_id: float = 95.0, 
                    no_rfb_sequence: bool = False, figure: bool = False, verbose: bool = False) -> None:
    
    sv0, sv1 = result_file['Serovar']
    isolate = result_file['Isolate']
    species = result_file['Species']

    log(f'Generating output for {isolate.name}.', verbose)
    # Write header row if file does not yet exist
    if not Path(output_path).is_file():
        with open(output_path, 'wt') as fh:
            fh.write('\t'.join(build_headers(result_file)) + '\n')
    # Simple console summary
    if sv0 is None and sv1 is None:
        warning(f'Cannot find gene alignments for {isolate.name} in all databases.\n')
    elif sv0 is None and sv1 is not None and sv1.confidence != 'Untypeable':
        print(f'{isolate.name}: {species.reference}, no good match in main database, best singleton match: {sv1.output_name(min_locus_id, min_locus_cov)}')
    else:
        sin_suffix = '' if (sv1 is None or sv1.confidence == 'Untypeable') else f' singleton: {sv1.output_name(min_locus_id, min_locus_cov)}'
        print(f"{isolate.name}: {species.reference} {sv0.output_name(min_locus_id, min_locus_cov)}{sin_suffix}")

    # Build and write TSV line
    line = [
        isolate.name, species.reference, species.hit_strength,
        str(species.distance), str(species.p_value), str(species.shared_hashes),
    ]
    for scheme in ['mlst_scheme_1', 'mlst_scheme_2', 'mlst_scheme_3']:
        line += list(result_file['MLST'][scheme].values())
    line += build_serovar_fields(sv0, sv1, min_locus_id, min_locus_cov)

    with open(output_path, 'at') as fh:
        fh.write('\t'.join(line) + '\n')

    # Write rfb locus FASTA sequence
    if not no_rfb_sequence:
        selected = select_rfb_serovar(sv0, sv1)
        if selected:
            write_rfb_sequence(selected, Path(output_path).parents[0])
        else:
            warning(f'Low sequence similarity for {isolate.name} in the database, no sequence output.')

    # If the identity and coverage are both below 50%, we consider the serovar prediction unreliable and do not output the cps sequence.
    # TODO: We may want to make a figure for the rfb locus in the input.
    '''
    if figure and result_file['Serovar'].best_match:
        from dna_features_viewer import GraphicFeature, GraphicRecord
        draw_gene_map(inputfile, serotype, orfs, sub_cps_length)
    '''