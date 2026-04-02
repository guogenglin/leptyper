# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:14:05 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

from __future__ import annotations

from pathlib import Path

def generate_output(result_file: dict, output_path: str, no_cps_sequence: bool = False, 
                    figure: bool = False, verbose: bool = False):
    if not Path(output_path).is_file():
        headers = ['Isolate', 'Species', 'Hit_strength', 'Distance', 'P-val', 'Shared_hashes', 
                   'mlst_scheme_1']
        headers += list(result_file['MLST']['mlst_scheme_1'].keys())[1:]
        headers += ['mlst_scheme_2']
        headers += list(result_file['MLST']['mlst_scheme_2'].keys())[1:]
        headers += ['mlst_scheme_3']
        headers += list(result_file['MLST']['mlst_scheme_3'].keys())[1:]
        headers += ['Serovar', 'Confidence', 'Problems', 'Identity', 'Coverage', 'Length discrepancy', 
                    'Expected genes in locus', 'Expected genes in locus, details', 'Missing expected genes', 
                    'Other genes in locus', 'Other genes in locus, details', 'Expected genes outside locus', 
                    'Expected genes outside locus, details', 'Other genes outside locus', 
                    'Other genes outside locus, details', 'Truncated genes, details', 'Extra genes, details']
        with open(output_path, 'wt') as file:
            file.write('\t'.join(headers))
            file.write('\n')


    simple_output = result_file['Isolate'].name + ': '  + result_file['Species'].reference + ' ' + result_file['Serovar'].output_name
    line = [result_file['Isolate'].name, result_file['Species'].reference, result_file['Species'].hit_strength, 
            str(result_file['Species'].distance), str(result_file['Species'].p_value), str(result_file['Species'].shared_hashes)]
    for scheme in ['mlst_scheme_1', 'mlst_scheme_2', 'mlst_scheme_3']:
        line += list(result_file['MLST'][scheme].values())
    line += [result_file['Serovar'].output_name, 
             result_file['Serovar'].confidence, 
             result_file['Serovar'].problems, 
             str(f"{result_file['Serovar'].percent_identity:.2f}" if result_file['Serovar'].percent_identity <= 100 else 100.00), 
             str(f"{result_file['Serovar'].percent_coverage:.2f}" if result_file['Serovar'].percent_coverage <= 100 else 100.00), 
             str(f"{result_file['Serovar'].__len__() - len(result_file['Serovar'].best_match)} bp" if len(result_file['Serovar'].pieces) == 1 else 'n/a'),
             str(f"{(n_inside := len({i.gene.name for i in result_file['Serovar'].expected_genes_inside_locus}))} / {(n_expected := len(result_file['Serovar'].best_match.genes))} ({100 * n_inside / n_expected:.2f}%)"),
            ';'.join(map(str, result_file['Serovar'].expected_genes_inside_locus)),
            ';'.join(result_file['Serovar'].missing_genes),
            str(f"{len(result_file['Serovar'].unexpected_genes_inside_locus)}"),
            ';'.join(map(str, result_file['Serovar'].unexpected_genes_inside_locus)),
            str(f"{(n_outside := len({i.gene.name for i in result_file['Serovar'].expected_genes_outside_locus}))} / {n_expected} ({100 * n_outside / n_expected:.2f}%)"),
            ';'.join(map(str, result_file['Serovar'].expected_genes_outside_locus)),
             str(f"{len(result_file['Serovar'].unexpected_genes_outside_locus)}"),
            ';'.join(map(str, result_file['Serovar'].unexpected_genes_outside_locus)),
            ';'.join(map(str, filter(lambda z: z.phenotype == "truncated", result_file['Serovar']))),
            ';'.join(map(str, result_file['Serovar'].extra_genes))
             ]

    print(simple_output)
    with open(output_path, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n')

    # If the identity and coverage are both below 50%, we consider the serovar prediction unreliable and do not output the cps sequence.
    # TODO: There has to be a rfb locus in the genome, we may construct a strain-specific rfb locus database for generate the ourput.
    if not no_cps_sequence and result_file['Serovar'].best_match and not (result_file['Serovar'].percent_coverage < 50 and result_file['Serovar'].percent_identity < 50):
        root_dir = Path(output_path).parents[0]
        if not Path(f'{root_dir}/rfb_locus').is_dir():
            Path(f'{root_dir}/rfb_locus').mkdir()
        with open(Path(f'{root_dir}/rfb_locus/{result_file['Serovar'].sample_name}_rfb_locus_sequence.fasta'), 'wt') as handle:
            handle.write("".join([i.format('fna') for i in result_file['Serovar'].pieces]))
    if figure and result_file['Serovar'].best_match:
        from dna_features_viewer import GraphicFeature, GraphicRecord
        draw_gene_map(inputfile, serotype, orfs, sub_cps_length)