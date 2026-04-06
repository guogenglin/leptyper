# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:11:07 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import tempfile
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .mash import lepto_species_mash
from .mlst import lepto_mlst
from .serotyping import lepto_serotyping, load_database
from .alignment import Contig, Assembly
from .utils import check_programs_shutil, check_cpus, check_assembly, gunzip_assembly, LeptyperError
from .log import log, formatted_description
from .output import generate_output


__version__ = '1.0.0'

def get_argument() -> argparse.ArgumentParser:
    '''
    Build the argument parser for the command-line interface.
    '''
    # Create the main parser
    parser = argparse.ArgumentParser(prog = 'leptyper', description = formatted_description, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Input and Output group
    io_group = parser.add_argument_group('Input and Output')
    io_group.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                          help = 'Input FASTA files')
    io_group.add_argument('-r', '--reference', type = str, default = 'Leptospira_lps_locus_reference.gbk', 
                          help = 'Reference rfb locus GenBank file')
    io_group.add_argument('-o', '--output', type = str, default = 'leptyper_output.txt', 
                          help = 'Output file')

    # Parameters group
    param_group = parser.add_argument_group('Parameters')
    param_group.add_argument('-t', '--threads', type = int, default = check_cpus(),
                             help = 'Number of alignment threads or 0 for all available CPUs')
    param_group.add_argument('--min-cov', type=float, default=50.0, 
                      help='Minimum gene %%coverage (blen/q_len*100) to be used for scoring (default: %(default)s)')
    param_group.add_argument('--n-best', type=int, default=2, 
                      help='Number of best loci from the 1st round of scoring to be\n'
                           'fully aligned to the assembly (default: %(default)s)')
    param_group.add_argument('--percent-expected', type=float, default=50.0,
                      help="Typeable if >= %% expected genes (default: %(default)s)")
    param_group.add_argument('--no_singleton',action = 'store_true',
                             help = 'Do not run alignment against the extra singleton database.')
    param_group.add_argument('--no_rfb_sequence', action = 'store_true',
                             help = 'Suppress output of rfb sequence file')
    param_group.add_argument('-f', '--figure', action = 'store_true',
                             help = 'Export the gene structure map of rfb locus')
    param_group.add_argument('-V', '--verbose', action = 'store_true',
                             help = 'Print debug messages to stderr')
    param_group.add_argument('--version', action = 'version', version = 'Leptyper v' + __version__,
                             help = 'Show version number and exit')

    return parser

def parse_input_assembly(assembly: Path, verbose: bool = False) -> Assembly:
    basename = assembly.stem
    log(f'Parsing input assembly: {assembly}', verbose = verbose)
    log(f'Assuming {basename} is in fasta format', verbose = verbose)
    # Parse fasta file
    contigs = {}

    try:
        with open(assembly, mode = 'rt') as f:
            for header, seq in SimpleFastaParser(f):
                parts = header.split(maxsplit = 1)
                name = parts[0]
                description = parts[1] if len(parts) > 1 else ''
                contigs[name] = Contig(name, description, Seq(seq))
        log(f'Parsed {len(contigs)} contigs from {basename}', verbose = verbose)

    except Exception as e:
        log(f'Failed parsing assembly: {assembly}', verbose = verbose)
        raise LeptyperError(f'Error parsing {assembly.name}: {e}') from e

    return Assembly(path = assembly, name = basename, contigs = contigs)
    

def main(argv: list[str] | None = None):
    '''
    Main function to run the Leptyper command-line interface.
    '''
    argv = sys.argv[1:] if argv is None else argv
    args = get_argument().parse_args(argv)
    # Check required external programs
    check_programs_shutil(['minimap2', 'mash'], verbose = args.verbose)
    db = load_database(args.reference, args.verbose)
    sin_db = None if args.no_singleton else load_database('singleton_database.gbk', args.verbose)
    for assembly in args.input:
        validated_assembly = check_assembly(assembly)
        if not validated_assembly:
            log(f'Skipping {assembly} due to validation failure.', args.verbose)
            print(f'Error: {assembly} is not a valid assembly file. Skipping.', file=sys.stderr)
            continue
        with tempfile.TemporaryDirectory() as temp_dir:
            # Prepare assembly
            assembly = gunzip_assembly(validated_assembly, temp_dir, args.verbose)
            assembly_obj = parse_input_assembly(assembly, args.verbose)
            assembly_obj.build_minimap2_index(temp_dir, args.verbose)

            # Collect results
            result_file = {
                'Isolate': assembly_obj,
                'Species': lepto_species_mash(assembly, args.verbose),
                'MLST': lepto_mlst(assembly_obj, args.verbose),
                'Serovar': lepto_serotyping(assembly_obj, db, sin_db, args.threads, args.min_cov, args.n_best, 
                                            args.percent_expected, args.no_singleton, args.verbose)
                                            }
            # Generate output
            generate_output(result_file, args.output, args.no_rfb_sequence, args.figure, args.verbose)

    
    
    
if __name__ == '__main__':
    raise SystemExit(main())
