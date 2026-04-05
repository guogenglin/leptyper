# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 15:40:29 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from io import TextIOBase
from warnings import catch_warnings
from typing import Optional, TextIO
from itertools import groupby, chain
from functools import cached_property
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from collections.abc import Generator, Iterable, Callable

from .log import log, warning
from .alignment import Alignment
from .utils import merge_ranges, LeptyperError

class Database:
    '''
    This class represents the database of loci and genes for serotyping.
    '''
    def __init__(self, name: str, loci: Optional[dict[str, Locus]] = None, genes: Optional[dict[str, Gene]] = None):
        self.name = name
        self.loci = loci or {}
        self.genes = genes or {}
        self.gene_threshold = 0.8
        self._expected_gene_counts = None

    def __repr__(self):
        return f"{self.name} ({len(self.loci)} Loci) ({len(self.genes)} Genes)"

    def __str__(self) -> str:
        return self.name

    def __len__(self) -> int:
        return len(self.loci)

    def __iter__(self):
        return iter(self.loci.values())

    def __getitem__(self, item: str | int) -> Locus:
        if isinstance(item, int):
            if not 0 <= item < len(self):
                raise LeptyperError(f'Index {item} out of range for database {self.name}')
            return list(self.loci.values())[item]
        else:
            raise LeptyperError(f'Could not find {item} in database {self.name}')

    @property
    def expected_gene_counts(self) -> np.ndarray:
        if self._expected_gene_counts is None:
            self._expected_gene_counts = np.array([len(l.genes) for l in self.loci.values()])
        return self._expected_gene_counts

    def format(self, format_spec):
        if format_spec in {'fna'}:
            return ''.join([locus.format(format_spec) for locus in self.loci.values()])
        elif format_spec in {'ffn', 'faa'}:
            return ''.join([gene.format(format_spec) for gene in self.genes.values()])
        else:
            raise LeptyperError(f'Invalid format specifier: {format_spec}')

    def add_locus(self, locus: Locus):
        """
        Adds a locus and its genes to the database. Checks that the locus and genes don't already exist in the database.
        """
        if locus.name in self.loci:
            raise LeptyperError(f'Locus {locus.name} already exists in database {self.name}.')
        self.loci[locus.name] = locus
        for gene in locus:
            if gene.name in self.genes:
                raise LeptyperError(f'Gene {gene} already exists in database {self.name}.')
            self.genes[gene.name] = gene


class Gene:
    '''
    This class represents a gene in the database. 
    It contains the gene name, location, strand, and sequence information.
    '''

    def __init__(self, name: str, start: int = 0, end: int = 0, strand: str = "+",
                 protein_seq: Optional[Seq] = None, dna_seq: Optional[Seq] = None, 
                 gene_name: Optional[str] = None, product: Optional[str] = None):
        self.name = name or ''
        self.start = start  # 0-based
        self.end = end
        self.strand = strand  # Either + or -
        self.gene_name = gene_name or ''
        self.product = product or ''  # Can also be description
        self.dna_seq = dna_seq or Seq('')
        self.protein_seq = protein_seq or Seq('')

    def __hash__(self):
        return hash(self.name)  # The name of the gene is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.dna_seq)

    def format(self, format_spec):
        if format_spec == 'ffn':
            if len(self.dna_seq) == 0:
                warning(f'No DNA sequence for {self}')
                return ""
            return f'>{self.name}\n{self.dna_seq}\n'
        if format_spec == 'faa':
            self.extract_translation()
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for {self.__repr__()}')
                return ""
            return f'>{self.name}\n{self.protein_seq}\n'
        raise LeptyperError(f'Invalid format specifier: {format_spec}')

    def extract_translation(self, **kwargs):
        """
        Extracts the protein sequence from the DNA sequence of the gene. Implemented as a method so unnecessary
        translations are not performed.
        :param table: NCBI translation table number
        :param cds: if True, only translates the CDS
        :param to_stop: if True, stops translation at the first stop codon
        :param gap: gap character
        :param stop_symbol: stop codon character
        """
        if len(self.protein_seq) == 0:  # Only translate if the protein sequence is not already stored
            if len(self.dna_seq) == 0:
                raise LeptyperError(f'No DNA sequence for reference {self}')
            with catch_warnings(record=True) as w:
                self.protein_seq = self.dna_seq.translate(**kwargs)
                # for i in w:
                #     warning(f"{i.message}: {self.__repr__()}")
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for reference {self}')

class Locus:
    def __init__(self, name: Optional[str] = None, seq: Seq | None = Seq(''), genes: Optional[dict[str, Gene]] = None,
                 index: int | None = 0):
        self.name = name or ''
        self.seq = seq or Seq('')
        self._length = len(self.seq)
        self.genes = genes or {}
        self.index = index

    @classmethod
    def from_seqrecord(cls, record: SeqRecord, locus_name: str):
        self = cls(name=locus_name, seq=Seq(record.seq))
        n = 1
        for feature in record.features:
            if feature.type == 'CDS':

                name = f"{locus_name}_{str(n).zfill(4)}" + (
                    f"_{gene_name}" if (gene_name := feature.qualifiers.get('gene', [''])[0]) else '')

                gene = Gene(
                    name=name, gene_name=gene_name, dna_seq=feature.extract(record.seq), start=int(feature.location.start),
                    end=int(feature.location.end), strand='+' if feature.location.strand == 1 else '-',
                    product=feature.qualifiers.get('product', [''])[0]
                )

                self.genes[name] = gene
                n += 1
        return self

    def __hash__(self):  # TODO: Check if this is used at all
        return hash(self.name)  # The name of the locus is unique and immutable

    def __repr__(self):
        return self.name

    def __len__(self):
        return self._length or len(self.seq)

    def __getitem__(self, item) -> Gene:
        if not (result := self.genes.get(item)):
            raise LeptyperError(f'Could not find {item} in locus {self.name}')
        return result

    def __iter__(self):
        return iter(self.genes.values())

    def format(self, format_spec):
        if format_spec == 'fna':
            if len(self.seq) == 0:
                warning(f'No DNA sequence for {self}')
                return ""
            return f'>{self.name}\n{self.seq}\n'
        if format_spec in {'ffn', 'faa'}:
            return ''.join([gene.format(format_spec) for gene in self])
        raise LeptyperError(f'Invalid format specifier: {format_spec}')

    def write(self, fna: str | PathLike | TextIO = None, ffn: str | PathLike | TextIO = None,
              faa: str | PathLike | TextIO = None):
        """Write the typing result to files or file handles."""
        for f, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')]:
            if f:
                if isinstance(f, TextIOBase):
                    f.write(self.format(fmt))
                elif isinstance(f, PathLike) or isinstance(f, str):
                    with open(path.join(f, f'{self.name.replace("/", "_")}.{fmt}'), 'wt') as handle:
                        handle.write(self.format(fmt))

def name_from_record(record: SeqRecord) -> str | None:
    '''
    Extracts the locus name from the note qualifier of the source feature of a genbank record.
    '''
    locus_name = set()

    if not (source := next((f for f in record.features if f.type == 'source'), None)):
        raise LeptyperError(f'Could not find source feature in genbank record: {record.id}')
    if "note" not in source.qualifiers:
        raise LeptyperError(f'Could not find note qualifier in source feature of genbank record: {record.id}')

    for note in source.qualifiers['note']:
        locus_name.add(note.replace('Serovar:', '').strip())

    if len(locus_name) > 1:
        raise LeptyperError(f'Found multiple locus names in record: {record.id}\n\tNote: {source.qualifiers["note"]}')

    return locus_name.pop() if len(locus_name) == 1 else None

def parse_database(db_path: Path, verbose: bool = False) -> Generator[Locus, None, None]:
    '''
    This function parses the genbank file and yields Locus objects. 
    It uses the name_from_record function to extract the locus name from the record.
    '''
    log(f'Parsing {db_path.stem}', verbose=verbose)
    for record in SeqIO.parse(db_path, 'genbank'):
            locus_name = name_from_record(record)
            if not locus_name:
                raise LeptyperError(f'Could not parse locus name from {record.id}')
            yield Locus.from_seqrecord(record, locus_name)

def load_database(reference: str, verbose: bool = False) -> Database:
    '''
    Load the database from the genbank file and return a Database object. 
    The genbank file must be in the same directory as this script and named 'Leptospira_lps_locus_reference.gbk'.
    '''
    db_relevant_path = Path(reference)
    db_abs_path = Path(__file__).parents[0] / 'database' / reference
    if not db_abs_path.is_file():
        db_path = db_relevant_path
    else:
        db_path = db_abs_path
    if not db_path.is_file():
        raise LeptyperError(f"Database file not found: {db_path}")
    db = Database(db_path.stem)
    for locus in parse_database(db_path, verbose=verbose):
        db.add_locus(locus)
    if not db.loci:  # Check that loci were properly loaded
        raise LeptyperError(f'No loci found in database {db.name}')
    for n, locus in enumerate(db.loci.values()):
        locus.index = n
    ''' The reason we need to assign a .index to each locus is so that we will need to arrange the results 
    in a numpy array for scoring, and we need to know which row corresponds to which locus. This is done in 
    the lepto_serotyping function.
    '''
    return db

def group_alns(alignments: Iterable[Alignment], key: str = 'q') -> Generator[tuple[str, Generator[Alignment]]]:
    '''
    Group alignments by a key (query gene name). Yields tuples of (key, generator of alignments).
    '''
    for q, group in groupby(sorted(alignments, key=lambda x: getattr(x, key)), key=lambda x: getattr(x, key)):
        yield q, (item for item in group)

class TypingResult:
    """
    This is a class to store the results of a typing analysis for a single sample. It is designed to be flexible
    enough to store results from both read and assembly typing, and can be easily reconstructed from JSON to use
    with the `convert` utility. It should not store any information that is not directly related to the typing
    such as contig sequences or read alignments.
    """

    def __init__(
            self, sample_name: str | None, db: Database | None, best_match: Locus = None,
            pieces: list[LocusPiece] = None, expected_genes_inside_locus: list[GeneResult] = None,
            expected_genes_outside_locus: list[GeneResult] = None, missing_genes: list[str] = None,
            unexpected_genes_inside_locus: list[GeneResult] = None,
            unexpected_genes_outside_locus: list[GeneResult] = None,
            extra_genes: list[GeneResult] = None):
        self.sample_name = sample_name or ""
        self.db = db
        self.best_match = best_match
        self.pieces = pieces or []  # Pieces of locus reconstructed from alignments
        self.expected_genes_inside_locus = expected_genes_inside_locus or []  # Genes from best_match
        self.expected_genes_outside_locus = expected_genes_outside_locus or []  # Genes from best_match
        self.missing_genes = missing_genes or []  # Genes from best_match that were not found
        self.unexpected_genes_inside_locus = unexpected_genes_inside_locus or []  # Genes from other loci
        self.unexpected_genes_outside_locus = unexpected_genes_outside_locus or []  # Genes from other loci
        self.extra_genes = extra_genes or []  # in db.extra_genes, ALWAYS outside locus (gene_result.piece == None)
        # Properties to cache the values, these are protected to prevent accidental modification
        self._percent_identity = None
        self._percent_coverage = None
        self._phenotype = None
        self._problems = None
        self._confidence = None

    def __repr__(self):
        return f"{self.sample_name} {self.best_match.name}"

    def __len__(self):
        return sum(len(i) for i in self.pieces) if self.pieces else 0

    def __iter__(self):
        return chain(
            self.expected_genes_inside_locus, self.unexpected_genes_inside_locus,
            self.expected_genes_outside_locus, self.unexpected_genes_outside_locus, self.extra_genes)

    def add_gene_result(self, gene_result: GeneResult):
        if gene_result.piece:  # If gene_result.piece is not None, the gene is inside the locus
            gene_result.piece.add_gene_result(gene_result)
            gene_type = f"{gene_result.gene_type}{'_inside_locus' if gene_result.gene_type.startswith(('expected', 'unexpected')) else ''}"
        else:  # If gene_result.piece is None, the gene is outside the locus
            gene_type = f"{gene_result.gene_type}{'_outside_locus' if gene_result.gene_type.startswith(('expected', 'unexpected')) else ''}"
        getattr(self, gene_type).append(gene_result)  # Add gene result to the appropriate list

    @property
    def output_name(self) -> str:
        return f'{self.best_match.name}' if (self.percent_identity >=95 and self.percent_coverage >= 95) else f'{self.best_match.name}?'
    
    @property
    def percent_identity(self) -> float:
        if self._percent_identity is None:
            self._percent_identity = (sum(i.percent_identity for i in x) / len(x)) if (
                x := self.expected_genes_inside_locus) else 0
        return self._percent_identity

    @property
    def percent_coverage(self) -> float:
        if self._percent_coverage is None:
            self._percent_coverage = sum(len(i) for i in x) / sum(len(i) for i in self.best_match.genes.values()) * 100 \
                if (x := self.expected_genes_inside_locus) else 0
        return self._percent_coverage

    @property
    def phenotype(self) -> str:
        if self._phenotype is None:
            gene_phenotypes = set()  # Init set to store gene phenotypes to be used as a key in the phenotypes dict
            for gene in self:
                if gene.gene_type in {'expected_genes', 'extra_genes'}:  # The reported phenotype only considers these
                    gene_phenotypes.add((gene.gene.name, gene.phenotype))  # or extra genes
            # NOTE: The best_match.phenotypes MUST be sorted from largest to smallest gene set to make sure any sets with
            # extra genes are tested first.
            self._phenotype = next(
                (p for gs, p in self.best_match.phenotypes if len(gs) == len(gene_phenotypes.intersection(gs))),
                self.best_match.type_label)  # If no phenotype is found, return the type label
        return self._phenotype

    @property
    def problems(self) -> str:
        if self._problems is None:
            self._problems = f'?{x}' if (x := len(self.pieces)) != 1 else ''
            self._problems += '-' if self.missing_genes else ''
            self._problems += '+' if self.unexpected_genes_inside_locus else ''
            self._problems += '*' if any(
                i.percent_coverage >= 90 and i.below_threshold for i in self.expected_genes_inside_locus) else ''
            self._problems += '!' if any(i.phenotype == "truncated" for i in self) else ''
        return self._problems

    @property
    def confidence(self) -> str:
        return self._confidence if self._confidence is not None else "Not calculated"

    def get_confidence(self, max_other_genes: int = 1, percent_expected_genes: float = 50):
        p = len(set(i.gene.name for i in self.expected_genes_inside_locus)) / len(self.best_match.genes) * 100
        other_genes = len(
            set(i.gene.name for i in self.unexpected_genes_inside_locus if not i.phenotype == "truncated"))
        if "*" in self.problems:
            self._confidence = "Untypeable"
        else:
            if len(self.pieces) == 1 and not self.missing_genes and not other_genes:
                self._confidence = "Typeable"
            elif other_genes <= max_other_genes and p >= percent_expected_genes:
                self._confidence = "Typeable"
            else:
                self._confidence = "Untypeable"

    @classmethod
    def from_dict(cls, d: dict, db: Database) -> TypingResult:
        if not (best_match := db.loci.get(d['best_match'])):
            raise TypingResultError(f"Best match {d['best_match']} not found in database")
        self = TypingResult(sample_name=d['sample_name'], db=db, best_match=best_match,
                            missing_genes=d['missing_genes'])
        # Set the cached properties
        self._percent_identity = float(d['percent_identity'])
        self._percent_coverage = float(d['percent_coverage'])
        self._phenotype = d['phenotype']
        self._problems = d['problems']
        self._confidence = d['confidence']
        # Add the pieces and create the gene results
        self.pieces = [LocusPiece.from_dict(i, result=self) for i in d['pieces']]
        pieces = {i.__repr__(): i for i in self.pieces}
        gene_results = {}  # This was previously a dict comp, but we need to check the gene is in the database, see #31
        for r in chain(d['expected_genes_inside_locus'], d['unexpected_genes_inside_locus'],
                       d['expected_genes_outside_locus'], d['unexpected_genes_outside_locus'],
                       d['extra_genes']):
            if not (gene := db.genes.get(r['gene'])) and not (gene := db.extra_genes.get(r['gene'])):
                raise TypingResultError(f"Gene {r['gene']} not found in database")
            x = GeneResult.from_dict(r, result=self, piece=pieces.get(r['piece']), gene=gene)
            gene_results[x.__repr__()] = x

        for gene_result in gene_results.values():
            self.add_gene_result(gene_result)
        return self

    def format(self, format_spec) -> str | dict:
        if format_spec == 'tsv':
            return '\t'.join(
                [
                    self.sample_name,
                    self.best_match.name,
                    self.phenotype,
                    self.confidence,
                    self.problems,
                    f"{self.percent_identity:.2f}%",
                    f"{self.percent_coverage:.2f}%",
                    f"{self.__len__() - len(self.best_match)} bp" if len(self.pieces) == 1 else 'n/a',
                    f"{(n_inside := len({i.gene.name for i in self.expected_genes_inside_locus}))} / {(n_expected := len(self.best_match.genes))} ({100 * n_inside / n_expected:.2f}%)",
                    ';'.join(map(str, self.expected_genes_inside_locus)),
                    ';'.join(self.missing_genes),
                    f"{len(self.unexpected_genes_inside_locus)}",
                    ';'.join(map(str, self.unexpected_genes_inside_locus)),
                    f"{(n_outside := len({i.gene.name for i in self.expected_genes_outside_locus}))} / {n_expected} ({100 * n_outside / n_expected:.2f}%)",
                    ';'.join(map(str, self.expected_genes_outside_locus)),
                    f"{len(self.unexpected_genes_outside_locus)}",
                    ';'.join(map(str, self.unexpected_genes_outside_locus)),
                    ';'.join(map(str, filter(lambda z: z.phenotype == "truncated", self))),
                    ';'.join(map(str, self.extra_genes))
                ]
            ) + "\n"
        if format_spec == 'fna':  # Return the nucleotide sequence of the locus
            return "".join([i.format(format_spec) for i in self.pieces])
        if format_spec in {'faa', 'ffn'}:  # Return the protein or nucleotide sequence of gene results
            return "".join([i.format(format_spec) for i in self])
        if format_spec == 'json':
            return dumps(
                {
                    'sample_name': self.sample_name, 'best_match': self.best_match.name, 'confidence': self.confidence,
                    'phenotype': self.phenotype, 'problems': self.problems,
                    'percent_identity': str(self.percent_identity),
                    'percent_coverage': str(self.percent_coverage), 'missing_genes': self.missing_genes
                } | {
                    attr: [i.format(format_spec) for i in getattr(self, attr)] for attr in {
                        'pieces', 'expected_genes_inside_locus', 'unexpected_genes_inside_locus',
                        'expected_genes_outside_locus', 'unexpected_genes_outside_locus', 'extra_genes'
                    }
                }) + "\n"
        raise LeptyperError(f"Unknown format specifier {format_spec}")

    def write(self,
              tsv: TextIO = None,
              json: TextIO = None,
              fna: str | PathLike | TextIO = None,
              ffn: str | PathLike | TextIO = None,
              faa: str | PathLike | TextIO = None,
              plot: str | PathLike = None,
              plot_fmt: str = 'png'):
        """Write the typing result to files or file handles."""
        [f.write(self.format(fmt)) for f, fmt in [(tsv, 'tsv'), (json, 'json')] if isinstance(f, TextIOBase)]
        for f, fmt in [(fna, 'fna'), (ffn, 'ffn'), (faa, 'faa')]:
            if f:
                if isinstance(f, TextIOBase):
                    f.write(self.format(fmt))
                elif isinstance(f, PathLike) or isinstance(f, str):
                    with open(path.join(f, f'{self.sample_name}_kaptive_results.{fmt}'), 'wt') as handle:
                        handle.write(self.format(fmt))
        if plot:
            ax = self.format(plot_fmt).plot(figure_width=18)[0]  # type: 'matplotlib.axes.Axes'
            ax.set_title(f"{self.sample_name} {self.best_match} ({self.phenotype}) - {self.confidence}")
            ax.figure.savefig(path.join(plot, f'{self.sample_name}_kaptive_results.{plot_fmt}'), bbox_inches='tight')
            ax.figure.clear()  # TODO: Check if this is necessary

class LocusPiece:
    def __init__(self, id: str = None, result: TypingResult = None, start: int | None = 0,
                 end: int | None = 0, strand: str = None, sequence: Seq = None,
                 expected_genes: list[GeneResult] = None, unexpected_genes: list[GeneResult] = None,
                 extra_genes: list[GeneResult] = None):
        self.id = id or ''  # TODO: rename as seq_id for clarity, actual id is self.__repr__()
        self.result = result
        self.start = start
        self.end = end
        self.strand = strand or "unknown"
        self.sequence = sequence or Seq("")
        self.expected_genes = expected_genes or []  # Genes from best_match
        self.unexpected_genes = unexpected_genes or []  # Genes that were found from other loci
        self.extra_genes = extra_genes or []  # Genes that were found outside the locus

    def __len__(self):
        return self.end - self.start

    def __iter__(self):
        return chain(self.expected_genes, self.unexpected_genes, self.extra_genes)

    def __str__(self):
        return self.id

    def __repr__(self):
        return f"{self.id}:{self.start}-{self.end}{self.strand}"

    @classmethod
    def from_dict(cls, d: dict, **kwargs) -> LocusPiece:
        return cls(id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'],
                   sequence=Seq(d['sequence']), **kwargs)

    def format(self, format_spec) -> str | dict:
        if format_spec == 'fna':
            if self.strand == '-':
                sequence = self.sequence.reverse_complement()
                return f">{self.result.sample_name}|{self.id}:{self.start + 1}-{self.end}_{self.strand}_rev\n{sequence}\n"
            else:
                sequence = self.sequence
                return f">{self.result.sample_name}|{self.id}:{self.start + 1}-{self.end}_{self.strand}\n{sequence}\n"
        if format_spec == 'json':
            return {'id': self.id, 'start': str(self.start), 'end': str(self.end), 'strand': self.strand,
                    'sequence': str(self.sequence)}
        raise ValueError(f"Unknown format specifier {format_spec}")

    def add_gene_result(self, gene_result: GeneResult):
        '''
        Pending the boundary of the rfb locus
        '''
        if gene_result.start < self.start:  # Update start and end if necessary
            self.start = gene_result.start
        if gene_result.end > self.end:
            self.end = gene_result.end
        getattr(self, gene_result.gene_type).append(gene_result)

def range_overlap(range1: tuple[int, int], range2: tuple[int, int], skip_sort: bool = False) -> int:
    """
    Returns the overlap between two ranges
    :param range1: Tuple of start and end positions
    :param range2: Tuple of start and end positions
    :param skip_sort: Skip sorting each range before calculating the overlap
    :return: Integer of overlap
    """
    start1, end1 = range1 if skip_sort else sorted(range1)
    start2, end2 = range2 if skip_sort else sorted(range2)
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    return max(0, overlap_end - overlap_start)

def cull(keep: Alignment, alignments: Iterable[Alignment],
         overlap_fraction: float = 0.1) -> Generator[Alignment]:
    '''
    Cull alignments that overlap with the keep alignment by more than the specified fraction of the alignment length.
    Setting of overlap_fraction because there are some genes in bacteria genome that have natural overlap, but they
    are different, this Cull process is designed to only remove alignments that are likely to be the same gene, 
    overlap lower than 0.1 blen means only ~100 bp overlap, and the length of a gene usually is around 1k bp.
    '''
    for a in alignments:
        if (a.ctg != keep.ctg or  # Different contig
                range_overlap((a.r_st, a.r_en), (keep.r_st, keep.r_en), skip_sort=True) / a.blen < overlap_fraction):
            yield a

def cull_all(alignments: list[Alignment]) -> list[Alignment]:
    kept_alignments = []
    sorted_alignments = sorted(list(alignments), key=lambda x: x.mlen, reverse=True)
    # sort based on mlen, from largest to smallest, so that the best alignment will be kept first
    while sorted_alignments:
        kept_alignments.append(sorted_alignments.pop(0))
        sorted_alignments = list(cull(kept_alignments[-1], sorted_alignments))
    return kept_alignments

def cull_filtered(pred: Callable, alignments: Iterable[Alignment]) -> Generator[Alignment]:
    """
    Cull and flatten alignments that don't overlap with alignments matching the predicate.
    E.g. list(cull_filtered(lambda i: i.q in query_genes, alignments))
    """
    keep, other = [], []
    [(keep if pred(a) else other).append(a) for a in alignments]
    # Only keep the alignments that match the gene in the best match locus
    # And cull all the other alignments.
    other = cull_all(other)  # Remove conflicting other alignments
    for i in keep:  # Remove other alignments overlapping best match gene alignments
        other = list(cull(i, other))
        yield i
    yield from other

class GeneResult:
    """
    Class to store alignment results for a single gene in a locus for either a ReadResult or a AssemblyResult.
    """

    def __init__(self, id: str, gene: Gene, result: TypingResult = None,
                 piece: LocusPiece | None = None, start: int | None = 0, end: int | None = 0, strand: str | None = None,
                 dna_seq: Seq | None = Seq(""), protein_seq: Seq | None = Seq(""), below_threshold: bool | None = False,
                 phenotype: str | None = "present", gene_type: str | None = None, partial: bool | None = False,
                 percent_identity: float | None = 0, percent_coverage: float | None = 0):
        self.id = id or ''  # TODO: rename as seq_id for clarity, actual id is self.__repr__()
        self.gene = gene
        self.result = result  # TODO: replace with sample_name (only TypingResult attribute used, needed for fa headers)
        self.start = start
        self.end = end
        self.strand = strand
        self.partial = partial
        self.piece = piece  # inside locus if not None
        self.dna_seq = dna_seq
        self.protein_seq = protein_seq
        self.below_threshold = below_threshold
        self.phenotype = phenotype
        self.gene_type = gene_type or ""
        self.percent_identity = percent_identity
        self.percent_coverage = percent_coverage

    def __repr__(self):
        return f"{self.gene.name} {self.id}:{self.start}-{self.end}{self.strand}"

    def __len__(self):
        return self.end - self.start

    def __str__(self) -> str:
        s = f'{self.gene.name},{self.percent_identity:.2f}%,{self.percent_coverage:.2f}%'
        s += ",partial" if self.partial else ""
        s += ',truncated' if self.phenotype == "truncated" else ""
        s += ",below_id_threshold" if self.below_threshold else ""
        return s

    @classmethod
    def from_dict(cls, d: dict, **kwargs) -> GeneResult:
        return cls(
            id=d['id'], start=int(d['start']), end=int(d['end']), strand=d['strand'], dna_seq=Seq(d['dna_seq']),
            protein_seq=Seq(d['protein_seq']), below_threshold=True if d['below_threshold'] == 'True' else False,
            phenotype=d['phenotype'], gene_type=d['gene_type'], partial=True if d['partial'] == 'True' else False,
            percent_identity=float(d['percent_identity']), percent_coverage=float(d['percent_coverage']), **kwargs
        )

    def format(self, format_spec) -> str | dict:
        if format_spec == 'ffn':
            if len(self.dna_seq) == 0:
                warning(f'No DNA sequence for {self}')
                return ""
            return (f'>{self.gene.name} {self.result.sample_name}|{self.id}:{self.start}-{self.end}{self.strand}\n'
                    f'{self.dna_seq}\n')
        if format_spec == 'faa':
            if len(self.protein_seq) == 0:
                warning(f'No protein sequence for {self.__repr__()}')
                return ""
            return (f'>{self.gene.name} {self.result.sample_name}|{self.id}:{self.start}-{self.end}{self.strand}\n'
                    f'{self.protein_seq}\n')
        if format_spec == 'json':
            return {
                'id': self.id, 'start': str(self.start), 'end': str(self.end), 'strand': self.strand,
                'dna_seq': str(self.dna_seq), 'protein_seq': str(self.protein_seq), 'partial': str(self.partial),
                'below_threshold': str(self.below_threshold), 'phenotype': self.phenotype, 'gene_type': self.gene_type,
                'percent_identity': str(self.percent_identity), 'percent_coverage': str(self.percent_coverage),
                'gene': self.gene.name, 'piece': self.piece.__repr__() if self.piece else '',
            }
        raise LeptyperError(f"Unknown format specifier {format_spec}")

    def compare_translation(self, truncation_tolerance: float = 95, **kwargs):
        """
        Extracts the translation from the DNA sequence of the gene result.
        Will also extract the translation from the gene if it is not already stored.
        """
        # One translation for gene, one for gene alignment result.
        self.gene.extract_translation(**kwargs)  # Extract the translation from the gene if it is not already stored
        if len(self.dna_seq) == 0:  # If the DNA sequence is empty, raise an error
            raise LeptyperError(f'No DNA sequence for {self.__repr__()}')
        with catch_warnings(record=True) as w:  # Catch Biopython warnings
            protein_seqs = [self.dna_seq[i:].translate(**kwargs) for i in range(3)]  # Translate in all 3 frames
        frame, self.protein_seq = max(enumerate(protein_seqs), key=lambda x: len(x[1]))  # Get the longest translation
        self.start += frame  # Update the start position to the frame with the longest translation
        if len(self.protein_seq) <= 1:  # If the protein sequence is still empty, raise a warning
            warning(f'No protein sequence for {self.__repr__()}')
        elif len(self.gene.protein_seq) > 1:  # If both sequences are not empty
            if alignments := PairwiseAligner(scoring='blastp', mode='local').align(self.gene.protein_seq, self.protein_seq):  # Align the sequences
                alignment = max(alignments, key=lambda x: x.score)  # Get the best alignment
                self.percent_identity = alignment.counts().identities / alignment.shape[1] * 100
                self.percent_coverage = (len(self.protein_seq) / len(self.gene.protein_seq)) * 100
                if (self.percent_coverage < truncation_tolerance and  # If the coverage is less than the tolerance
                        # not partial, and not unexpected gene outside locus
                        not self.partial and not (self.gene_type == "unexpected_genes" and not self.piece)):
                    self.phenotype = "truncated"  # Set the phenotype to truncated
            else:
                warning(f'Error aligning {self.__repr__()}')

def lepto_serotyping(assembly_obj, db: Database, sin_db: Database | None, threads: int = 0, min_cov: float = 50, 
                     n_best: int = 2, percent_expected_genes: float = 50, verbose: bool = False) -> TypingResult | None:

    # ALIGN GENES ------------------------------------------------------------------------------------------------------
    # Init scores array with 6 columns: AS, mlen, blen, q_len, genes_found, genes_expected
    scores, alignments = np.zeros((len(db), 6)), []
    # Group alignments by query gene (Alignment.q)
    for q, alns in group_alns(assembly_obj.map(db.format('ffn'), threads, verbose=verbose)):
        # eg. q, alns = "Copenhagen_0001", <generator object group_alns.<locals>.group_generator at 0x7f8c9c9d8e50>
        alignments.extend(alns := list(alns))  # Add all alignments to the list, and convert generator to list
        # Use the best alignment for each gene for scoring, if the coverage is above the minimum
        if ((best := max(alns, key=lambda x: x.mlen)).blen / best.q_len) * 100 >= min_cov:
            scores[db.loci[q.split('_', 1)[0]].index] += [
                best.tags['AS'], best.mlen, best.blen, best.q_len, 1, 0]
            # Sum the alignment score, mlen, blen, q_len for the best alignment of the gene to the locus.
        # For each gene, add: AS, mlen, blen, q_len, genes_found (1), genes_expected (0 but will update later)
    if scores.max() == 0:  # TODO：This gene alignment based method may not suitable for Leptospira, also, we can detect the flanking genes further instead.
        return warning(f'No gene alignments sufficient for typing {assembly_obj.name}\n'
                       f'Have you used the appropriate database for your species?')
    
    scores[:, 5] = db.expected_gene_counts  # Add expected genes (which is the total gene number of every loci)
    
    scores = scores[:, 0] * (scores[:, 4] / scores[:, 5])  # AS * (genes_found / genes_expected)

    best_loci = [db[int(i)] for i in
                 np.argsort(scores)[::-1][:min(n_best, len(scores))]]  # Get the best loci to fully align

    # ALIGN LOCUS ------------------------------------------------------------------------------------------------------
    scores, idx = np.zeros((len(best_loci), 4)), {l.name: i for i, l in enumerate(best_loci)}  # Init scores and index
    locus_alignments = {l.name: [] for l in best_loci}  # Init dict to store alignments for each locus
    # Group alignments by locus
    # Since last time we only looked at the best alignment for each gene, now we will align the full-length sequences of the loci
    best_match = None
    for locus, alns in group_alns(assembly_obj.map(''.join(i.format('fna') for i in best_loci), threads, verbose=verbose)):
        for a in alns:  # For each alignment of the locus
            # Sometimes other locus will have multiple irrelevant matches across the genome, which will affect the result,
            # But if there is a very good match (coverage and identity above 95%), we can be confident that this is the best match, 
            # and we can skip the rest of the alignments, which will save time and avoid noise from other matches.
            if a.coverage >= 95 and a.identity >= 95:
                best_match = best_loci[idx[locus]]
            scores[idx[locus]] += [a.tags['AS'], a.mlen, a.blen, a.q_len]  # Add alignment metrics to the scores
            # All match will be included, there maybe some noise, but this could avoid the rfb locus be separated into multiple pieces in different contigs.
            locus_alignments[locus].append(a)  # Add the alignment to the locus alignments

    if best_match is None:
        best_match = best_loci[np.argmax(scores[:, 0])]  # Get the best match based on the highest score
    
    result = TypingResult(assembly_obj.name, db, best_match)  # Create the result object

    # RECONSTRUCT LOCUS ------------------------------------------------------------------------------------------------
    # The original tool used largest locus for merge_ranges tolerance, but this may be too large for 
    # some loci, so we will use the actual length of the best match locus as the tolerance, which is more 
    # reasonable.
    pieces = {  # Init dict to store pieces for each contig
        ctg: [LocusPiece(ctg, result, s, e) for s, e in  # Create pieces for each merged contig range
              merge_ranges([(a.r_st, a.r_en) for a in alns], best_match._length)]
        for ctg, alns in group_alns(locus_alignments[best_match.name], key='ctg')  # Group by contig
    }  # We can't add strand as the pieces may be merged from multiple alignments, we will determine from the genes

    # GET GENE RESULTS -------------------------------------------------------------------------------------------------
    for a in cull_filtered(lambda i: i.q in best_match.genes, alignments):  # For each non-overlapping gene alignment
        if gene := best_match.genes.get(a.q):  # Get gene reference from database and gene type
            gene_type = "expected_genes"
        else:
            gene = db.genes.get(a.q)
            gene_type = "unexpected_genes"
        # Get Piece if gene range overlaps with a piece, return a LocusPiece object or None
        piece = next(filter(lambda p: range_overlap((p.start, p.end), (a.r_st, a.r_en)) > 0,
                            pieces.get(a.ctg, [])), None)
        # pieces.get(a.ctg, []) will get the pieces for the contig of the alignment
        # Create gene result and extract sequence from assembly
        gene_result = GeneResult(a.ctg, gene, result, piece, a.r_st, a.r_en, a.strand, gene_type=gene_type,
                                 partial=a.partial, dna_seq=assembly_obj.seq(a.ctg, a.r_st, a.r_en, a.strand))
        # Evaluate the gene in protein space by comparing the translation to the reference gene
        gene_result.compare_translation(table=11, to_stop=True)  # This will also trigger the protein alignment
        gene_result.below_threshold = gene_result.percent_identity < db.gene_threshold  # Check if below threshold
        if not piece and gene_result.below_threshold:  # If below protein identity threshold
            continue  # Skip this gene, probably a homologue in another part of the genome
        result.add_gene_result(gene_result)  # Add the gene result to the result to get neighbouring genes
        # previous_result = gene_result  # Set the previous gene result to the current result

    # FINALISE PIECES --------------------------------------------------------------------------------------------------
    for ctg, pieces in pieces.items():  # Add sequences to pieces and add them to the result
        for piece in pieces:
            if piece.expected_genes:  # If the piece has expected genes
                piece.strand = "+" if any(i.strand == i.gene.strand for i in piece.expected_genes) else "-"
                # i.stand is the strand from alignment, and i.gene.strand is the strand from the database
                # Piece strand is consensus of expected gene strands
                piece.sequence = assembly_obj.seq(ctg, piece.start, piece.end, piece.strand)
                result.pieces.append(piece)  # Add the piece to the result

    # FINALISE RESULT -------------------------------------------------------------------------------------------------
    # Sort the pieces by the sum of the expected gene order to get the expected order of the pieces
    result.pieces.sort(key=lambda x: min(i.gene.start for i in x.expected_genes))
    [l.sort(key=lambda x: gene.start) for l in (
        result.expected_genes_inside_locus, result.expected_genes_outside_locus, result.unexpected_genes_inside_locus,
        result.unexpected_genes_outside_locus)]
    result.missing_genes = list(set(best_match.genes) - {
        i.gene.name for i in chain(result.expected_genes_inside_locus, result.expected_genes_outside_locus)
    })
    result.get_confidence(1, percent_expected_genes)
    log(f"Finished typing {result}", verbose=verbose)
    return result