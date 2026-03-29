# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 16:46:12 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

import uuid
import subprocess
from Bio.Seq import Seq
from pathlib import Path
from typing import Generator
from subprocess import Popen, PIPE

from .log import log, warning
from .utils import LeptyperError

class Contig:
    """
    Represents a single contig in an assembly, with a name, description, and sequence.
    """
    def __init__(self, name: str, description: str, seq: Seq):
        self.name = name
        self.description = description
        self.seq = seq

    def __len__(self) -> int:
        return len(self.seq)

    def __repr__(self) -> str:
        return f"<Contig {self.name}, length={len(self)}>"

class Alignment:
    '''
    Represents a single alignment between a query and a contig, with coordinates, strand, and optional tags.
    '''
    def __init__(
            self, q: str | None = None, q_len: int | None = 0, q_st: int | None = 0,
            q_en: int | None = 0, strand: str | None = None, ctg: str | None = None,
            ctg_len: int | None = 0, r_st: int | None = 0, r_en: int | None = 0,
            mlen: int | None = 0, blen: int | None = 0, mapq: int | None = 0,
            tags: dict | None = None):
        self.q = q or 'unknown'  # Query sequence name
        self.q_len = q_len or 0  # Query sequence length, None -> 0
        self.q_st = q_st or 0  # Query start coordinate (0-based)
        self.q_en = q_en or 0  # Query end coordinate (0-based)
        self.strand = strand or 'unknown'  # ‘+’ if query/target on the same strand; ‘-’ if opposite
        self.ctg = ctg or 'unknown'  # Target sequence name
        self.ctg_len = ctg_len or 0  # Target sequence length, None -> 0
        self.r_st = r_st or 0  # Target start coordinate on the original strand (0-based)
        self.r_en = r_en or 0  # Target end coordinate on the original strand (0-based)
        self.mlen = mlen or 0  # Number of matching bases in the alignment
        self.blen = blen or 0  # Number bases, including gaps, in the alignment
        self.mapq = mapq or 0  # Mapping quality (0-255 with 255 for missing)
        self.tags = tags or {}  # {tag: value} pairs

    @classmethod
    def from_paf_line(cls, line: str):
        '''
        Parse a PAF line into an Alignment object. Expects at least 12 columns, with optional tags after.
        '''
        fields = line.rstrip().split('\t')

        if len(fields) < 12:
            raise LeptyperError(f'Invalid PAF line (<12 columns): {line}')
        
        try:
            def parse_tag(tag: str):
                key, type_, value = tag.split(':', 2)
                if type_ == 'i':
                    return key, int(value)
                if type_ == 'f':
                    return key, float(value)
                return key, value
            
            tags = dict(parse_tag(t) for t in fields[12:])
            
            return Alignment(q = fields[0], q_len = int(fields[1]), q_st = int(fields[2]), q_en = int(fields[3]), 
                             strand = fields[4], ctg = fields[5], ctg_len = int(fields[6]), r_st = int(fields[7]), 
                             r_en = int(fields[8]), mlen = int(fields[9]), blen = int(fields[10]), 
                             mapq = int(fields[11]), tags = tags)
        
        except Exception as e:
            raise LeptyperError(f'Error parsing PAF line: {line}') from e

    @property
    def identity(self) -> float:
        ''' Returns percent identity (0-100). Requires --eqx in minimap2. '''
        return 100.0 * self.mlen / self.blen if self.blen > 0 else 0.0

    @property
    def coverage(self) -> float:
        ''' Returns query coverage percentage (0-100). '''
        return 100.0 * (self.q_en - self.q_st) / self.q_len if self.q_len > 0 else 0.0
    
    def is_exact(self) -> bool:
        ''' Returns True if the alignment is a 100% match across the full query length. '''
        return self.mlen == self.blen and self.q_en - self.q_st == self.q_len

    @property
    def score(self) -> int:
        ''' Returns the alignment score (AS tag) if available. '''
        return self.tags.get('AS', 0)
    
    def __repr__(self) -> str:
        return f'{self.q}:{self.q_st}-{self.q_en} {self.ctg}:{self.r_st}-{self.r_en} {self.strand}'

    def __len__(self) -> int:
        return self.q_en - self.q_st

    def __getattr__(self, item):
        # First look in tags if not a normal attribute
        if item in self.tags:
            return self.tags[item]
        raise AttributeError(f'{self.__class__.__name__} object has no attribute {item}')

    @property
    def partial(self) -> bool:
        '''
        Determine if the alignment is partial.
        True if query does not cover entire contig or is clipped at ends.
        '''
        if self.q_len <= 0 or self.ctg_len <= 0:
            return False  # Cannot determine, treat as full

        # Alignment shorter than query or at contig edges
        aligned_len = self.q_en - self.q_st
        if aligned_len < self.q_len:
            if (self.r_st == 0 or self.r_en == self.ctg_len):
                return True

        # Check strand-specific overhangs
        query_start_diff = self.q_st if self.strand == '+' else self.q_len - self.q_en
        query_end_diff = self.q_len - self.q_en if self.strand == '+' else self.q_st
        ref_start_diff = self.r_st
        ref_end_diff = self.ctg_len - self.r_en

        if query_start_diff > ref_start_diff or query_end_diff > ref_end_diff:
            return True

        return False
    
class Assembly:
    '''
    Represents an assembly with multiple contigs. Provides methods to build minimap2 index and map queries.
    '''
    def __init__(self, path: Path, name: str | None = None, contigs: dict[str, Contig] | None = None):
        self.path = Path(path).resolve()
        self.name = name or self.path.stem
        self.contigs = contigs or {}
        self.mmi_index: Path | None = None

    def __repr__(self) -> str:
        return f"<Assembly {self.name}, {len(self.contigs)} contigs>"

    def __len__(self) -> int:
        return sum(len(c) for c in self.contigs.values())

    def seq(self, ctg: str, start: int, end: int, strand: str = '+') -> Seq:
        '''
        Retrieve the sequence for a given contig and coordinates, accounting for strand.
        '''
        if ctg not in self.contigs:
            raise LeptyperError(f"Contig {ctg} not found in assembly {self.name}")
        
        seq = self.contigs[ctg].seq[start:end]
        return seq if strand == '+' else seq.reverse_complement()

    def build_minimap2_index(self, temp_dir: str, verbose: bool = False) -> Path:
        if self.mmi_index and self.mmi_index.exists():
            log(f'Using existing minimap2 index: {self.mmi_index}', verbose = verbose)
            return self.mmi_index

        temp_path = Path(temp_dir)
        self.mmi_index = (temp_path / (uuid.uuid4().hex + '.mmi')).resolve()

        cmd = ['minimap2', '-d', str(self.mmi_index), str(self.path)]

        log(f'Building minimap2 index: {" ".join(cmd)}', verbose=verbose)

        p = subprocess.run(cmd, capture_output = True, text = True)

        if p.returncode != 0:
            log(p.stderr, verbose = True)
            raise LeptyperError(f'minimap2 failed indexing {self.name}:\n{p.stderr}')

        log(f'Index built successfully: {self.mmi_index}', verbose = verbose)

        return self.mmi_index
    
    def map(self, query: Path | str | None, threads: int = 1, extra_args: str = '', verbose: bool = False, 
            use_stdin: bool = True) -> Generator[Alignment, None, None]:
        '''
        Map a query sequence to the assembly using minimap2 and yield Alignment objects.
        '''
        if not self.mmi_index:
            raise LeptyperError(f'Index not built for assembly {self.name}')
        
        if use_stdin:
            cmd = f"minimap2 -c {extra_args} -t {threads} '{self.mmi_index}' -"
        else:
            cmd = f"minimap2 -c {extra_args} -t {threads} '{self.mmi_index}' '{query}'"

        log(f"Running minimap2: {cmd}", verbose=verbose)

        proc = Popen(cmd, shell = True, stdin = PIPE, stdout = PIPE, stderr = PIPE, text = True)
        
        stdout, stderr = proc.communicate(str(query) if use_stdin else None)

        if stderr:
            log(stderr.strip(), verbose = verbose)

        for line in stdout.splitlines():
            try:
                yield Alignment.from_paf_line(line)
            except LeptyperError as e:
                warning(f"Skipping malformed alignment line: {line}\n{e}")