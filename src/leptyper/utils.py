# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:29:22 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

import os
import uuid
import gzip
import shutil
from typing import Any
from pathlib import Path
from operator import itemgetter
from collections.abc import Generator

from .log import log

class LeptyperError(RuntimeError):
    pass

def check_cpus(cpus: Any = None, max_cpus: int = 32, verbose: bool = False) -> int:
    '''
    Check the number of CPUs to use for parallel processing. If cpus is None, use all available.
    '''
    avail_cpus = os.cpu_count() or max_cpus
    if isinstance(cpus, str):
        cpus = int(cpus) if cpus.isdigit() else avail_cpus
    elif isinstance(cpus, float):
        cpus = int(cpus)
    else:
        cpus = avail_cpus
    cpus = min(cpus, avail_cpus, max_cpus)
    log(f'Using {cpus=}', verbose)
    return cpus

def check_programs_shutil(progs: list[str], verbose: bool = False):
    '''
    Check if programs are installed and executable using shutil.which.
    '''
    missing = []
    for prog in progs:
        path = shutil.which(prog)
        if path:
            print(f'{prog}: {path}')
        else:
            missing.append(prog)
    if not missing:
        log(f'All required programs found: {", ".join(progs)}', verbose = verbose)
    if missing:
        raise LeptyperError(f'Missing programs: {", ".join(missing)}')
    
def check_assembly(assembly: str | Path) -> Path | None:
    '''
    Check that an assembly file exists and is valid.
    '''
    assembly_path = Path(assembly)

    if not assembly_path.exists():
        raise LeptyperError(f'{assembly_path} does not exist')

    if not assembly_path.is_file() or assembly_path.stat().st_size == 0:
        raise LeptyperError(f'{assembly_path} is not a valid non-empty file')

    return assembly_path.resolve()

def get_compression_type(filename: Path) -> str:
    '''
    Determine the compression type of a file based on its magic bytes. Supports gzip, bzip2, and zip.
    '''
    filename = Path(filename)
    if not filename.is_file():
        raise LeptyperError(f"{filename} does not exist or is not a file")
    magic_dict = {
        'gz': b'\x1f\x8b\x08',
        'bz2': b'BZh',
        'zip': b'PK\x03\x04',
    }

    max_len = max(len(m) for m in magic_dict.values())
    with filename.open('rb') as f:
        file_start = f.read(max_len)

    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            if file_type in ['bz2', 'zip']:
                raise LeptyperError(f"{file_type} format is not supported; use gzip instead")
            return file_type

    return 'plain'

def gunzip_assembly(assembly: Path, temp_dir: str, verbose: bool = False) -> Path:
    '''
    If the assembly is gzipped, gunzip it to a temporary file and yield the path. 
    Otherwise, yield the original path.
    '''
    temp_path = Path(temp_dir)
    compression_type = get_compression_type(assembly)

    if compression_type == 'gz':
        unzipped_path = temp_path / f"{uuid.uuid4().hex}.fasta"
        log(f'Gunzip to temporary file: {unzipped_path}', verbose = verbose)

        with gzip.open(assembly, 'rb') as f_in, unzipped_path.open('wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        log('Decompression completed', verbose = verbose)
        return unzipped_path
    
    log(f'Assembly is not compressed: {assembly}', verbose = verbose)
    return assembly

def merge_ranges(ranges: list[tuple[int, int]], tolerance: int = 0, skip_sort: bool = False
                 ) -> Generator[tuple[int, int], None, None]:
    '''
    Concated those alignment pieces that are close.
    '''
    if not ranges:
        return None
    if len(ranges) == 1:
        yield ranges[0]
        return None
    current_range = (ranges := ranges if skip_sort else sorted(ranges, key=itemgetter(0)))[0]  # Start with the first range
    for start, end in ranges[1:]:  # Iterate through the ranges
        if start - tolerance <= current_range[1]:  # Overlap, merge the ranges
            current_range = (current_range[0], max(current_range[1], end))
        else:  # No overlap, add the current range to the merged list and start a new range
            yield current_range  # Yield the current range
            current_range = (start, end)   # Start a new range
    yield current_range  # Yield the last range