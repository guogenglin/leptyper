# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:26:42 2026

@author: Genglin Guo
e-mail : gg599@drexel.edu
"""

from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from subprocess import Popen, PIPE

from .log import log
from .utils import LeptyperError


@dataclass(frozen=True)
class MashHit:
    reference: str
    distance: float
    p_value: float | None
    shared_hashes: str | None  # e.g. "123/1000"
    hit_strength: str | None = None  # e.g. "strong", "weak"


def lepto_species_mash(assembly: Path, verbose: bool = False) -> MashHit | None:
    '''
    Run: mash dist <sketch> <assembly>

    Mash output columns are typically:
      ref  query  dist  p-value  shared-hashes
    '''

    sketch = Path(__file__).parents[0] / 'database' / 'leptospira_ref.msh'
    log(f'Using mash sketch: {sketch}', verbose = verbose)

    if not sketch.is_file():
        raise LeptyperError(f'Mash sketch not found: {sketch}')

    cmd = ["mash", "dist", str(sketch), str(assembly)]
    log(f'Running command: {" ".join(cmd)}', verbose = verbose)

    proc = Popen(cmd, stdout = PIPE, stderr = PIPE, text = True)
    stdout, stderr = proc.communicate()

    if proc.returncode != 0:
        raise LeptyperError(f'mash dist failed (return code {proc.returncode}):\n{stderr.strip()}')

    log('mash dist completed successfully', verbose = verbose)

    best: MashHit | None = None
    best_dist = 1.0

    for line in stdout.splitlines():
        if not line.strip():
            continue
        
        cols = line.split('\t')
        if len(cols) < 3:
            continue

        species = cols[0].split('/')[0]
        try:
            dist = float(cols[2])
        except ValueError:
            continue

        # parse optional columns
        pval = float(cols[3]) if len(cols) >= 4 else None
        shared = cols[4] if len(cols) >= 5 else None

        if dist < best_dist:
            best_dist = dist

            # determine hit strength and reference label
            if dist <= 0.02:
                hit_strength = 'strong'
                reference_label = species
            elif dist <= 0.04:
                hit_strength = 'weak'
                reference_label = species
            else:
                hit_strength = ''
                reference_label = 'unknown'

            best = MashHit(
                reference = reference_label,
                distance = dist,
                p_value = pval,
                shared_hashes = shared,
                hit_strength = hit_strength
            )
        
    if best:
        log(
            f'Best hit: {best.reference} '
            f'(dist={best.distance:.4f}, strength={best.hit_strength})',
            verbose = verbose
        )
    else:
        log('No valid mash hits found', verbose = verbose)
        
    return best
