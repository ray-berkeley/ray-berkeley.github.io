import csv
import hashlib
import io
import math

from nmrstarlib.nmrstarlib import NMRStarFile, StarFile


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def crosspeaks(payload: bytes, settings: dict) -> dict:
    require(0 < len(payload) <= 16 * 1024**2, 'Choose an NMR-STAR file up to 16 MB.')
    entry = StarFile.read(io.StringIO(payload.decode('utf-8-sig')), 'browser upload')
    require(isinstance(entry, NMRStarFile), 'Expected an NMR-STAR file with saveframes.')
    loops = [value for frame in entry.values() if isinstance(frame, dict)
             for name, value in frame.items() if name.startswith('loop_')
             and any(tag.startswith('Atom_chem_shift.') for tag in value[0])]
    require(bool(loops), 'No Atom_chem_shift loop found.')
    index = int(settings['loop'])
    require(1 <= index <= len(loops), f'Choose a shift list between 1 and {len(loops)}.')
    fields, records = loops[index - 1]
    tags = ['Entity_assembly_ID', 'Entity_ID', 'Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Val']
    tags = [f'Atom_chem_shift.{tag}' for tag in tags]
    require(len(fields) == len(set(fields)), 'Duplicate field names in the shift list.')
    require(all(tag in fields for tag in tags), 'Shift loop requires entity, residue, atom, and value fields.')
    direct, indirect = (settings[name].strip().upper() for name in ('direct', 'indirect'))
    require(bool(direct) and bool(indirect) and direct != indirect, 'Choose two different atom names.')
    residues = {}
    for record in records:
        assembly, entity, index, residue, atom, value = (record[tag] for tag in tags)
        if atom not in (direct, indirect) or value in ('.', '?'):
            continue
        key = (assembly, entity, index, residue)
        require(all(part not in ('.', '?', '') for part in key), 'Assignment identity is incomplete; resolve missing assembly, entity, or residue fields.')
        shifts = residues.setdefault(key, {})
        require(atom not in shifts, f'Duplicate assignment for {key}, {atom}; resolve ambiguity before exporting.')
        shifts[atom] = float(value)
        require(math.isfinite(shifts[atom]), 'Nonfinite chemical shift.')
    rows = [(*key, direct, values[direct], indirect, values[indirect])
            for key, values in residues.items() if direct in values and indirect in values]
    require(bool(rows), 'No residues contain both selected atoms in this shift list.')
    columns = ['Assembly', 'Entity', 'Residue index', 'Residue', 'Direct atom', 'Direct ppm', 'Indirect atom', 'Indirect ppm']
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(columns)
    writer.writerows(rows)
    return dict(columns=columns, rows=rows, count=len(rows), loops=len(loops), csv=output.getvalue(),
                source_sha256=hashlib.sha256(payload).hexdigest(), parser='nmrstarlib 2.1.1 (Python lexer)')
