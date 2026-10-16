#!/usr/bin/env python3
"""Append generated PSLibrary potentials with classified full reference pairs."""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from prepare_upf import prepare_qe
from upf_data import values


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True,
                        help='Completed generate_pslibrary.py result.json')
    parser.add_argument('--collection', type=Path, required=True)
    parser.add_argument('--skip-failed', action='store_true',
                        help='Record unavailable generator cases without publishing a UPF')
    args = parser.parse_args()
    generation = json.loads(args.manifest.read_text())
    if not generation['completed']:
        raise ValueError('Reference generation is still running')
    failures = [r for r in generation['entries'] if r['status'] != 'generated']
    if failures and not args.skip_failed:
        raise ValueError('Incomplete reference exports: ' + repr(failures))
    catalog = json.loads((args.collection / 'catalog.json').read_text())
    existing = {r['file']: r for r in catalog['entries']}
    for source in sorted(generation['entries'], key=lambda r: (r['functional'], r['element'])):
        if source['status'] != 'generated':
            continue
        original, reference = Path(source['original']), Path(source['generated'])
        for key, path in [('original', original), ('generated', reference),
                          ('input', Path(source['input']))]:
            if digest(path) != source[key + '_sha256']:
                raise ValueError('Changed generator output: ' + str(path))
        raw = original.read_bytes()
        tree = ET.fromstring(raw)
        header = tree.find('PP_HEADER')
        element, functional = header.get('element').strip(), header.get('functional').strip().lower()
        if (element, functional) != (source['element'], source['functional']):
            raise ValueError('Generated element/functional differs from requested recipe')
        q = float(header.get('z_valence'))
        data, details = prepare_qe(raw, reference.read_bytes())
        variant = 'high' if source['recipe_job'] == 'us_ps_high.job' else 'low'
        filename = f'{element}.us.{functional}.z_{q:g}.ld1.psl.v1.0.0-{variant}.upf'
        relative = Path('data') / functional / 'additional/pslibrary' / filename
        target = args.collection / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.read_bytes() != data:
            raise ValueError('Refusing to replace different prepared data: ' + str(target))
        target.write_bytes(data)
        record = {k: source[k] for k in ('element', 'functional', 'source_functional',
                                         'source_library', 'sssp', 'recipe_job',
                                         'recipe_name', 'recipe_sha256', 'original_sha256')}
        record.update(details, status='prepared', file=str(relative),
                      sha256=sha256(data).hexdigest(), bytes=len(data), z_valence=q,
                      generator='Quantum ESPRESSO atomic 6.3MaX',
                      generator_sha256=generation['generator_sha256'],
                      generator_input='Embedded PP_INFO/PP_INPUTFILE with lsave_wfc enabled',
                      generator_input_sha256=source['input_sha256'],
                      reference_source_sha256=digest(reference),
                      pslibrary_commit=generation['source_commit'],
                      nlcc=any(values(tree.find('PP_NLCC'))) if tree.find('PP_NLCC') is not None else False,
                      model_core_tau_added=False, complete_projector_pairs=True,
                      original_operator_sections_unchanged=True,
                      cp2k_native_data_check_passed=False,
                      full_cp2k_calculation_qualified=False)
        old = existing.get(str(relative))
        if old is not None:
            if old['sha256'] != record['sha256']:
                raise ValueError('Catalog records different prepared data')
            # Reimporting verified data must not discard its native qualification.
            record['cp2k_native_data_check_passed'] = old['cp2k_native_data_check_passed']
        existing[str(relative)] = record
        print(element, functional, 'prepared', details['route'], flush=True)
    catalog['entries'] = list(existing.values())
    catalog['summary'] = dict(
        datasets=len(existing), prepared=len(existing), failed=0,
        elements=len({r['element'] for r in existing.values()}),
        routes=dict(Counter(str(r['route']) for r in existing.values())),
        upf_bytes=sum(r['bytes'] for r in existing.values()))
    catalog['status'] = 'prepared_pending_native_checks'
    unavailable = {(r['element'], r['functional']): r
                   for r in catalog.get('pslibrary_unavailable', [])}
    for source in generation['entries']:
        key = (source['element'], source['functional'])
        if source['status'] == 'generated':
            unavailable.pop(key, None)
        else:
            unavailable[key] = {k: source[k] for k in (
                'element', 'functional', 'source_functional', 'recipe_job',
                'recipe_name', 'recipe_sha256', 'failure')}
    catalog['pslibrary_unavailable'] = list(unavailable.values())
    (args.collection / 'catalog.json').write_text(json.dumps(catalog, indent=2) + '\n')
    print(json.dumps(catalog['summary']), flush=True)


if __name__ == '__main__':
    main()
