#!/usr/bin/env python3
"""Replay selected GBRV recipes with the reference-wave export patch."""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
import concurrent.futures
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
import threading


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--retry-failed', action='store_true')
    args = parser.parse_args()
    root = args.root
    catalog = json.loads((root / 'catalog.json').read_text())
    records = [r for r in catalog['entries'] if r['source_library'] == 'GBRV' and r.get('sssp', True)]
    source = root / 'sources/gbrv-recipes/Work_GBRV'
    output = root / 'work/gbrv'
    output.mkdir(exist_ok=True, parents=True)
    report = {'completed': False, 'entries': []}
    if args.retry_failed and (output / 'result.json').exists():
        previous = json.loads((output / 'result.json').read_text())['entries']
        report['entries'] = [r for r in previous if r['status'] == 'generated']
        retry = {r['original_sha256'] for r in previous if r['status'] != 'generated'}
        records = [r for r in records if r['original_sha256'] in retry]
    lock = threading.Lock()

    def run_case(row):
        name = Path(row['original_path']).name
        functional = 'pbesol' if '.pbesol.' in name else 'pbe'
        version = re.search(r'\.gbrv\.(v[0-9.]+)\.upf$', name).group(1)
        element = row['element']
        recipe_name = element.lower() + '_pbe_' + version
        recipe_dirs = list(source.glob('*/' + recipe_name))
        record = dict(row, functional=functional)
        try:
            if len(recipe_dirs) != 1:
                raise ValueError('No unique original recipe ' + recipe_name)
            recipe = recipe_dirs[0]
            case = output / (element + '-' + functional + '-' + version)
            case.mkdir(exist_ok=True)
            for path in recipe.glob('*.adat'):
                shutil.copyfile(path, case / path.name)
            ae = list(case.glob('*_ae_ref.adat'))
            if not ae:
                makefile = (recipe / 'Makefile').read_text()
                config = re.search(r'^CFG\s*=\s*(\S+)', makefile, re.M)
                if config:
                    ae = list(case.glob('*_ae_' + config.group(1) + '.adat'))
            ps = list(case.glob('*_ps.adat'))
            if len(ae) != 1 or len(ps) != 1:
                raise ValueError('Ambiguous AE or PS recipe')
            binary = root / 'sources' / ('vanderbilt-' + functional) / 'runatom.x'
            for input_file, stage in [(ae[0], 'ae'), (ps[0], 'ps')]:
                arguments = [str(binary), input_file.name, stage + '.out', 'atom.ae', stage + '.atwf', stage + '.logd', stage + '.uspp']
                with (case / (stage + '.log')).open('w') as out:
                    code = subprocess.run(arguments, cwd=case, stdout=out, stderr=subprocess.STDOUT, timeout=120).returncode
                if code:
                    if stage == 'ps' and (case / 'reference-waves.dat').exists():
                        record['generator_stopped_after_reference_export'] = True
                        record['generator_exit_status'] = code
                        lines = (case / 'ps.out').read_text(errors='replace').splitlines()
                        record['generator_diagnostic'] = next((line.strip() for line in reversed(lines) if 'error:' in line.lower()), lines[-1].strip())
                    else:
                        raise ValueError(stage + ' generator failed with status ' + str(code))
            export = case / 'reference-waves.dat'
            if not export.exists():
                raise ValueError('Missing complete reference export')
            record.update(status='generated', directory=str(case), export=str(export), export_sha256=sha256(export.read_bytes()).hexdigest(), generator_sha256=sha256(binary.read_bytes()).hexdigest(), recipe_directory=str(recipe.relative_to(source)), ae_input=str(ae[0]), ps_input=str(ps[0]))
        except (ValueError, subprocess.TimeoutExpired) as exc:
            record.update(status='failed', failure=str(exc))
        with lock:
            report['entries'].append(record)
            (output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
            print(element, functional, record['status'], record.get('failure',''), flush=True)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for _ in pool.map(run_case, records):
            pass
    report['completed'] = True
    (output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
