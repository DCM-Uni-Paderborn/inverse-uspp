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
    parser.add_argument('--manifest', type=Path, help='Explicit GBRV dataset entries')
    parser.add_argument('--output', type=Path, help='Separate generator work directory')
    args = parser.parse_args()
    root = args.root
    catalog = json.loads((args.manifest or root / 'catalog.json').read_text())
    records = [r for r in catalog['entries'] if r['source_library'] == 'GBRV'
               and (args.manifest or r.get('sssp', True))]
    source = root / 'sources/gbrv-recipes/Work_GBRV'
    output = args.output or root / 'work/gbrv'
    output.mkdir(exist_ok=True, parents=True)
    report = {'completed': False, 'entries': []}
    if args.retry_failed and (output / 'result.json').exists():
        previous = json.loads((output / 'result.json').read_text())['entries']
        key = lambda r: r.get('case_name', r['original_sha256'])
        requested = {key(r): r for r in records}
        report['entries'] = [dict(r, **requested[key(r)])
                             for r in previous if r['status'] == 'generated'
                             and key(r) in requested
                             and sha256(Path(r['export']).read_bytes()).hexdigest() == r['export_sha256']]
        completed = {key(r) for r in report['entries']}
        records = [r for r in records if key(r) not in completed]
    lock = threading.Lock()

    def run_case(row):
        name = Path(row['original_path']).name
        functional = row.get('functional') or ('pbesol' if '.pbesol.' in name else 'pbe')
        version = row.get('version') or re.search(r'\.gbrv\.(v[0-9.]+)\.upf$', name).group(1)
        element = row['element']
        recipe_name = row.get('recipe_name') or element.lower() + '_pbe_' + version
        recipe_dirs = list(source.glob('*/' + recipe_name))
        record = dict(row, functional=functional)
        try:
            if len(recipe_dirs) != 1:
                raise ValueError('No unique original recipe ' + recipe_name)
            recipe = recipe_dirs[0]
            case = output / row.get('case_name', element + '-' + functional + '-' + version)
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
            if row.get('repair_maxit_format'):
                lines = ps[0].read_text().splitlines(keepends=True)
                fields = lines[2].split()
                intended = int(fields[3])
                lines[2] = (f'{float(fields[0].replace("d", "e")):10.1e}'
                            f'{float(fields[1].replace("d", "e")):10.1e}'
                            f'{float(fields[2]):10.5f}{intended:5d}'
                            '    thresh,tol,damp,maxit (2e10.1,f10.5,i5)\n')
                ps[0].write_text(''.join(lines))
                record['recipe_format_repair'] = 'Honor the source PS input maxit=400 in its documented i5 field; the original misplaced field reads as 4.'
            if functional == 'lda':
                # exfact=0 selects Ceperley-Alder/Perdew-Zunger in USPP 7.3.6.
                lines = ae[0].read_text().splitlines(keepends=True)
                if float(lines[4][15:25]) != 5:
                    raise ValueError('Unexpected source recipe XC selector')
                line = lines[4].rstrip('\r\n')
                lines[4] = line[:15] + f'{0.0:10.5f}' + line[25:] + '\n'
                ae[0].write_text(''.join(lines))
            binary = root / 'sources' / ('vanderbilt-' + ('pbe' if functional == 'lda' else functional)) / 'runatom.x'
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
