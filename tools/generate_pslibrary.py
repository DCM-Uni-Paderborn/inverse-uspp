#!/usr/bin/env python3
"""Generate scalar USPP potentials and full reference pairs from PSLibrary recipes."""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
import concurrent.futures
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import threading
import xml.etree.ElementTree as ET

FUNCTIONALS = {'pz': 'LDA', 'pbe': 'PBE', 'pbesol': 'PBESOL', 'pw91': 'PW91',
               'bp': 'BP', 'revpbe': 'REVPBE', 'wc': 'WC'}


def recipes(source):
    result = {}
    for job in ['us_ps_high.job', 'us_ps_low.job']:
        text = (source / job).read_text()
        for name, body in re.findall(r'^cat > ([^ ]+) << EOF\n(.*?)^EOF$', text, re.M | re.S):
            element = name.split('.')[0]
            result.setdefault(element, []).append((job, name, body))
    return result


def prepare_input(text, functional, full):
    for key, value in {'fct': functional, 'nrel': '1', 'gfun': FUNCTIONALS[functional]}.items():
        text = text.replace('$' + key, value)
    text, count = re.subn(r"(?i)(file_pseudopw\s*=\s*)['\"][^'\"]+['\"]", r"\1'potential.upf'", text)
    assert count == 1
    assert not re.search(r'(?i)\blsave_wfc\s*=', text)
    text, count = re.subn(r'(?i)&inputp\b', '&inputp\n lsave_wfc=' + ('.true.' if full else '.false.') + ',', text)
    assert count == 1
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--generator', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sssp-catalog', type=Path,
                        help='Skip elements already represented by SSSP USPPs')
    parser.add_argument('--functionals', nargs='+', choices=FUNCTIONALS,
                        default=['pbe', 'pbesol'])
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    catalog = json.loads(args.sssp_catalog.read_text()) if args.sssp_catalog else {'entries': []}
    existing = {row['element'] for row in catalog['entries'] if row.get('sssp', True)}
    choices = recipes(args.source)
    commit = subprocess.check_output(['git', '-C', str(args.source), 'rev-parse', 'HEAD'], text=True).strip()
    report = {'source_commit': commit, 'generator_sha256': sha256(args.generator.read_bytes()).hexdigest(),
              'functionals': args.functionals, 'controller_pid': os.getpid(),
              'entries': [], 'completed': False}
    lock = threading.Lock()
    environment = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')

    def run_case(element, functional):
        record = {'element': element, 'functional': 'lda' if functional == 'pz' else functional,
                  'source_functional': functional, 'source_library': 'PSLibrary', 'sssp': False}
        for variant, (job, name, recipe) in enumerate(choices[element]):
            directory = args.output / (element + '-' + functional + '-' + str(variant))
            directory.mkdir(exist_ok=True)
            record.update(recipe_job=job, recipe_name=name.replace('$fct', functional), recipe_sha256=sha256(recipe.encode()).hexdigest())
            failed = False
            for stage, full in [('original', False), ('waves', True)]:
                target = directory / stage
                target.mkdir(exist_ok=True)
                text = prepare_input(recipe, functional, full)
                (target / 'input.in').write_text(text)
                try:
                    with (target / 'input.in').open() as inp, (target / 'generator.log').open('w') as out:
                        rc = subprocess.run([str(args.generator)], cwd=target, env=environment, stdin=inp, stdout=out, stderr=subprocess.STDOUT, timeout=180).returncode
                    if rc or not (target / 'potential.upf').exists():
                        record['failure'] = stage + ' generator failed with status ' + str(rc)
                        failed = True
                        break
                    tree = ET.parse(target / 'potential.upf')
                    header = tree.find('PP_HEADER')
                    if header is None or header.get('element', '').strip() != element:
                        raise ValueError('Wrong generated element')
                    if not float(header.get('z_valence', '0')) > 0:
                        raise ValueError('Nonpositive generated valence charge')
                    if stage == 'waves' and tree.find('PP_FULL_WFC') is None:
                        raise ValueError('No full reference waves')
                except (ValueError, subprocess.TimeoutExpired, ET.ParseError) as exc:
                    record['failure'] = str(exc)
                    failed = True
                    break
            if not failed:
                record.update(status='generated', original=str(directory / 'original/potential.upf'), generated=str(directory / 'waves/potential.upf'), input=str(directory / 'waves/input.in'))
                record['input_sha256'] = sha256((directory / 'waves/input.in').read_bytes()).hexdigest()
                record['original_sha256'] = sha256((directory / 'original/potential.upf').read_bytes()).hexdigest()
                record['generated_sha256'] = sha256((directory / 'waves/potential.upf').read_bytes()).hexdigest()
                record.pop('failure', None)
                break
        else:
            record['status'] = 'failed'
        with lock:
            report['entries'].append(record)
            (args.output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
            print(element, functional, record['status'], flush=True)

    todo = [(element, functional) for element in choices if element not in existing
            for functional in dict.fromkeys(args.functionals)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for _ in pool.map(lambda task: run_case(*task), todo):
            pass
    report['completed'] = True
    (args.output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
