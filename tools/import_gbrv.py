#!/usr/bin/env python3
"""Add the complete current GBRV UPF libraries with auxiliary reference pairs."""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
from collections import Counter
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
from urllib.request import urlopen
import xml.etree.ElementTree as ET

from prepare_upf import prepare_gbrv
from upf_data import values


ARCHIVES = {
    'pbe': '4ebfeaa4c9d916c1ac2da6e89ec888fdffce096714f1956db01fea5ec9c9e4f0',
    'pbesol': '9254d3a9796f9506a5de835ec1823085b1f5893db1188a5ac0360ffd30b3c37d',
    'lda': '327b0ab395660bd6fa11f1d1004be12abdd0fc5077c709cb26d21634ef7ed0db',
}
NA_RECIPE = {
    'Makefile': '06f0ace04713e54dddacce6966c18c3f05e59be0ba3d67a6bf9929a4d21e73ee',
    'Na_ae_ref.adat': 'f172622ab0bcf13d765171aaa50455ce6a55acaa4abc94a14886b6b9a9738dc6',
    'Na_ps.adat': '2330f49880eb0817634cda86da3859d2dc2a127ed08544ba42c7e032d9e74ce4',
}
EXTRA_UPFS = {
    'lda': {'bi_lda_v1.uspp.F.UPF': '0517346438f71e9f0ab8cc51af7568b5d2743c17f899d84fc195a0f7dd657db9'},
    'pbesol': {'bi_pbesol_v1.uspp.F.UPF': '729ddbbc3a461329c3d46ae377edc71ea6ac4178a90c471b16f7d10d082ef4e2'},
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def download(url, path, expected):
    if not path.exists():
        with urlopen(url, timeout=90) as response:
            data = response.read()
        if sha256(data).hexdigest() != expected:
            raise ValueError('Source checksum mismatch: ' + url)
        path.write_bytes(data)
    if digest(path) != expected:
        raise ValueError('Existing source checksum mismatch: ' + str(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--collection', type=Path)
    parser.add_argument('--converter', type=Path, required=True,
                        help='convert_gbrv executable linked with QE modules')
    parser.add_argument('--functionals', nargs='+', choices=ARCHIVES,
                        default=['pbe', 'pbesol', 'lda'])
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--retry-failed', action='store_true')
    parser.add_argument('--reconvert', action='store_true',
                        help='Recreate derived UPF v2 files from retained v1 sources')
    args = parser.parse_args()
    root = args.root.resolve()
    collection = args.collection or root / 'collection'
    catalog = json.loads((collection / 'catalog.json').read_text())
    work = root / 'work/gbrv-complete'
    sources = root / 'sources/gbrv-complete'
    work.mkdir(parents=True, exist_ok=True)
    sources.mkdir(parents=True, exist_ok=True)
    # The aggregate recipe archive predates the separately published Na v1.5.
    sodium = root / 'sources/gbrv-recipes/Work_GBRV/011-Na/na_pbe_v1.5'
    sodium.mkdir(parents=True, exist_ok=True)
    for name, expected in NA_RECIPE.items():
        download('https://www.physics.rutgers.edu/gbrv/011-Na/na_pbe_v1.5/' + name,
                 sodium / name, expected)
    records = []
    archives = []
    for functional in dict.fromkeys(args.functionals):
        name = 'all_' + functional + '_UPF_v1.5.tar.gz'
        url = 'https://www.physics.rutgers.edu/gbrv/' + name
        archive = sources / name
        download(url, archive, ARCHIVES[functional])
        archives.append(dict(functional=functional, url=url, sha256=digest(archive)))
        with tarfile.open(archive) as handle:
            members = [m for m in handle.getmembers()
                       if m.isfile() and m.name.upper().endswith('.UPF')]
            raw_files = {Path(m.name).name: handle.extractfile(m).read() for m in members}
        for filename, expected in EXTRA_UPFS.get(functional, {}).items():
            directory = sources / functional / filename
            directory.mkdir(parents=True, exist_ok=True)
            original = directory / filename
            download('https://www.physics.rutgers.edu/gbrv/' + filename, original, expected)
            raw_files[filename] = original.read_bytes()
        for filename, raw in sorted(raw_files.items()):
            directory = sources / functional / filename
            directory.mkdir(parents=True, exist_ok=True)
            original = directory / filename
            if original.exists() and original.read_bytes() != raw:
                raise ValueError('Conflicting original source: ' + filename)
            original.write_bytes(raw)
            converted = original.with_name(original.name + '.UPF')
            if args.reconvert or not converted.exists():
                with (directory / 'conversion.log').open('w') as log:
                    subprocess.run([str(args.converter), filename], cwd=directory,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
            # CP2K's line reader needs a separate PP_INFO closing tag. Retain
            # the original author's information instead of synthesized text.
            converted_data = converted.read_bytes()
            info = re.search(rb'<PP_INFO>(.*?)</PP_INFO>', raw, re.S).group(1)
            if b']]>' in info:
                raise ValueError('Unsupported CDATA terminator in original information')
            converted_data = re.sub(rb'<PP_INFO>.*?</PP_INFO>',
                                    lambda m: b'<PP_INFO>\n<![CDATA[' + info + b']]>\n</PP_INFO>',
                                    converted_data, count=1, flags=re.S)
            converted.write_bytes(converted_data)
            header = ET.fromstring(converted.read_bytes()).find('PP_HEADER')
            version = re.search(r'(v[0-9.]+)\.uspp', filename).group(1)
            recipe = filename.split('.uspp')[0].replace('_' + functional + '_', '_pbe_')
            records.append(dict(
                element=header.get('element').strip(), functional=functional,
                version=version, recipe_name=recipe, case_name=filename.split('.uspp')[0],
                source_library='GBRV', sssp=False, original_path=str(converted),
                original_sha256=digest(converted), upstream_upf_v1_sha256=digest(original),
                upstream_filename=filename,
                upstream_url='https://www.physics.rutgers.edu/gbrv/' + filename,
                upstream_archive_sha256=None if filename in EXTRA_UPFS.get(functional, {}) else digest(archive),
                z_valence=float(header.get('z_valence')),
                repair_maxit_format=recipe == 'ta_pbe_v1'))
    manifest = work / 'inputs.json'
    manifest.write_text(json.dumps({'entries': records}, indent=2) + '\n')
    environment = os.environ.copy()
    environment.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    arguments = [sys.executable, str(Path(__file__).with_name('generate_gbrv.py')),
                 '--root', str(root), '--manifest', str(manifest),
                 '--output', str(work / 'exports'), '--workers', str(args.workers)]
    if args.retry_failed:
        arguments.append('--retry-failed')
    subprocess.run(arguments,
                   env=environment, check=True)
    exports = json.loads((work / 'exports/result.json').read_text())
    failures = [r for r in exports['entries'] if r['status'] != 'generated']
    if failures:
        raise ValueError('Incomplete reference exports: ' + repr([
            (r['upstream_filename'], r.get('failure')) for r in failures]))
    existing = {r['file']: r for r in catalog['entries']}
    prepared = []
    for source in sorted(exports['entries'], key=lambda r: r['upstream_filename']):
        original = Path(source['original_path'])
        reference = Path(source['export'])
        raw = original.read_bytes()
        tree = ET.fromstring(raw)
        directory = Path(source['directory'])
        atomic_z = float(re.search(r'\bz\s*=\s*([0-9.]+)',
                                  (directory / 'ae.out').read_text(errors='replace')).group(1))
        core = float(re.search(r'core charge\s*=\s*([0-9.]+)',
                              (directory / 'ps.out').read_text(errors='replace')).group(1))
        if abs(atomic_z - core - source['z_valence']) > 1e-4:
            raise ValueError('Reference/original valence partition differs')
        data, details = prepare_gbrv(raw, reference)
        relative = Path('data') / source['functional'] / 'additional/gbrv' / source['upstream_filename']
        target = collection / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.read_bytes() != data:
            raise ValueError('Refusing to replace different prepared data: ' + str(target))
        target.write_bytes(data)
        recipe = Path('inputs/gbrv') / source['recipe_name'] / source['functional']
        (collection / recipe).mkdir(parents=True, exist_ok=True)
        for key, filename in [('ae_input', 'ae.adat'), ('ps_input', 'ps.adat')]:
            content = Path(source[key]).read_bytes()
            destination = collection / recipe / filename
            if destination.exists() and destination.read_bytes() != content:
                raise ValueError('Conflicting retained atomic input: ' + str(destination))
            destination.write_bytes(content)
        record = {k: source[k] for k in (
            'element', 'functional', 'source_library', 'sssp', 'original_sha256',
            'upstream_upf_v1_sha256', 'upstream_filename', 'upstream_url', 'upstream_archive_sha256',
            'z_valence', 'generator_sha256')}
        record.update(details, status='prepared', file=str(relative), sha256=sha256(data).hexdigest(),
                      original_path=str(original.relative_to(root)), bytes=len(data),
                      reference_source_sha256=digest(reference),
                      generator='Vanderbilt ultrasoft 7.3.6 with reference export',
                      generator_input=str(recipe),
                      format_conversion='QE 6.3MaX read_upf_v1, set_pseudo_upf and write_upf; original inner augmentation polynomials expanded to QIJL without resampling',
                      converter_sha256=digest(args.converter),
                      nlcc=any(values(tree.find('PP_NLCC'))) if tree.find('PP_NLCC') is not None else False,
                      model_core_tau_added=False, complete_projector_pairs=True,
                      original_operator_sections_unchanged=True,
                      cp2k_native_data_check_passed=False,
                      full_cp2k_calculation_qualified=False)
        if source['functional'] == 'lda':
            record['reference_recipe_change'] = 'AE input exfact=0 selects Ceperley-Alder/Perdew-Zunger; original UPF operators are retained.'
        for key in ('generator_stopped_after_reference_export', 'generator_exit_status',
                    'generator_diagnostic', 'recipe_format_repair'):
            if key in source:
                record[key] = source[key]
        previous = existing.get(str(relative))
        if previous is not None and previous['sha256'] != record['sha256']:
            raise ValueError('Catalog already records different prepared data')
        existing[str(relative)] = record
        prepared.append(record)
        print(source['upstream_filename'], 'prepared', flush=True)
    catalog['entries'] = list(existing.values())
    catalog['summary'] = dict(
        datasets=len(existing), prepared=len(existing), failed=0,
        elements=len({r['element'] for r in existing.values()}),
        routes=dict(Counter(str(r['route']) for r in existing.values())),
        upf_bytes=sum(r['bytes'] for r in existing.values()))
    catalog['status'] = 'prepared_pending_native_checks'
    catalog['gbrv_libraries'] = dict(version='1.5', archives=archives,
                                   counts=dict(Counter(r['functional'] for r in prepared)))
    (collection / 'catalog.json').write_text(json.dumps(catalog, indent=2) + '\n')
    print(json.dumps(catalog['summary']), flush=True)


if __name__ == '__main__':
    main()
