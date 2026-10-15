#!/usr/bin/env python3
"""Assemble preserved USPPs and reference-wave provenance from generator results."""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
import xml.etree.ElementTree as ET
from prepare_upf import prepare_qe, prepare_gbrv


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--sssp-root',type=Path,required=True)
    parser.add_argument('--replay-root',type=Path,required=True)
    args=parser.parse_args()
    root=args.root
    source_catalog=json.loads((root/'catalog.json').read_text())
    output=root/'collection'
    output.mkdir(exist_ok=True)
    additional=json.loads((root/'work/pslibrary-additional/result.json').read_text())
    gbrv=json.loads((root/'work/gbrv/result.json').read_text())
    original_replays=json.loads((args.replay_root/'attempt3/result.json').read_text())
    by_hash={r['original_sha256']:r for r in original_replays['records']}
    gbrv_by_hash={r['original_sha256']:r for r in gbrv['entries']}
    manifest={'schema_version':1,'status':'preparing','scope':'Original scalar USPP operators with compatible replay or explicitly auxiliary reference waves. Atomic reconstruction checks are not complete system/SCF/property qualification.','entries':[]}
    generator=Path('/home/kuehne88/work/uspp-inverse-data/qe-6.3MaX/atomic/src/ld1.x')
    qe_hash=sha256(generator.read_bytes()).hexdigest()

    def save():
        good=[r for r in manifest['entries'] if r['status']=='prepared']
        manifest['summary']={'datasets':len(manifest['entries']),'prepared':len(good),'elements':len({r['element'] for r in good}),'routes':dict(Counter(str(r['route']) for r in good)),'failed':len(manifest['entries'])-len(good),'upf_bytes':sum(r['bytes'] for r in good)}
        path=output/'catalog.json'
        temporary=path.with_suffix('.writing')
        temporary.write_text(json.dumps(manifest,indent=2)+'\n')
        temporary.replace(path)

    def prepare(row,original,reference,kind,source_record):
        raw=original.read_bytes()
        digest=sha256(raw).hexdigest()
        record=dict(row)
        record['original_sha256']=digest
        try:
            tree=ET.fromstring(raw)
            header=tree.find('PP_HEADER')
            element=header.get('element').strip()
            functional=header.get('functional').strip().lower()
            q=float(header.get('z_valence'))
            if row['element']!=element:raise ValueError('Original element mismatch')
            if kind=='qe':
                result,details=prepare_qe(raw,reference.read_bytes())
                recipe=tree.find('.//PP_INPUTFILE')
                record['generator_input']='Original PP_INFO/PP_INPUTFILE with lsave_wfc enabled' if recipe is not None else source_record['input']
                record['generator']='Quantum ESPRESSO atomic 6.3MaX'
                record['generator_sha256']=qe_hash
            else:
                directory=Path(source_record['directory'])
                ae_output=(directory/'ae.out').read_text(errors='replace')
                ps_output=(directory/'ps.out').read_text(errors='replace')
                atomic_z=float(re.search(r'\bz\s*=\s*([0-9.]+)',ae_output).group(1))
                core=float(re.search(r'core charge\s*=\s*([0-9.]+)',ps_output).group(1))
                if abs(atomic_z-core-q)>1e-4:raise ValueError('Generator/original valence partition mismatch')
                result,details=prepare_gbrv(raw,reference)
                input_dir=Path('inputs/gbrv')/source_record['recipe_directory'].split('/')[-1]
                (output/input_dir).mkdir(parents=True,exist_ok=True)
                for key,filename in [('ae_input','ae.adat'),('ps_input','ps.adat')]:
                    src=Path(source_record[key])
                    target=output/input_dir/filename
                    if target.exists() and target.read_bytes()!=src.read_bytes():raise ValueError('Conflicting generator recipe')
                    target.write_bytes(src.read_bytes())
                record['generator_input']=str(input_dir)
                record['generator']='Vanderbilt ultrasoft 7.3.6 with reference export'
                record['generator_sha256']=source_record['generator_sha256']
                if source_record.get('generator_stopped_after_reference_export'):
                    record['generator_stopped_after_reference_export']=True
                    record['generator_exit_status']=source_record['generator_exit_status']
                    record['generator_diagnostic']=source_record['generator_diagnostic'].strip()
            if record['sssp']:
                filename=original.name
                group='sssp'
            else:
                filename=f'{element}.us.{functional}.z_{q:g}.ld1.psl.v1.0.0-high.upf'
                group='additional'
                record['recipe_job']=source_record['recipe_job']
                record['recipe_name']=source_record['recipe_name']
                record['recipe_sha256']=source_record['recipe_sha256']
                record['pslibrary_commit']=additional['source_commit']
            relative=Path('data')/functional/group/filename
            target=output/relative
            target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists() and target.read_bytes()!=result:raise ValueError('Output filename collision')
            target.write_bytes(result)
            core=tree.find('PP_NLCC')
            core_nonzero=core is not None and any(float(v.replace('D','E').replace('d','e'))!=0 for v in (core.text or '').split())
            record.update(details,status='prepared',functional=functional,z_valence=q,file=str(relative),sha256=sha256(result).hexdigest(),reference_source_sha256=sha256(reference.read_bytes()).hexdigest(),bytes=len(result),nlcc=core_nonzero,model_core_tau_added=False,complete_projector_pairs=True,original_operator_sections_unchanged=True,full_cp2k_calculation_qualified=False)
            record.pop('failure',None)
            record.pop('directory',None)
        except Exception as exc:
            record.update(status='failed',failure=type(exc).__name__+': '+str(exc))
        manifest['entries'].append(record)
        save()
        print(record['element'],record.get('functional'),record['status'],record.get('route'),record.get('failure',''),flush=True)

    for old in source_catalog['entries']:
        row=dict(old,sssp=True)
        original=args.sssp_root/row['original_path']
        if sha256(original.read_bytes()).hexdigest()!=row['original_sha256']:raise ValueError('Original dataset hash changed')
        if row['source_library']=='PSLibrary':
            source=by_hash[row['original_sha256']]
            directory=args.replay_root/'attempt2'/(row['element']+'-'+('pbesol' if '.pbesol.' in original.name else 'pbe')+'-'+row['original_sha256'][:10])
            generated=directory/'generated.upf'
            if sha256(generated.read_bytes()).hexdigest()!=source['generated_sha256']:raise ValueError('Reference export hash changed')
            prepare(row,original,generated,'qe',source)
        else:
            source=gbrv_by_hash[row['original_sha256']]
            if source['status']!='generated':
                manifest['entries'].append(dict(row,status='failed',failure=source.get('failure','No generator export')))
                save()
                continue
            prepare(row,original,Path(source['export']),'vanderbilt',source)
    for source in additional['entries']:
        row={k:source[k] for k in ['element','functional','source_library','sssp']}
        if source['status']!='generated':
            manifest['entries'].append(dict(row,status='failed',failure=source['failure']))
            save()
        else:
            prepare(row,Path(source['original']),Path(source['generated']),'qe',source)
    manifest['status']='prepared' if manifest['summary']['failed']==0 else 'partially_prepared'
    save()
    print(json.dumps(manifest['summary']),flush=True)


if __name__=='__main__':
    main()
