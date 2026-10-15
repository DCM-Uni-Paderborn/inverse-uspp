#!/usr/bin/env python3
"""Preserve original UPF operators while attaching classified reference pairs.

Route 3 requires the recorded replay compatibility bound. Route 4 explicitly
marks an auxiliary reconstruction. The original byte sections are verified
independently after insertion. Atomic invertibility is checked with the same
radial Simpson/trapezoid convention as CP2K; it does not establish a complete
multi-atom inverse, SCF convergence or the absence of reconstruction residuals.
"""
# SPDX-License-Identifier: GPL-2.0-or-later
import argparse
import copy
from hashlib import sha256
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import numpy as np
from scipy.interpolate import CubicSpline
from upf_data import children, compare, flag, top_level_spans, values


def array(node):
    if node is None:
        raise ValueError('Missing radial data')
    return np.asarray(values(node), dtype=float)


def weights(tree):
    r = array(tree.find('PP_MESH/PP_R'))
    rab = array(tree.find('PP_MESH/PP_RAB'))
    if len(r) < 3 or len(rab) != len(r) or np.any(np.diff(r) <= 0) or np.any(rab <= 0):
        raise ValueError('Invalid original radial mesh or quadrature')
    n = len(r)
    ns = n if n % 2 else n - 1
    w = np.zeros(n)
    w[:ns:2] = 2 / 3
    w[1:ns-1:2] = 4 / 3
    w[0] = w[ns-1] = 1 / 3
    if ns < n:
        w[-2:] += 0.5
    return r, w * rab


def projectors(tree):
    n = int(tree.find('PP_HEADER').get('number_of_proj'))
    mesh = len(array(tree.find('PP_MESH/PP_R')))
    beta = np.zeros((mesh,n))
    angular = []
    for i in range(n):
        node = tree.find('PP_NONLOCAL/PP_BETA.' + str(i+1))
        b = array(node)
        if len(b) > mesh:
            raise ValueError('Oversized original projector')
        beta[:len(b),i] = b
        angular.append(int(node.get('angular_momentum')))
    return beta, np.asarray(angular)


def complete_pairs(tree, full):
    header = tree.find('PP_HEADER')
    count = int(header.get('number_of_proj'))
    mesh = int(header.get('mesh_size'))
    if count < 1 or mesh < 3:
        raise ValueError('Empty original data')
    expected = {f'PP_{kind}WFC.{i}' for kind in ['AE','PS'] for i in range(1,count+1)}
    nodes = children(full)
    if nodes.keys() != expected or int(full.get('number_of_wfc',count)) != count:
        raise ValueError('Incomplete or unexpected projector-indexed reference pairs')
    fields = []
    for kind in ['AE','PS']:
        columns = []
        for i in range(1,count+1):
            node = nodes[f'PP_{kind}WFC.{i}']
            beta = tree.find(f'PP_NONLOCAL/PP_BETA.{i}')
            if int(node.get('index',i)) != i or int(node.get('l')) != int(beta.get('angular_momentum')):
                raise ValueError('Reference-wave/projector angular mismatch')
            data = array(node)
            if len(data) != mesh:
                raise ValueError('Wrong reference-wave mesh size')
            columns.append(data)
        fields.append(np.column_stack(columns))
    return fields


def reconstruction_check(tree, full):
    r,w = weights(tree)
    beta,angular = projectors(tree)
    hard,soft = complete_pairs(tree,full)
    delta = hard-soft
    blocks=[]
    for l in sorted(set(angular)):
        ids = np.where(angular==l)[0]
        b = beta[:,ids]
        d = delta[:,ids]
        k = np.eye(len(ids)) + b.T @ (w[:,None]*d)
        condition = float(np.linalg.cond(k,1))
        if not np.isfinite(condition) or 1/condition < 64*np.finfo(float).eps:
            raise ValueError('Singular or insufficiently conditioned atomic inverse')
        error=0.0
        for exponent in [0.2,0.8,2.4]:
            g = r**(int(l)+1)*np.exp(-exponent*r*r)
            coeff = np.linalg.solve(k,b.T@(w*g))
            u = g-d@coeff
            restored = u+d@(b.T@(w*u))
            error=max(error,float(np.linalg.norm(restored-g)/np.linalg.norm(g)))
        if error>1e-8:
            raise ValueError('Atomic inverse round-trip failed')
        blocks.append({'l':int(l),'size':len(ids),'condition_1':condition,'roundtrip_relative':error})
    norm = float(np.linalg.norm(delta))
    return {'atomic_inverse_passed':True,'blocks':blocks,'maximum_condition_1':max(v['condition_1'] for v in blocks),'maximum_roundtrip_relative':max(v['roundtrip_relative'] for v in blocks),'difference_norm':norm,'scope':'Original projector arrays, radial finite-rank atomic inverse and three Gaussian radial probes. Multi-atom inverse, SCF and properties are separate.'}


def insert(original_data, full, route, original_hash, wave_hash):
    tree=ET.fromstring(original_data)
    header=tree.find('PP_HEADER')
    if tree.tag!='UPF' or tree.get('version')!='2.0.1' or header is None:
        raise ValueError('Requires UPF 2.0.1')
    if not flag(header,'is_ultrasoft') or flag(header,'is_paw') or flag(header,'has_so'):
        raise ValueError('Requires scalar ultrasoft original data')
    if tree.find('PP_FULL_WFC') is not None or flag(header,'has_wfc'):
        raise ValueError('Original already contains full-wave data')
    complete_pairs(tree,full)
    full=copy.deepcopy(full)
    full.set('reconstruction_role','compatible_generator_replay' if route==3 else 'auxiliary_reconstruction')
    full.set('reconstruction_route',str(route))
    full.set('original_sha256',original_hash)
    full.set('reference_source_sha256',wave_hash)
    spans,end=top_level_spans(original_data)
    a,b=spans['PP_HEADER']
    old=original_data[a:b]
    if 'has_wfc' in header.attrib:
        updated,count=re.subn(rb'\bhas_wfc\s*=\s*([\x27\"])[^\x27\"]*\1',b'has_wfc="true"',old,count=1)
        if count!=1:raise ValueError('Cannot update full-wave availability')
    else:
        match=re.match(rb'<PP_HEADER\b(?:[^>\"\x27]|\"[^\"]*\"|\x27[^\x27]*\x27)*>',old)
        if match is None:raise ValueError('Invalid original header')
        pos=match.end()-(2 if match.group().endswith(b'/>') else 1)
        updated=old[:pos]+b' has_wfc="true"'+old[pos:]
    note=b'\n<!-- Inverse reference supplement. Original potential/operator/core sections are unchanged. '
    note+=b'Compatible generator replay. -->\n' if route==3 else b'Auxiliary reconstruction, not recovered historical original waves. All original-functional residuals are required. -->\n'
    data=original_data[:a]+updated+original_data[b:end]+note+ET.tostring(full,encoding='utf-8')+b'\n'+original_data[end:]
    result=ET.fromstring(data)
    compare(tree,result)
    new_spans,_=top_level_spans(data)
    for tag,(start,stop) in spans.items():
        if tag!='PP_HEADER':
            x,y=new_spans[tag]
            if original_data[start:stop]!=data[x:y]:raise ValueError('Original section changed: '+tag)
    return data


def prepare_qe(original_data,generated_data,tolerance=1e-6):
    original=ET.fromstring(original_data)
    generated=ET.fromstring(generated_data)
    for tag in ['PP_R','PP_RAB']:
        if not np.array_equal(array(original.find('PP_MESH/'+tag)),array(generated.find('PP_MESH/'+tag))):
            raise ValueError('Original and generated radial conventions differ')
    h1,h2=original.find('PP_HEADER'),generated.find('PP_HEADER')
    for key in ['element','z_valence','functional','number_of_proj']:
        if key=='z_valence':
            same=float(h1.get(key))==float(h2.get(key))
        else:same=h1.get(key).strip().lower()==h2.get(key).strip().lower()
        if not same:raise ValueError('Different reference data: '+key)
    full=generated.find('PP_FULL_WFC')
    if full is None:raise ValueError('Missing exported full-wave pairs')
    info={'compatibility_bound':tolerance}
    try:
        info['maximum_replay_relative']=compare(original,generated,tolerance)
        route=3
    except ValueError as exc:
        info['compatibility_failure']=str(exc)
        route=4
    check=reconstruction_check(original,full)
    data=insert(original_data,full,route,sha256(original_data).hexdigest(),sha256(generated_data).hexdigest())
    return data,dict(info,route=route,reconstruction=check)


def gbrv_pairs(original_data,export):
    tree=ET.fromstring(original_data)
    r,w=weights(tree)
    beta,angular=projectors(tree)
    lines=Path(export).read_text().splitlines()
    nr,nb=map(int,lines[0].split())
    ls=np.asarray([int(x) for x in lines[1].split()])
    if nb!=beta.shape[1] or not np.array_equal(ls,angular):
        raise ValueError('GBRV reference/projector layout differs')
    table=np.loadtxt(lines[3:])
    if table.shape!=(nr,2+3*nb):raise ValueError('Incomplete Vanderbilt export')
    rr=table[:,0]
    if r[0]<rr[0] or r[-1]<rr[-1] or np.any(np.diff(rr)<=0):
        raise ValueError('Incompatible reference radial coverage')
    raw_hard,raw_soft,raw_beta=np.split(table[:,2:],3,axis=1)
    mapped=[]
    for waves in [raw_hard,raw_soft,raw_beta]:
        spline=CubicSpline(rr,waves,axis=0)
        data=np.zeros((len(r),nb))
        inside=r<=rr[-1]
        data[inside]=spline(r[inside])
        mapped.append(data)
    hard,soft,regenerated_beta=mapped
    # Common tails do not alter the finite reconstruction difference. They are
    # an explicit C2 continuation, not recovered historical reference tails.
    common=CubicSpline(rr,0.5*(raw_hard+raw_soft),axis=0)
    x=r[r>rr[-1]]-rr[-1]
    if len(x):
        rate=max(0.5,1/rr[-1])
        v=common(rr[-1]);d=common(rr[-1],1);dd=common(rr[-1],2)
        tail=(v[None,:]+x[:,None]*(d+rate*v)[None,:]+0.5*x[:,None]**2*(dd+2*rate*d+rate**2*v)[None,:])*np.exp(-rate*x[:,None])
        u=np.clip((r[r>rr[-1]]-0.9*r[-1])/(0.1*r[-1]),0,1)
        tail*= (1-10*u**3+15*u**4-6*u**5)[:,None]
        hard[r>rr[-1]]=tail
        soft[r>rr[-1]]=tail
    alignment=np.zeros((nb,nb))
    max_error=0.0
    for l in sorted(set(angular)):
        ids=np.where(angular==l)[0]
        left=regenerated_beta[:,ids]*np.sqrt(w[:,None])
        right=beta[:,ids]*np.sqrt(w[:,None])
        a=np.linalg.lstsq(left,right,rcond=None)[0]
        if not np.isfinite(np.linalg.cond(a)) or np.linalg.cond(a)>1e10:
            raise ValueError('Singular original/generated projector alignment')
        alignment[np.ix_(ids,ids)]=a
        error=np.linalg.norm(left@a-right)/np.linalg.norm(right)
        max_error=max(max_error,float(error))
    inverse_transpose=np.linalg.inv(alignment).T
    hard=hard@inverse_transpose
    soft=soft@inverse_transpose
    full=ET.Element('PP_FULL_WFC',{'number_of_wfc':str(nb)})
    for i in range(nb):
        for kind,waves in [('AE',hard),('PS',soft)]:
            node=ET.SubElement(full,f'PP_{kind}WFC.{i+1}',{'type':'real','size':str(len(r)),'columns':'4','index':str(i+1),'l':str(int(angular[i]))})
            node.text='\n'+'\n'.join(' '.join(f'{x:.16e}' for x in waves[j:j+4,i]) for j in range(0,len(r),4))+'\n'
    info={'projector_alignment_relative':max_error,'projector_alignment':alignment.tolist(),'reference_export_radius':float(rr[-1]),'outer_continuation':'Identical C2 decaying continuation outside the generator export. The original reconstruction difference is unchanged there.','reference_mapping':'Generator pairs transformed with U^-T after qdiag, interpolated onto the original radial mesh, then aligned with the original projector basis.'}
    return full,info


def prepare_gbrv(original_data,export):
    full,info=gbrv_pairs(original_data,export)
    tree=ET.fromstring(original_data)
    check=reconstruction_check(tree,full)
    # The complete historical operator replay has not been established. A
    # projector fit is not sufficient to classify a historical dataset match.
    data=insert(original_data,full,4,sha256(original_data).hexdigest(),sha256(Path(export).read_bytes()).hexdigest())
    return data,dict(info,route=4,compatibility_failure='Complete historical operator replay is not established. Original arrays are retained.',reconstruction=check)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('original',type=Path)
    p.add_argument('reference',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--format',choices=['qe','vanderbilt'],default='qe')
    a=p.parse_args()
    if a.output.exists():raise SystemExit('Output already exists')
    raw=a.original.read_bytes()
    result,info=prepare_qe(raw,a.reference.read_bytes()) if a.format=='qe' else prepare_gbrv(raw,a.reference)
    a.output.write_bytes(result)
    print(json.dumps(info,indent=2))


if __name__=='__main__':
    main()
