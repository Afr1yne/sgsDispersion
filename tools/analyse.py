#!/usr/bin/env python3
"""Read real OpenFOAM VTK particles; retain every result, without success filtering."""
import csv, hashlib, json, math, re, statistics, sys, os
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT=Path(os.environ['FOAM_RUN'])/'sgsDispersion'
NAMES=('RANS_DRW','LES_resolved','LES_SGS')
def read_vtp(path):
    tree=ET.parse(path)
    result={}
    for array in tree.findall('.//DataArray'):
        values=[float(x) for x in (array.text or '').split()]
        components=int(array.attrib.get('NumberOfComponents','1'))
        result[array.attrib['Name']]=values if components==1 else [values[i:i+components] for i in range(0,len(values),components)]
    return result
def moment(values):
    mean=statistics.fmean(values)
    return mean, math.sqrt(statistics.fmean([(x-mean)**2 for x in values]))
def foam_values(path,components=1):
    text=path.read_text()
    match=re.search(r'internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    if not match:
        uniform=re.search(r'internalField\s+uniform\s+([^;]+);',text)
        if not uniform:return []
        return [float(x) for x in uniform.group(1).strip(' ()').split()]
    return [float(x) for x in re.findall(r'[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?',match.group(2))]

def analyse(suite):
    base=ROOT/suite
    manifest=json.loads((base/'manifest.json').read_text())
    rows=[]; last={}; particles={}; checks={}
    for name in NAMES:
        case=base/name
        assert not re.search(r'\bfunctions\s*\{',(case/'system/controlDict').read_text()), 'Unexpected function objects'
        log=(case/'log.kinematicParcelFoam').read_text()
        assert '\nEnd\n' in log, f'{name} did not finish'
        assert 'FOAM FATAL' not in log
        assert 'Mesh OK' in (case/'log.checkMesh').read_text()
        courant=[float(v) for v in re.findall(r'Courant Number mean: \S+ max: (\S+)',log)]
        checks[name]={'finished':True,'mesh_ok':True,'max_courant':max(courant),
                      'warnings':log.count('FOAM Warning')}
        folder=case/'VTK/lagrangian/kinematicCloud'
        series=json.loads((folder/'kinematicCloud.vtp.series').read_text())
        for entry in series['files']:
            data=read_vtp(folder/entry['name'])
            assert all(math.isfinite(v) for xyz in data['Points'] for v in xyz)
            for d in manifest['diameters_m']:
                ids=[i for i,v in enumerate(data['d']) if abs(v-d)<d*1e-4 and data['active'][i]>0]
                if not ids:continue
                coords=[data['Points'][i] for i in ids]
                x,sx=moment([p[0] for p in coords]); y,sy=moment([p[1] for p in coords]); z,sz=moment([p[2] for p in coords])
                velocities=[data['UTurb'][i] for i in ids]
                rms=[math.sqrt(statistics.fmean([u[j]**2 for u in velocities])) for j in range(3)]
                ages=[data['age'][i] for i in ids]
                row={'case':name,'time_s':entry['time'],'age_s':statistics.fmean(ages),
                     'time_after_scheduled_release_s':entry['time']-manifest['release_s'],
                     'particle_age_min_s':min(ages),'particle_age_max_s':max(ages),
                     'diameter_um':round(d*1e6),'active_count':len(ids),
                     'mean_x_m':x,'mean_y_m':y,'mean_z_m':z,
                     'sigma_x_m':sx,'sigma_y_m':sy,'sigma_z_m':sz,
                     'sgs_rms_x_m_s':rms[0],'sgs_rms_y_m_s':rms[1],'sgs_rms_z_m_s':rms[2]}
                rows.append(row)
                last[(name,round(d*1e6))]=row
                particles[(name,round(d*1e6))]=coords
    # All saved LES carrier states must match. A particle-only closure must not change
    # the one-way-coupled carrier realization or the inflow random sequence.
    les0,les1=base/'LES_resolved',base/'LES_SGS'
    same={}; differences=[]
    times=sorted((p for p in les0.iterdir() if p.is_dir() and re.fullmatch(r'\d+(\.\d+)?',p.name)),key=lambda p:float(p.name))
    assert {p.name for p in times} == {p.name for p in les1.iterdir() if p.is_dir() and re.fullmatch(r'\d+(\.\d+)?',p.name)}, 'LES saved time directories differ'
    for time in times:
        for field in ('U','p','k','nut'):
            a,b=time/field,les1/time.name/field
            assert a.is_file() and b.is_file(), f'Missing carrier field {time.name}/{field}'
            if a.is_file() and b.is_file():
                equal=hashlib.sha256(a.read_bytes()).digest()==hashlib.sha256(b.read_bytes()).digest()
                same[f'{time.name}/{field}']=equal
                if not equal:differences.append(f'{time.name}/{field}')
    if differences:raise RuntimeError('LES carrier realizations differ: '+','.join(differences[:8]))
    assert same, 'No carrier comparison was performed'
    rans_zero=all(r['sgs_rms_x_m_s']==r['sgs_rms_y_m_s']==r['sgs_rms_z_m_s']==0 for r in rows if r['case']=='RANS_DRW')
    no_sgs=all(r['sgs_rms_x_m_s']==r['sgs_rms_y_m_s']==r['sgs_rms_z_m_s']==0 for r in rows if r['case']=='LES_resolved')
    corrected_nonzero=any(r['sgs_rms_y_m_s']>0 for r in rows if r['case']=='LES_SGS')
    assert no_sgs and corrected_nonzero, 'Dispersion model integration failed'
    latest=base/'LES_SGS'/str(int(manifest['end_s']))
    sigma=foam_values(latest/'kinematicCloudSGSVariance')
    alpha=foam_values(latest/'kinematicCloudSGSFraction')
    assert sigma and alpha and min(sigma)>=0 and min(alpha)>=0 and max(alpha)<=1
    duration=foam_values(latest/'kinematicCloudSGSWindowDuration')
    resolvedK=foam_values(latest/'kinematicCloudSGSResolvedK')
    assert duration and resolvedK and min(resolvedK)>=0
    if suite=='abl': assert min(duration)>=399.999999, 'Internal window not filled'
    out=ROOT/'results'/suite;out.mkdir(parents=True,exist_ok=True)
    with (out/'dispersion.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary={'manifest':manifest,'run_checks':checks,'identical_LES_carrier_fields':len(same),
             'all_LES_carrier_fields_identical':all(same.values()),
             'RANS_UTurb_zero_for_all_saved_particles':rans_zero,
             'LES_without_SGS_UTurb_zero':no_sgs,'SGS_nonzero_fluctuation':corrected_nonzero,
             'sigma2_range': [min(sigma),max(sigma)],'alpha_range':[min(alpha),max(alpha)],
             'internal_window_duration_s':[min(duration),max(duration)],
             'internal_resolved_k_range':[min(resolvedK),max(resolvedK)],
             'no_controlDict_functions_in_all_cases':True,
             'final_statistics':list(last.values())}
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    (out/'final_particles.json').write_text(json.dumps({f'{n}_{d}':p for (n,d),p in particles.items()}),encoding='utf-8')
    lines=['# Actual development-run results','',f'Suite: {suite}; release {manifest["release_s"]} s; end {manifest["end_s"]} s.',
           '',f'Identical LES carrier states: {len(same)} U/p/k/nut files.','',
           '| Case | diameter (um) | count | sigma y (m) | sigma z (m) |',
           '|---|---:|---:|---:|---:|']
    lines += [f'| {r["case"]} | {r["diameter_um"]} | {r["active_count"]} | {r["sigma_y_m"]:.4f} | {r["sigma_z_m"]:.4f} |' for r in last.values()]
    lines += ['', f'Final OpenFOAM particle age: {list(last.values())[0]["age_s"]} s; scheduled-release elapsed time: {manifest["end_s"]-manifest["release_s"]} s.',
              'These are pulse-cloud widths about the instantaneous cloud centroid.',
              'They do not include the meandering of a continuously released plume.',
              'RANS uses frozen prescribed profiles; LES carrier flow evolves.',
              'No reference experiment, grid convergence, or multi-seed convergence is established.']
    (out/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':analyse(sys.argv[1] if len(sys.argv)>1 else 'smoke')
