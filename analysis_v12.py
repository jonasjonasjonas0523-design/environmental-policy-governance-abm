"""Post-principal analysis tail (all report-only; nothing alters the frozen
41x41 preregistered results):
1. A_min sensitivity at {0.02, 0.10} x A*_cum (= {8, 40}): recount of CE-defined
   cells (mean A >= threshold) per plane x instrument on existing raw outputs.
2. Higher-resolution S/T-only boundary check on the reference plane at 81x81
   (next doubling step of the pilot ladder 6-11-21-41), principal seeds, scored
   with the pilot's own nearest-matched crossing metric (all Delta_ST zero
   crossings on 41 reference lines per axis; tol = 0.5 x finer spacing).
   POST-PRINCIPAL ROBUSTNESS ONLY; the 41x41 grid did not formally converge and
   no convergence is claimed here.
3. Connected-component crossover-band extraction on the final joint bands:
   4-connectivity components of anchored transition cells adjacent to BOTH a
   supported-positive and a supported-negative anchored cell; per component:
   cell count, p/c bounding box, width = n_cells / max cell-extent along the
   longer axis (mean thickness; definition documented here).
4. FINAL_READOUTS.md consolidation.
"""
import numpy as np, pandas as pd, json, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from scipy import ndimage
from model_v12 import CONFIG as CFG, batch, run_standard, run_charge
from pipeline import all_crossings, cross_disp

RES='results_v12'
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
n = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(100000,100000+n))
Rg=41; pp=np.linspace(0,1,Rg)
PLANES=[('reference',CFG['ref']['k'],CFG['ref']['h']),('k010',0.10,CFG['ref']['h']),
        ('k055',0.55,CFG['ref']['h']),('h050',CFG['ref']['k'],0.50),('h095',CFG['ref']['k'],0.95)]
pairs=[('S-T',0,1),('S-V',0,2),('T-V',1,2)]

def a_min_sensitivity():
    rows=[]
    for pl,_,_ in PLANES:
        out=np.load(f'{RES}/principal/plane_{pl}_out.npy')   # [n,cells,3,2]
        mA=out[:,:,:,0].mean(0)                              # [cells,3]
        for ji,j in enumerate('STV'):
            r=dict(plane=pl,instrument=j)
            for tag,th in (('A_min_8',8.0),('A_min_20',20.0),('A_min_40',40.0)):
                r[tag+'_defined']=int((mA[:,ji]>=th).sum())
            rows.append(r)
    df=pd.DataFrame(rows); df.to_csv(f'{RES}/a_min_sensitivity.csv',index=False)
    return df

def hires_check():
    f=f'{RES}/hires_boundary_check.json'
    if os.path.exists(f): return json.load(open(f))
    b=batch(seeds); R2=81; pp2=np.linspace(0,1,R2)
    d2=np.zeros((n,R2,R2,2)); t0=time.time()
    for ip in range(R2):
        for ic in range(R2):
            g=dict(p=pp2[ip],c=pp2[ic],k=CFG['ref']['k'],h=CFG['ref']['h'])
            oS=run_standard(b,g,TAU); oT=run_charge(b,g,TAU)
            d2[:,ip,ic,0]=oS['A']-oT['A']; d2[:,ip,ic,1]=oS['C_social']-oT['C_social']
    dt=time.time()-t0
    out41=np.load(f'{RES}/principal/plane_reference_out.npy')
    d1=(out41[:,:,0,:]-out41[:,:,1,:]).reshape(n,Rg,Rg,2)
    res={}
    for mi,m in enumerate(('A','C_social')):
        m1=d1[...,mi].mean(0); s1=d1[...,mi].std(0,ddof=1)/np.sqrt(n)
        m2=d2[...,mi].mean(0); s2=d2[...,mi].std(0,ddof=1)/np.sqrt(n)
        cr1=all_crossings(m1,s1,pp); cr2=all_crossings(m2,s2,pp2)
        disp,matched,mismatch=cross_disp(cr1,cr2)
        res[m]=dict(max_boundary_shift=disp,matched_lines=matched,mismatched_lines=mismatch)
    tol=0.5/(R2-1)
    shifts=[v['max_boundary_shift'] for v in res.values() if v['max_boundary_shift'] is not None]
    rep=dict(from_res=41,to_res=81,seeds='principal 100000..100152 (n=153)',
             per_metric=res,tol=tol,
             max_shift=float(np.max(shifts)) if shifts else None,
             within_tol=bool(shifts and max(shifts)<tol),
             seconds=round(dt,1),
             status=('POST-PRINCIPAL ROBUSTNESS ONLY: the preregistered 41x41 grid was a '
                     'fallback to the finest candidate and did NOT formally converge; this '
                     'S/T-only 81x81 check does not replace the 41x41 result and no formal '
                     'convergence is claimed.'))
    json.dump(rep,open(f,'w'),indent=1)
    return rep

def crossover_components():
    xref=pd.read_csv(f'{RES}/principal/crossover_summary_final_all_planes.csv')
    struct=np.array([[0,1,0],[1,1,1],[0,1,0]])   # 4-connectivity
    rows=[]
    for pl,_,_ in PLANES:
        anc=np.load(f'{RES}/principal/plane_{pl}_anchor.npz')['union'].astype(bool)
        bl=pd.read_parquet(f'{RES}/principal/plane_{pl}_bands_final_joint.parquet')
        for pn,_,_ in pairs:
            for m in ('A','C_social'):
                sup=xref[(xref['plane']==pl)&(xref['pair']==pn)&(xref['outcome']==m)]
                if not bool(sup['supported_crossover'].iloc[0]): continue
                r=bl[(bl['pair']==pn)&(bl['metric']==m)].sort_values('cell')
                lo=r['lo'].values.reshape(Rg,Rg); hi=r['hi'].values.reshape(Rg,Rg)
                pos=(lo>0)&anc; neg=(hi<0)&anc; tr=(lo<=0)&(hi>=0)&anc
                labarr,ncomp=ndimage.label(tr,structure=struct)
                for ci in range(1,ncomp+1):
                    comp=labarr==ci
                    grown=ndimage.binary_dilation(comp,structure=struct)
                    if not ((grown&pos).any() and (grown&neg).any()): continue
                    ii,jj=np.where(comp)
                    ep=ii.max()-ii.min()+1; ec=jj.max()-jj.min()+1
                    rows.append(dict(plane=pl,pair=pn,metric=m,component=ci,
                        n_cells=int(comp.sum()),
                        p_min=float(pp[ii.min()]),p_max=float(pp[ii.max()]),
                        c_min=float(pp[jj.min()]),c_max=float(pp[jj.max()]),
                        extent_p_cells=int(ep),extent_c_cells=int(ec),
                        width_cells=float(comp.sum()/max(ep,ec))))
    df=pd.DataFrame(rows)
    df.to_csv(f'{RES}/crossover_components.csv',index=False)
    return df

def final_readouts(amin,hires,comps):
    lines=['# FINAL READOUTS (v1.2, all planes, final paper-wide joint bands)\n\n'
           'Scope: 41x41 p-x-c planes at reference (k=0.30,h=0.75) and slices k in {0.10,0.55}, '
           'h in {0.50,0.95}; n=153 principal seeds; paper-wide joint max-t family '
           '(crit_A=3.1597, crit_C=3.5874); anchored = G_anchor 95 percent HDR union.\n\n'
           '## Per-plane classification (counts of 1681 cells)\n\n'
           '| plane | anchored | T_dominant | mixed/no-dominance | unresolved | A-leader T (anchored) | C-leader S (anchored) |\n'
           '|---|---|---|---|---|---|---|\n']
    for pl,_,_ in PLANES:
        s=pd.read_csv(f'{RES}/principal/plane_{pl}_summary_final.csv')
        c=s.drop_duplicates(subset=['p','c'])
        anc=c['anchored'].astype(bool)
        lines.append(f"| {pl} | {int(anc.sum())} | {int((c.cell_label=='T_dominant').sum())} | "
            f"{int((c.cell_label=='mixed').sum())} | {int((c.cell_label=='unresolved').sum())} | "
            f"{float((c.A_leader[anc]=='T').mean()):.2f} | {float((c.C_leader[anc]=='S').mean()):.2f} |\n")
    lines.append('\nZero cells carry S_dominant or V_dominant labels on any plane. '
        'H1 6/6 pair-metric combinations, H2 3/3 pairs, H3 SUPPORTED (within-plane; reference, '
        'h050, h095) per H_verdicts_final.json.\n')
    lines.append('\n## Supported crossover bands (connected components, 4-connectivity)\n\n'
        'Component = anchored transition cells adjacent to both supported signs; width = '
        'n_cells / max cell extent along the longer axis.\n\n')
    lines.append(comps.to_markdown(index=False)+'\n')
    lines.append('\n## A_min sensitivity (report-only; PARAMETER_PROVENANCE row)\n\n'
        'CE-defined cells (mean A >= A_min) per plane x instrument at A_min in {8, 20, 40} '
        '(= {0.02, 0.05, 0.10} x A*_cum). Classification never uses CE; identical labels at '
        'all three thresholds by construction.\n\n')
    p=amin.pivot_table(index='plane',columns='instrument',
        values=['A_min_8_defined','A_min_20_defined','A_min_40_defined'])
    lines.append(p.to_markdown()+'\n')
    hs=hires
    lines.append('\n## Higher-resolution S/T boundary check (81x81, post-principal robustness only)\n\n'
        f"Max nearest-matched crossing shift 41->81: {hs['max_shift']:.4f} "
        f"(tol {hs['tol']:.5f}, within: {hs['within_tol']}); "
        f"A: shift {hs['per_metric']['A']['max_boundary_shift']}, matched {hs['per_metric']['A']['matched_lines']}, "
        f"mismatched {hs['per_metric']['A']['mismatched_lines']}; "
        f"C_social: shift {hs['per_metric']['C_social']['max_boundary_shift']}, matched "
        f"{hs['per_metric']['C_social']['matched_lines']}, mismatched {hs['per_metric']['C_social']['mismatched_lines']}. "
        f"{hs['status']}\n")
    lines.append('\n## Diagnostics and robustness (pointers)\n\n'
        '- Ablations m1-m5 + 6 implicated pairs: results_v12/ablations/ (identity gate PASS; '
        'm1 evasion moves only T-channels, m3 rho=0 only V-channels with a vanishing V-pair '
        'class, m2 collection friction moves all three).\n'
        '- Robustness F1/F2/F4 (49 executed members incl. kappa post-processing and the m1 '
        'reuse): results_v12/robustness/ROBUSTNESS_REPORT.md; strict Sec-8 core nearly empty '
        'on abatement orderings (F2 leverage corners + F1 behavioral flips responsible), '
        'F4 population core large; uniform-nu and nu-visibility recorded infeasible-to-'
        'calibrate under the frozen single-moment rule (nu_refits.json).\n'
        '- F3 governance measurement (9 variants): results_v12/f3/f3_report.md; H2 3/3 and H3 '
        'SUPPORTED under every variant.\n'
        '- Sobol, archetype readouts, V5-prime: respective results_v12 subfolders.\n')
    open(f'{RES}/FINAL_READOUTS.md','w').write(''.join(lines))

if __name__=='__main__':
    amin=a_min_sensitivity(); print('A_min table written')
    hires=hires_check(); print(f"hires: max shift {hires['max_shift']} tol {hires['tol']} within {hires['within_tol']}")
    comps=crossover_components(); print(f'crossover components: {len(comps)} qualifying')
    final_readouts(amin,hires,comps); print('FINAL_READOUTS.md written')
