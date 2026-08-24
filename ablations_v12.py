"""Mechanism ablations m1-m5 on the reference plane (PRE_REG Sec 7), v1.2.
Same seeds/grid as the principal reference plane. Each ablation gets its own
per-plane joint max-t family (rng 77 weights, identical construction to the
staged reference analysis) and classification; comparisons are therefore
per-plane-family vs per-plane-family (baseline = staged reference plane).
  m1 evasion off (Etilde=E)                 -- run_charge truthful branch forced
  m2 collection friction off (c=1,k=0 in collection roles; evasion untouched)
  m3 rho = 0                                -- cfg nu0=0 (nu_i == 0)
  m4 theta -> 0                             -- cfg theta0=-700, lam=0
  m5 k = 0 globally                         -- cell k overridden to 0
Outputs per ablation: mechanism-contribution surfaces Dm = Delta_base - Delta_abl
(3 pairs x 2 metrics), label map, per-pair sign maps, band records. Trigger:
if any single ablation moves a supported classification boundary by more than
one grid cell (Hausdorff distance between boundary-cell sets, Chebyshev metric,
over the label map and each pair|metric sign map; an appearing/vanishing
boundary counts as infinite shift), the six paired ablations among
{m1,m2,m4,m5} are run with the same treatment (predeclared rule).
Gate first: with no switches active, the ablation module must reproduce
model_v12.run_all bitwise on a probe cell.
"""
import numpy as np, json, os, sys, time
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch, run_all
from model_v12_ablations import run_all_abl

RES='results_v12'; CFG=CONFIG
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
os.makedirs(f'{RES}/ablations', exist_ok=True)
n  = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(100000,100000+n))
Rg = 41; pp = np.linspace(0,1,Rg)
pairs = [('S-T',0,1),('S-V',0,2),('T-V',1,2)]
B = CFG['B']; KR,HR = CFG['ref']['k'], CFG['ref']['h']

ABL = dict(
    m1=dict(kind='flag', abl=dict(no_evasion=True)),
    m2=dict(kind='flag', abl=dict(no_collection_friction=True)),
    m3=dict(kind='cfg',  cfg=dict(nu0=0.0)),
    m4=dict(kind='cfg',  cfg=dict(theta0=-700.0, lam=0.0)),
    m5=dict(kind='kzero'),
)
PAIRED = [('m1','m2'),('m1','m4'),('m1','m5'),('m2','m4'),('m2','m5'),('m4','m5')]

def merged(names):
    abl={}; cfg=dict(CFG); kzero=False
    for nm in names:
        s=ABL[nm]
        if s['kind']=='flag': abl.update(s['abl'])
        elif s['kind']=='cfg': cfg.update(s['cfg'])
        elif s['kind']=='kzero': kzero=True
    return abl, cfg, kzero

def plane_out(tag, names, b):
    f = f'{RES}/ablations/{tag}_out.npy'
    if os.path.exists(f): return np.load(f)
    abl,cfg,kzero = merged(names)
    out = np.zeros((n,Rg*Rg,3,2))
    t0=time.time()
    for ic in range(Rg*Rg):
        p,c = pp[ic//Rg], pp[ic%Rg]
        g = dict(p=p,c=c,k=(0.0 if kzero else KR),h=HR)
        r = run_all_abl(b,g,TAU,cfg,**abl) if abl else run_all(b,g,TAU,cfg)
        for ji,j in enumerate('STV'):
            out[:,ic,ji,0]=r[j]['A']; out[:,ic,ji,1]=r[j]['C_social']
    np.save(f,out); print(f'  {tag}: plane in {time.time()-t0:.0f}s', flush=True)
    return out

def bands(out, W):
    D = np.stack([out[:,:,i,:]-out[:,:,j,:] for _,i,j in pairs],2)
    Dhat = D.mean(0); mb=(W@D.reshape(n,-1)).reshape(B,Rg*Rg,3,2); se=mb.std(0,ddof=1)
    crit=np.zeros(2)
    for mi in range(2):
        s=np.where(se[:,:,mi]>0,se[:,:,mi],np.inf)
        t=np.abs((mb[:,:,:,mi]-Dhat[None,:,:,mi])/s[None]).reshape(B,-1).max(1)
        crit[mi]=float(np.quantile(t,0.95))
    return Dhat, se, Dhat-crit[None,None,:]*se, Dhat+crit[None,None,:]*se, crit

def labels_and_signs(lo,hi):
    from principal_v12 import labels_and_leaders
    lab,Aled,Cled = labels_and_leaders(lo,hi,Rg*Rg)
    signs = np.zeros((Rg*Rg,3,2),dtype=int)   # +1 pos, -1 neg, 0 transition
    signs[lo>0]=1; signs[hi<0]=-1
    return lab, signs

def boundary_cells(m2d):
    bset=set()
    for i in range(Rg):
        for j in range(Rg):
            v=m2d[i,j]
            if (i+1<Rg and m2d[i+1,j]!=v) or (j+1<Rg and m2d[i,j+1]!=v):
                bset.add((i,j))
                if i+1<Rg and m2d[i+1,j]!=v: bset.add((i+1,j))
                if j+1<Rg and m2d[i,j+1]!=v: bset.add((i,j+1))
    return bset

def hausdorff(b1,b2):
    if not b1 and not b2: return 0.0
    if not b1 or not b2: return np.inf
    a1=np.array(sorted(b1)); a2=np.array(sorted(b2))
    d12 = max(np.abs(a1[:,None,:]-a2[None,:,:]).max(2).min(1).max(), 0)
    d21 = max(np.abs(a2[:,None,:]-a1[None,:,:]).max(2).min(1).max(), 0)
    return float(max(d12,d21))

def shift_report(lab0,s0,lab1,s1):
    fams={'label_map':(lab0.reshape(Rg,Rg),lab1.reshape(Rg,Rg))}
    for pi,(pn,_,_) in enumerate(pairs):
        for mi,m in enumerate(('A','C_social')):
            fams[f'{pn}|{m}']=(s0[:,pi,mi].reshape(Rg,Rg),s1[:,pi,mi].reshape(Rg,Rg))
    return {k:hausdorff(boundary_cells(a),boundary_cells(bb)) for k,(a,bb) in fams.items()}

def main():
    b = batch(seeds)
    # gate: no-switch ablation module == frozen model bitwise on probe cell
    g0=dict(p=0.5,c=0.5,k=KR,h=HR)
    r1=run_all(b,g0,TAU); r2=run_all_abl(b,g0,TAU)
    ok=all(np.array_equal(r1[j][m],r2[j][m]) for j in 'STV' for m in ('A','C_social','nM','nI'))
    print(f'GATE ablation-module identity (no switches): {"PASS" if ok else "FAIL"}'); assert ok
    W = np.random.default_rng(77).multinomial(n,[1/n]*n,size=B)/n
    base = np.load(f'{RES}/principal/plane_reference_out.npy')
    lo0=np.load(f'{RES}/principal/staged_reference_lo.npy'); hi0=np.load(f'{RES}/principal/staged_reference_hi.npy')
    Dhat0 = np.stack([base[:,:,i,:]-base[:,:,j,:] for _,i,j in pairs],2).mean(0)
    lab0,s0 = labels_and_signs(lo0,hi0)
    rows=[]; trigger=False; shifts_all={}
    for tag in ABL:
        out = plane_out(tag,[tag],b)
        Dh,se,lo,hi,crit = bands(out,W)
        lab,s = labels_and_signs(lo,hi)
        np.savez(f'{RES}/ablations/{tag}_bands.npz',Dhat=Dh,se=se,lo=lo,hi=hi,crit=crit,lab=lab,signs=s)
        Dm = Dhat0-Dh
        np.save(f'{RES}/ablations/{tag}_Dm.npy',Dm)
        sh = shift_report(lab0,s0,lab,s); shifts_all[tag]=sh
        mx = max(sh.values()); trigger |= (mx>1)
        for pi,(pn,_,_) in enumerate(pairs):
            for mi,m in enumerate(('A','C_social')):
                rows.append(dict(ablation=tag,pair=pn,metric=m,
                    mean_Dm=float(Dm[:,pi,mi].mean()),max_abs_Dm=float(np.abs(Dm[:,pi,mi]).max()),
                    boundary_shift=sh[f'{pn}|{m}']))
        print(f'  {tag}: crit_A={crit[0]:.3f} crit_C={crit[1]:.3f} max boundary shift={mx if np.isfinite(mx) else "inf"}', flush=True)
    json.dump(dict(single_shifts=shifts_all,
                   trigger_paired=bool(trigger),
                   rule='paired ablations among {m1,m2,m4,m5} if any supported classification boundary moves >1 grid cell (Hausdorff, Chebyshev)'),
              open(f'{RES}/ablations/trigger.json','w'),indent=1,default=float)
    print(f'PAIRED TRIGGER: {trigger}')
    if trigger:
        for a1,a2 in PAIRED:
            tag=f'{a1}_{a2}'
            out=plane_out(tag,[a1,a2],b)
            Dh,se,lo,hi,crit=bands(out,W)
            lab,s=labels_and_signs(lo,hi)
            np.savez(f'{RES}/ablations/{tag}_bands.npz',Dhat=Dh,se=se,lo=lo,hi=hi,crit=crit,lab=lab,signs=s)
            np.save(f'{RES}/ablations/{tag}_Dm.npy',Dhat0-Dh)
            sh=shift_report(lab0,s0,lab,s)
            for pi,(pn,_,_) in enumerate(pairs):
                for mi,m in enumerate(('A','C_social')):
                    rows.append(dict(ablation=tag,pair=pn,metric=m,
                        mean_Dm=float((Dhat0-Dh)[:,pi,mi].mean()),
                        max_abs_Dm=float(np.abs((Dhat0-Dh)[:,pi,mi]).max()),
                        boundary_shift=sh[f'{pn}|{m}']))
            print(f'  paired {tag}: done (max shift {max(sh.values())})', flush=True)
    pd.DataFrame(rows).to_csv(f'{RES}/ablations/ablation_summary.csv',index=False)
    print('ablations written')

if __name__=='__main__':
    main()
