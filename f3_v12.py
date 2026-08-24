"""F3 governance-measurement robustness (PRE_REG Sec 8; ARCHETYPE_PLACEMENT_MEMO
Sec 6). Variants: (i) independence copula; (ii) copula entries +/-0.2 (PSD repair
by eigenvalue clipping + unit-diagonal renormalization, documented); (iii) all
marginal means +/-0.10 re-centered, concentration a+b fixed, clipped to
[1e-3, 1-1e-3]; (iv) contemporary-A4 means p=0.30, c=0.40, k=0.30 (h unchanged);
(v) drop A5, recompute G_anchor; plus the +/-20 percent mean stress pair
(reported, not sole treatment). 'Alternative Table 2 bounds' has no operable
alternative specification anywhere in the freeze (Table 2 spans are the locked
[0,1] governance spans); the memo Sec 6 set + stress constitute the executed
family — reading documented here and in the final audit.

Per variant: per-plane anchor masks via the SAME construction as
principal_v12.anchor_masks (rng-2026 stream => CRN: identical standard-normal
draws across variants); archetype readouts via run_all at the variant's first
500 joint draws per archetype (baseline draw-outcome cache reused where an
archetype's marginals AND copula are unchanged: contemporary-A4, drop-A5);
H2 per pair and H3 (per-plane n_pos/n_neg/n_transition rule of the final
crossover summary) re-verdicts on the EXISTING final joint bands.
Gate: baseline anchors (saved plane_*_anchor.npz) must reproduce
H_verdicts_final.json H2/H3 and the final crossover flags before any variant.
Modes: run [budget_s] | report
"""
import numpy as np, pandas as pd, json, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from scipy.stats import norm, beta as beta_dist, gaussian_kde
from model_v12 import CONFIG as CFG, batch, run_all

RES='results_v12'; os.makedirs(f'{RES}/f3', exist_ok=True)
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
n  = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(100000,100000+n))
NDRAW=500; Rg=41; pp=np.linspace(0,1,Rg)
PLANES=[('reference',CFG['ref']['k'],CFG['ref']['h']),('k010',0.10,CFG['ref']['h']),
        ('k055',0.55,CFG['ref']['h']),('h050',CFG['ref']['k'],0.50),('h095',CFG['ref']['k'],0.95)]
pairs=[('S-T',0,1),('S-V',0,2),('T-V',1,2)]
EPS=1e-3

def mshift(ab, d=None, f=None):
    a,b = ab; m = a/(a+b); s = a+b
    m2 = m+d if d is not None else m*f
    m2 = min(max(m2,EPS),1-EPS)
    return (m2*s,(1-m2)*s)

def psd_repair(R):
    R = 0.5*(R+R.T)
    w,V = np.linalg.eigh(R)
    if w.min() < 1e-8:
        w = np.maximum(w,1e-6); R = V@np.diag(w)@V.T
        Dg = np.sqrt(np.diag(R)); R = R/np.outer(Dg,Dg); R = 0.5*(R+R.T)
    return R

def variant_params(tag):
    arch = {k:{d:tuple(v) for d,v in ab.items()} for k,ab in CFG['archetypes'].items()}
    R = np.array(CFG['copula_R'], dtype=float)
    if tag=='independence': R = np.eye(4)
    elif tag in ('cop_p02','cop_m02'):
        d = 0.2 if tag=='cop_p02' else -0.2
        off = ~np.eye(4,dtype=bool)
        R2 = R.copy(); R2[off] = np.clip(R2[off]+d,-0.95,0.95)
        R = psd_repair(R2)
    elif tag in ('means_p010','means_m010'):
        d = 0.10 if tag=='means_p010' else -0.10
        arch = {k:{dm: mshift(ab[dm],d=d) for dm in ab} for k,ab in arch.items()}
    elif tag in ('stress_p20','stress_m20'):
        f = 1.2 if tag=='stress_p20' else 0.8
        arch = {k:{dm: mshift(ab[dm],f=f) for dm in ab} for k,ab in arch.items()}
    elif tag=='cont_A4':
        s = {dm: sum(arch['A4'][dm]) for dm in ('p','c','k')}
        for dm,m2 in (('p',0.30),('c',0.40),('k',0.30)):
            arch['A4'][dm] = (m2*s[dm],(1-m2)*s[dm])
    elif tag=='drop_A5':
        arch.pop('A5')
    else: raise ValueError(tag)
    return arch, R

def gen_masks_draws(arch, R, need_masks=True, only=None):
    """Mirrors principal_v12.anchor_masks rng-2026 consumption exactly."""
    L = np.linalg.cholesky(R)
    rng = np.random.default_rng(2026)
    masks = {pl:{} for pl,_,_ in PLANES}; draws={}
    for name,ab in arch.items():
        z = rng.standard_normal((100000,4))@L.T; u = norm.cdf(z)
        D = np.column_stack([beta_dist.ppf(u[:,i],*ab[dim])
                             for i,dim in enumerate(('p','c','k','h'))])
        sub_idx = rng.choice(100000,3000,replace=False)
        draws[name] = D[:NDRAW]
        if not need_masks or (only is not None and name not in only): continue
        kde = gaussian_kde(D.T)
        thr = np.quantile(kde(D[sub_idx].T),0.05)
        P,C = np.meshgrid(pp,pp,indexing='ij')
        for pl,KR,HR in PLANES:
            pts = np.column_stack([P.ravel(),C.ravel(),np.full(P.size,KR),np.full(P.size,HR)])
            masks[pl][name] = (kde(pts.T)>=thr).reshape(Rg,Rg)
    return masks, draws

_bands_cache=None
def final_bands():
    global _bands_cache
    if _bands_cache is None:
        _bands_cache={}
        for pl,_,_ in PLANES:
            bl = pd.read_parquet(f'{RES}/principal/plane_{pl}_bands_final_joint.parquet')
            d={}
            for pn,_,_ in pairs:
                for m in ('A','C_social'):
                    r = bl[(bl['pair']==pn)&(bl['metric']==m)].sort_values('cell')
                    d[(pn,m)]=(r['lo'].values.reshape(Rg,Rg), r['hi'].values.reshape(Rg,Rg))
            _bands_cache[pl]=d
    return _bands_cache

def verdicts_from_anchors(union_by_plane):
    fb = final_bands()
    H2 = {pn: False for pn,_,_ in pairs}
    xr = []
    for pl,_,_ in PLANES:
        anc = union_by_plane[pl]
        for pn,_,_ in pairs:
            for m in ('A','C_social'):
                lo,hi = fb[pl][(pn,m)]
                pos = (lo>0)&anc; neg = (hi<0)&anc; tr = (lo<=0)&(hi>=0)&anc
                if pos.any() or neg.any(): H2[pn]=True
                xr.append(dict(plane=pl,pair=pn,outcome=m,n_pos=int(pos.sum()),
                               n_neg=int(neg.sum()),n_transition=int(tr.sum()),
                               supported_crossover=bool(pos.any() and neg.any() and tr.any())))
    xdf = pd.DataFrame(xr)
    H3 = bool(xdf.supported_crossover.any())
    return H2, H3, xdf

def readout_rows(OUT, tag):
    rows=[]
    for name,arr in OUT.items():
        for mi,m in enumerate(('A','C_social')):
            vals = arr[:,:,mi]
            best = vals.argmax(1) if m=='A' else vals.argmin(1)
            Pl = [float((best==j).mean()) for j in range(3)]
            modal = int(np.argmax(Pl))
            rows.append(dict(variant=tag, archetype=name, metric=m,
                P_lead_S=Pl[0],P_lead_T=Pl[1],P_lead_V=Pl[2],
                modal_leader='STV'[modal], P_displaced=1-Pl[modal]))
    return rows

VTAGS=['independence','cop_p02','cop_m02','means_p010','means_m010',
       'cont_A4','drop_A5','stress_p20','stress_m20']
REUSE_BASE={'cont_A4':{'A1','A2','A3','A5'}, 'drop_A5':{'A1','A2','A3','A4'}}

def gate():
    gf=f'{RES}/f3/gate.json'
    if os.path.exists(gf): return
    union={}
    for pl,_,_ in PLANES:
        z=np.load(f'{RES}/principal/plane_{pl}_anchor.npz')
        union[pl]=z['union'].astype(bool)
    H2,H3,xdf = verdicts_from_anchors(union)
    ref=json.load(open(f'{RES}/principal/H_verdicts_final.json'))
    xref=pd.read_csv(f'{RES}/principal/crossover_summary_final_all_planes.csv')
    mrg=xdf.merge(xref,on=['plane','pair','outcome'],suffixes=('','_ref'))
    ok_x=bool((mrg.supported_crossover==mrg.supported_crossover_ref).all())
    h2ref={k:bool(v) for k,v in ref['H2'].items()} if isinstance(ref.get('H2'),dict) else None
    ok_h2 = (H2==h2ref) if h2ref is not None else None
    ok_h3 = (H3==bool(ref.get('H3_supported', ref.get('H3'))))
    ok = ok_x and (ok_h2 in (True,None)) and ok_h3
    json.dump(dict(gate='F3 baseline reproduction from saved anchors',
                   crossover_flags_match=ok_x, H2_match=ok_h2, H3_match=ok_h3,
                   status='PASS' if ok else 'FAIL'), open(gf,'w'), indent=1)
    print(f'F3 GATE baseline reproduction: {"PASS" if ok else "FAIL"} '
          f'(xflags {ok_x}, H2 {ok_h2}, H3 {ok_h3})')
    assert ok

def run_variant(tag, b, OUT_base):
    arch,R = variant_params(tag)
    only = None
    if tag in REUSE_BASE:
        changed = set(arch)-REUSE_BASE[tag]
        only = changed
    masks,draws = gen_masks_draws(arch,R,only=only)
    union={}
    for pl,_,_ in PLANES:
        if tag in REUSE_BASE:
            z=np.load(f'{RES}/principal/plane_{pl}_anchor.npz',allow_pickle=True)
            mdict=z['masks'].item()
            for nm in REUSE_BASE[tag]&set(arch):
                masks[pl][nm]=np.asarray(mdict[nm],bool)
        u=np.zeros((Rg,Rg),bool)
        for m_ in masks[pl].values(): u|=m_
        union[pl]=u
    OUT={}
    for name,D in draws.items():
        if tag in REUSE_BASE and name in REUSE_BASE[tag] and name in OUT_base:
            OUT[name]=OUT_base[name]; continue
        arr=np.zeros((NDRAW,3,2))
        for di,(p,c,k,h) in enumerate(D):
            r=run_all(b,dict(p=float(p),c=float(c),k=float(k),h=float(h)),TAU)
            for ji,j in enumerate('STV'):
                arr[di,ji,0]=r[j]['A'].mean(); arr[di,ji,1]=r[j]['C_social'].mean()
        OUT[name]=arr
    H2,H3,xdf=verdicts_from_anchors(union)
    np.savez(f'{RES}/f3/{tag}.npz',
             **{f'union_{pl}':union[pl] for pl,_,_ in PLANES},
             OUT=np.array(OUT,dtype=object), H2=np.array(json.dumps(H2)),
             H3=np.array(H3), anch_frac=np.array([float(union[pl].mean()) for pl,_,_ in PLANES]))
    xdf.to_csv(f'{RES}/f3/{tag}_crossover.csv',index=False)
    pd.DataFrame(readout_rows(OUT,tag)).to_csv(f'{RES}/f3/{tag}_readouts.csv',index=False)
    return H2,H3,union

def mode_run(budget=240):
    t0=time.time(); gate()
    b=batch(seeds)
    OUT_base=np.load(f'{RES}/archetypes/draw_outcomes.npz',allow_pickle=True)['OUT'].item()
    for tag in VTAGS:
        if os.path.exists(f'{RES}/f3/{tag}.npz'): continue
        if time.time()-t0>budget:
            print('BUDGET pause; remaining:',[t for t in VTAGS if not os.path.exists(f'{RES}/f3/{t}.npz')])
            return
        ts=time.time(); H2,H3,union=run_variant(tag,b,OUT_base)
        print(f'  {tag}: {time.time()-ts:.0f}s  H2={sum(H2.values())}/3  H3={H3}  '
              f'anch_ref={union["reference"].mean():.3f}', flush=True)
    print('F3 DONE')

def mode_report():
    ref=json.load(open(f'{RES}/principal/H_verdicts_final.json'))
    base_ro=pd.read_csv(f'{RES}/archetypes/archetype_readouts.csv')
    rows=[]; ver={}
    for tag in VTAGS:
        z=np.load(f'{RES}/f3/{tag}.npz',allow_pickle=True)
        H2=json.loads(str(z['H2'])); H3=bool(z['H3'])
        ver[tag]=dict(H2=H2,H3_supported=H3,
                      anchored_fraction={pl:float(f) for (pl,_,_),f in zip(PLANES,z['anch_frac'])})
        ro=pd.read_csv(f'{RES}/f3/{tag}_readouts.csv'); rows.append(ro)
    allro=pd.concat(rows)
    json.dump(dict(baseline=dict(H2=ref.get('H2_detail'),H3_supported=ref.get('H3_supported')),
                   variants=ver,
                   note=('Alternative-Table-2-bounds member: no operable alternative specification '
                         'exists in the freeze; memo Sec 6 set + stress pair constitute the executed '
                         'family (documented in FINAL_COMPUTATIONAL_AUDIT).')),
              open(f'{RES}/f3/f3_verdicts.json','w'),indent=1)
    with open(f'{RES}/f3/f3_report.md','w') as f:
        f.write('# F3 governance-measurement robustness (memo Sec 6 set + stress pair)\n\n'
                'Anchors and archetype readouts recomputed per variant on the EXISTING final joint '
                'bands and the frozen seed-2026 draw stream (CRN across variants); H2/H3 re-verdicts '
                'per the memo ("readouts and the H2/H3 verdicts are re-reported under each").\n\n')
        f.write('| variant | H2 S-T/S-V/T-V | H3 | anchored frac (ref plane) |\n|---|---|---|---|\n')
        for tag in VTAGS:
            v=ver[tag]; h2=v['H2']
            f.write(f"| {tag} | {int(h2['S-T'])}/{int(h2['S-V'])}/{int(h2['T-V'])} | "
                    f"{'SUPPORTED' if v['H3_supported'] else 'not supported'} | "
                    f"{v['anchored_fraction']['reference']:.3f} |\n")
        f.write('\n## Modal-leader stability (A metric)\n\n')
        piv=allro[allro.metric=='A'].pivot_table(index='archetype',columns='variant',
              values='modal_leader',aggfunc='first')
        base_m=base_ro[base_ro.metric=='A'].set_index('archetype').modal_leader
        piv.insert(0,'baseline',base_m)
        f.write(piv.to_markdown()+'\n\n## Modal-leader stability (C_social)\n\n')
        piv2=allro[allro.metric=='C_social'].pivot_table(index='archetype',columns='variant',
              values='modal_leader',aggfunc='first')
        base_c=base_ro[base_ro.metric=='C_social'].set_index('archetype').modal_leader
        piv2.insert(0,'baseline',base_c)
        f.write(piv2.to_markdown()+'\n')
    print('f3 report written')

if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'run'
    if mode=='run': mode_run(int(sys.argv[2]) if len(sys.argv)>2 else 240)
    else: mode_report()
