"""Robustness families F1/F2/F4 driver (PRE_REG Sec 8; ranges PARAMETER_PROVENANCE;
v1.2 members AMENDMENTS Sec 3.2). Each simulated variant re-runs the reference
plane (41x41, principal seeds 100000..100152) under run_all_R / run_all with the
variant's cfg and decision rules, then gets its own per-plane joint max-t family
(rng-77 weights, identical construction to the staged reference analysis) and
classification. Baseline for all comparisons = staged reference plane (same
convention as the ablations; per-plane family vs per-plane family). kappa
variants are exact post-processing of the stored seed-level reference outputs
(C_social = C_ab + kappa'.counts; A unchanged). evasion-off REUSES ablation m1
(identical computation, same seeds). Robust core = per pair|metric intersection
of supported classifications across baseline and ALL variants (Sec VIII-A);
fragile features attributed to the families that break them.
nu refits (frozen single-moment rule): nu_0 solved so EVER_all = 0.15 alone at
the adopted rung p=0.29 with eta fixed (0 uniform; 1 visibility), on the frozen
seed-20260821 400k stream via the calibration script's own ever_stats.
Modes: run [budget_s] (self-checkpointing) | core
"""
import numpy as np, json, os, sys, time
import pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from model_v12 import CONFIG, batch, run_all
from model_v12_robust import run_all_R
from ablations_v12 import bands, labels_and_signs, pairs

RES='results_v12'; CFG=CONFIG
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
os.makedirs(f'{RES}/robustness', exist_ok=True)
n  = json.load(open('results/pilot/replication_choice.json'))['chosen_n']
seeds = list(range(100000,100000+n))
Rg=41; pp=np.linspace(0,1,Rg); KR,HR = CFG['ref']['k'], CFG['ref']['h']
B = CFG['B']

def refit_nu0(eta):
    import importlib.util
    spec = importlib.util.spec_from_file_location('cal','freeze_v12/v1_2_calibration_script.py')
    cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
    lo,hi = 0.0,3.0
    for _ in range(48):
        mid = 0.5*(lo+hi)
        sa,_,_ = cal.ever_stats(mid, eta, 0.29, cal.gF, cal.eF, cal.topF)
        if sa < 0.15: lo = mid
        else: hi = mid
    nu0 = 0.5*(lo+hi)
    sa,st,_ = cal.ever_stats(nu0, eta, 0.29, cal.gF, cal.eF, cal.topF)
    return nu0, sa, st

def get_refits():
    f = f'{RES}/robustness/nu_refits.json'
    if os.path.exists(f):
        d = json.load(open(f))
        if 'feasible' in d.get('uniform',{}): return d
    import importlib.util
    spec = importlib.util.spec_from_file_location('cal','freeze_v12/v1_2_calibration_script.py')
    cal = importlib.util.module_from_spec(spec); spec.loader.exec_module(cal)
    out={}
    for tag,eta in (('uniform',0.0),('visibility',1.0)):
        sa0,_,_ = cal.ever_stats(0.0, eta, 0.29, cal.gF, cal.eF, cal.topF)
        if sa0 >= 0.15 - 0.03:   # EVER nondecreasing in nu0 => no admissible solution
            out[tag]=dict(eta=eta, feasible=False, EVER_all_at_nu0_zero=sa0,
                target=0.15, tol=0.03,
                mechanism=('threat-channel ignition at the adopted rung (V1_2_CALIBRATION_MEMO '
                 'Sec 4): with eta forced to this value the participant set is no longer '
                 'concentrated in the largest emitters, recognized shortfall stays high, '
                 'theta reaches its ceiling, and EVER_all exceeds the 0.15 target for ALL '
                 'nu_0 >= 0. Member recorded infeasible-to-calibrate under the frozen '
                 'single-moment rule; a boundary nu_0=0 run would duplicate ablation m3.'))
        else:
            nu0,sa,st = refit_nu0(eta)
            out[tag]=dict(eta=eta, feasible=True, nu0=nu0, EVER_all=sa, EVER_top=st)
    out['rule']=('frozen single-moment rule: eta fixed, nu_0 bisected on [0,3] so '
                 'EVER_all = 0.15 alone at adopted rung p=0.29, 400k seed-20260821 '
                 'stream, calibration-script ever_stats verbatim')
    json.dump(out, open(f,'w'), indent=1); return out

def batch_corr(seeds_, cfg, r):
    ds=[]
    for s in seeds_:
        rng = np.random.default_rng(int(s)); N,T = cfg['N'],cfg['T']
        zg = rng.normal(0.0,1.0,N); zE = rng.normal(0.0,1.0,N)
        ds.append(dict(gamma=np.exp(cfg['mu_g']+cfg['sg']*zg),
                       E0=np.exp(cfg['mu_E']+cfg['sE']*(r*zg+np.sqrt(1-r*r)*zE)),
                       U_mon=rng.random((T,N)), U_imp=rng.random((T,N)),
                       U_imp2=rng.random((T,N))))
    return {k: np.stack([d[k] for d in ds]) for k in ds[0]}

def L(x): return float(np.log(x/(1-x)))

def registry():
    rf = get_refits()
    Rv = {}
    # ---- F1 behavioral ----
    Rv['cara_r1']  = dict(fam='F1', v=dict(rule='cara', r=1.0))
    Rv['cara_r10'] = dict(fam='F1', v=dict(rule='cara', r=10.0))
    Rv['satisficing'] = dict(fam='F1', v=dict(rule='satisficing'))
    Rv['normative'] = dict(fam='F1', v=dict(omega=0.25*TAU))
    Rv['power15'] = dict(fam='F1', v=dict(cost='power', phi=1.5))
    Rv['power3']  = dict(fam='F1', v=dict(cost='power', phi=3.0))
    Rv['pwl_cost'] = dict(fam='F1', v=dict(cost='pwl'))
    Rv['convex_evasion'] = dict(fam='F1', v=dict(convex_evasion=True))
    # evasion_off: alias of ablation m1 (reuse; identical computation/seeds)
    # ---- F2 instrument design ----
    Rv['mS1'] = dict(fam='F2', cfg=dict(m_S=1.0)); Rv['mS4'] = dict(fam='F2', cfg=dict(m_S=4.0))
    Rv['mT05'] = dict(fam='F2', cfg=dict(m_T=0.5)); Rv['mT2'] = dict(fam='F2', cfg=dict(m_T=2.0))
    Rv['Astar010'] = dict(fam='F2', cfg=dict(A_star_frac=0.10, A_cum_nominal=200.0, A_min=10.0))
    Rv['Astar030'] = dict(fam='F2', cfg=dict(A_star_frac=0.30, A_cum_nominal=600.0, A_min=30.0))
    Rv['T5'] = dict(fam='F2', cfg=dict(T=5), tbatch=True)
    Rv['T20'] = dict(fam='F2', cfg=dict(T=20), tbatch=True)
    Rv['theta0_lo'] = dict(fam='F2', cfg=dict(theta0=L(0.05)))
    Rv['theta0_hi'] = dict(fam='F2', cfg=dict(theta0=L(0.15)))
    Rv['lam15'] = dict(fam='F2', cfg=dict(lam=1.5)); Rv['lam45'] = dict(fam='F2', cfg=dict(lam=4.5))
    Rv['probit_theta'] = dict(fam='F2', v=dict(theta_form='probit'))
    Rv['linear_theta'] = dict(fam='F2', v=dict(theta_form='linear'))
    Rv['rep_standard'] = dict(fam='F2', v=dict(rule='rep_standard'))
    Rv['thr_uniform'] = dict(fam='F2', v=dict(threshold='uniform'))
    Rv['thr_cost'] = dict(fam='F2', v=dict(threshold='cost'))
    Rv['tgt_uniform'] = dict(fam='F2', v=dict(target='uniform'))
    Rv['tgt_cost'] = dict(fam='F2', v=dict(target='cost'))
    if rf['visibility'].get('feasible'):
        Rv['nu_visibility'] = dict(fam='F2', v=dict(nu0=rf['visibility']['nu0'], eta=1.0))
    Rv['one_period_M'] = dict(fam='F2', v=dict(one_period_M=True))
    if rf['uniform'].get('feasible'):
        Rv['uniform_nu'] = dict(fam='F2', v=dict(nu0=rf['uniform']['nu0'], eta=0.0))
    Rv['top_decile'] = dict(fam='F2', v=dict(nu0=0.2920, eta=2.0000))
    Rv['tail_cap'] = dict(fam='F2', v=dict(nu_cap=2.0))
    Rv['collective'] = dict(fam='F2', v=dict(collective=True))
    Rv['within_tick'] = dict(fam='F2', v=dict(within_tick=True))
    # ---- F4 population structure ----
    Rv['N50'] = dict(fam='F4', cfg=dict(N=50), tbatch=True)
    Rv['N1000'] = dict(fam='F4', cfg=dict(N=1000), tbatch=True)
    Rv['sg05'] = dict(fam='F4', cfg=dict(sg=0.5, mu_g=-0.125), tbatch=True)
    Rv['sg15'] = dict(fam='F4', cfg=dict(sg=1.5, mu_g=-1.125), tbatch=True)
    Rv['sE05'] = dict(fam='F4', cfg=dict(sE=0.5, mu_E=-0.125), tbatch=True)
    Rv['sE15'] = dict(fam='F4', cfg=dict(sE=1.5, mu_E=-1.125), tbatch=True)
    Rv['corr03'] = dict(fam='F4', corr=-0.3, tbatch=True)
    Rv['delta090'] = dict(fam='F4', cfg=dict(delta=0.90))
    Rv['delta097'] = dict(fam='F4', cfg=dict(delta=0.97))
    return Rv

_base_batch = None
def get_batch(spec, cfg):
    global _base_batch
    if 'corr' in spec: return batch_corr(seeds, cfg, spec['corr'])
    if spec.get('tbatch'): return batch(seeds, cfg)
    if _base_batch is None: _base_batch = batch(seeds)
    return _base_batch

def run_variant(tag, spec, W):
    cfg = dict(CFG); cfg.update(spec.get('cfg', {}))
    v = spec.get('v', None)
    b = get_batch(spec, cfg)
    out = np.zeros((n, Rg*Rg, 3, 2))
    t0=time.time()
    for ic in range(Rg*Rg):
        g = dict(p=pp[ic//Rg], c=pp[ic%Rg], k=KR, h=HR)
        r = run_all_R(b,g,TAU,cfg,v) if v else run_all(b,g,TAU,cfg)
        for ji,j in enumerate('STV'):
            out[:,ic,ji,0]=r[j]['A']; out[:,ic,ji,1]=r[j]['C_social']
    dt = time.time()-t0
    Dh,se,lo,hi,crit = bands(out, W)
    lab,s = labels_and_signs(lo,hi)
    np.savez(f'{RES}/robustness/{tag}_bands.npz', Dhat=Dh,se=se,lo=lo,hi=hi,
             crit=crit, lab=lab, signs=s, seconds=dt)
    np.save(f'{RES}/robustness/{tag}_out.npy', out)
    return dt, crit

def mode_run(budget=250):
    t_start = time.time()
    Rv = registry()
    W = np.random.default_rng(77).multinomial(n,[1/n]*n,size=B)/n
    # baseline-identity gate for run_all_R (once)
    gate_f = f'{RES}/robustness/identity_gate.json'
    if not os.path.exists(gate_f):
        b = get_batch({}, CFG)
        g0 = dict(p=0.5,c=0.5,k=KR,h=HR)
        r1 = run_all(b,g0,TAU); r2 = run_all_R(b,g0,TAU,CFG,{})
        ok = all(np.array_equal(r1[j][m], r2[j][m])
                 for j in 'STV' for m in ('A','C_social','nM','nI'))
        json.dump(dict(gate='run_all_R identity at v={}', status='PASS' if ok else 'FAIL'),
                  open(gate_f,'w'), indent=1)
        print(f'GATE run_all_R identity (v empty): {"PASS" if ok else "FAIL"}'); assert ok
    done=[]
    for tag,spec in Rv.items():
        if os.path.exists(f'{RES}/robustness/{tag}_bands.npz'): continue
        if time.time()-t_start > budget:
            print(f'BUDGET pause; remaining: '
                  f'{[t for t in Rv if not os.path.exists(f"{RES}/robustness/{t}_bands.npz")]}')
            return
        dt,crit = run_variant(tag, spec, W)
        done.append(tag)
        print(f'  {tag} [{spec["fam"]}]: {dt:.0f}s critA={crit[0]:.2f} critC={crit[1]:.2f}', flush=True)
    print('ALL SIM VARIANTS DONE')

def kappa_variants(W):
    df = pd.read_parquet(f'{RES}/principal/principal_raw_results.parquet')
    df = df.sort_values(['seed','p','c','instrument'], kind='mergesort')
    order = {'S':0,'T':1,'V':2}
    ii = df.instrument.map(order).to_numpy()
    df = df.assign(_i=ii).sort_values(['seed','p','c','_i'], kind='mergesort')
    arr = lambda col: df[col].to_numpy(np.float64).reshape(n, Rg*Rg, 3)
    A, C_ab = arr('A'), arr('C_ab')
    nM, nA, nI = arr('nM'), arr('nA'), arr('nI')
    base = np.load(f'{RES}/principal/plane_reference_out.npy')
    assert np.allclose(A, base[:,:,:,0], rtol=0, atol=2e-2), 'parquet/out alignment'
    outv={}
    kM,kA,kI = CFG['kappa_M'],CFG['kappa_A'],CFG['kappa_I']
    variants = dict(kappa_x5=(kM*5,kA*5,kI*5), kappa_d5=(kM/5,kA/5,kI/5),
                    kappaM_x5=(kM*5,kA,kI), kappaM_d5=(kM/5,kA,kI),
                    kappaA_x5=(kM,kA*5,kI), kappaA_d5=(kM,kA/5,kI),
                    kappaI_x5=(kM,kA,kI*5), kappaI_d5=(kM,kA,kI/5))
    for tag,(m_,a_,i_) in variants.items():
        C = C_ab + m_*nM + a_*nA + i_*nI
        out = np.stack([base[:,:,:,0], C], axis=3)
        Dh,se,lo,hi,crit = bands(out, W)
        lab,s = labels_and_signs(lo,hi)
        np.savez(f'{RES}/robustness/{tag}_bands.npz', Dhat=Dh,se=se,lo=lo,hi=hi,
                 crit=crit, lab=lab, signs=s, seconds=0.0,
                 note='post-processed from stored seed-level reference outputs (float32 fields)')
        outv[tag]=(float(crit[0]),float(crit[1]))
    return outv

def mode_core():
    Rv = registry()
    W = np.random.default_rng(77).multinomial(n,[1/n]*n,size=B)/n
    kv = kappa_variants(W)
    print('kappa variants:', {k:tuple(round(x,2) for x in v) for k,v in kv.items()})
    lo0 = np.load(f'{RES}/principal/staged_reference_lo.npy')
    hi0 = np.load(f'{RES}/principal/staged_reference_hi.npy')
    lab0, s0 = labels_and_signs(lo0, hi0)
    fams = {t:sp['fam'] for t,sp in Rv.items()}
    fams.update({t:'F4' for t in ('kappa_x5','kappa_d5','kappaM_x5','kappaM_d5',
                                  'kappaA_x5','kappaA_d5','kappaI_x5','kappaI_d5')})
    fams['evasion_off']='F1'
    all_tags = list(fams)
    S_stack={}; L_stack={}
    rows=[]; budget=[]
    for t in all_tags:
        f = (f'{RES}/ablations/m1_bands.npz' if t=='evasion_off'
             else f'{RES}/robustness/{t}_bands.npz')
        z = np.load(f)
        S_stack[t]=z['signs']; L_stack[t]=z['lab']
        agree = float((z['signs']==s0).mean())
        rows.append(dict(variant=t, family=fams[t],
                         crit_A=float(z['crit'][0]), crit_C=float(z['crit'][1]),
                         sign_agreement=agree,
                         supported_cells=int((z['signs']!=0).sum())))
        budget.append(dict(variant=t, family=fams[t],
                           seconds=float(z['seconds']) if 'seconds' in z.files else None,
                           reused='m1' if t=='evasion_off' else None))
    Sv = np.stack([S_stack[t] for t in all_tags])          # [V,cells,3,2]
    core_pos = (s0==1) & (Sv==1).all(0)
    core_neg = (s0==-1) & (Sv==-1).all(0)
    core = np.where(core_pos,1,np.where(core_neg,-1,0))
    Lv = np.stack([L_stack[t] for t in all_tags])
    core_lab = np.where((lab0>0) & (Lv==lab0[None]).all(0), lab0, 0)
    frag = (s0!=0) & (core==0)
    breakers = {}
    for pi,(pn,_,_) in enumerate(pairs):
        for mi,m in enumerate(('A','C_social')):
            mask = frag[:,pi,mi]
            if not mask.any():
                breakers[f'{pn}|{m}'] = {}
                continue
            cnt={}
            for vi,t in enumerate(all_tags):
                nb = int(((Sv[vi,:,pi,mi]!=s0[:,pi,mi]) & mask).sum())
                if nb: cnt[t]=nb
            breakers[f'{pn}|{m}'] = dict(sorted(cnt.items(), key=lambda x:-x[1])[:6])
    np.savez(f'{RES}/robustness/robust_core.npz', core=core, core_lab=core_lab,
             baseline_signs=s0, baseline_lab=lab0, tags=np.array(all_tags))
    df = pd.DataFrame(rows).sort_values(['family','variant'])
    df.to_csv(f'{RES}/robustness/robustness_summary.csv', index=False)
    fam_budget={}
    for r in budget:
        fam_budget.setdefault(r['family'], dict(variants=0, seconds=0.0, reused=[]))
        fam_budget[r['family']]['variants'] += 1
        if r['seconds']: fam_budget[r['family']]['seconds'] += r['seconds']
        if r['reused']: fam_budget[r['family']]['reused'].append(r['variant'])
    fam_budget['F3']=dict(variants=9, seconds=None,
                          note='anchor/readout recomputation on existing surfaces (f3_v12.py)')
    json.dump(dict(per_variant=budget, per_family=fam_budget,
                   plane='reference 41x41, n=153 principal seeds, per-plane joint bands rng-77'),
              open(f'{RES}/robustness/budget_log.json','w'), indent=1, default=float)
    with open(f'{RES}/robustness/ROBUSTNESS_REPORT.md','w') as f:
        f.write('# Robustness families F1/F2/F4: robust core on the reference plane\n\n')
        f.write(f'Variants evaluated: {len(all_tags)} (+9 F3 measurement variants '
                'reported separately in f3_report.md). Baseline = staged reference-plane '
                'classification. Robust core = supported sign preserved across baseline '
                'and ALL variants (PRE_REG Sec 8 intersection rule).\n\n')
        for pi,(pn,_,_) in enumerate(pairs):
            for mi,m in enumerate(('A','C_social')):
                b0 = int((s0[:,pi,mi]!=0).sum()); bc = int((core[:,pi,mi]!=0).sum())
                f.write(f'- {pn} | {m}: baseline supported {b0} cells -> robust core {bc} '
                        f'({bc/max(b0,1):.0%} retained); top breakers: '
                        f'{breakers[f"{pn}|{m}"]}\n')
        nl0 = int((lab0>0).sum()); nlc = int((core_lab>0).sum())
        f.write(f'- label map: baseline classified {nl0} -> identical under all variants {nlc}\n\n')
        f.write(df.to_markdown(index=False))
    print('core written; per-pair retention:',
          {f'{pn}|{m}': f'{int((core[:,pi,mi]!=0).sum())}/{int((s0[:,pi,mi]!=0).sum())}'
           for pi,(pn,_,_) in enumerate(pairs) for mi,m in enumerate(("A","C_social"))})

if __name__=='__main__':
    mode = sys.argv[1] if len(sys.argv)>1 else 'run'
    if mode=='run': mode_run(int(sys.argv[2]) if len(sys.argv)>2 else 250)
    elif mode=='core': mode_core()
