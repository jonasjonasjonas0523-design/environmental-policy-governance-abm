"""Revision-stage reproduction gate (2026-10-08): verify this workspace reproduces
the frozen v1.2 record exactly before any reviewer-requested new analysis is run."""
import numpy as np, json, os, sys, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))
from model_v12 import CONFIG, batch, run_all, solve_tau_star
from principal_v12 import band_stage, crit_from, labels_and_leaders, pairs

RES='results_v12'
report={}

# ---- Gate A: tau* exact re-solve ----
t0=time.time()
ts = solve_tau_star(10**6, 424242)
frozen = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
report['A_tau_star'] = dict(resolved=ts, frozen=frozen, exact=bool(ts==frozen), secs=round(time.time()-t0,1))
print('GATE A tau*:', 'EXACT' if ts==frozen else f'DIFF {ts!r} vs {frozen!r}')

# ---- Gate D: reference plane + staged family reproduction ----
n=153; seeds=list(range(100000,100000+n)); Rg=41; pp=np.linspace(0,1,Rg)
TAU=frozen; B=CONFIG['B']; KR,HR = CONFIG['ref']['k'], CONFIG['ref']['h']
f_out='revision2026/ref_plane_out.npy'
b = batch(seeds)
if os.path.exists(f_out):
    out=np.load(f_out)
else:
    out=np.zeros((n,Rg*Rg,3,2)); t0=time.time()
    for ic in range(Rg*Rg):
        g=dict(p=pp[ic//Rg], c=pp[ic%Rg], k=KR, h=HR)
        r=run_all(b,g,TAU)
        for ji,j in enumerate('STV'):
            out[:,ic,ji,0]=r[j]['A']; out[:,ic,ji,1]=r[j]['C_social']
        if ic%400==0: print(f'  cell {ic}/{Rg*Rg} ({time.time()-t0:.0f}s)', flush=True)
    np.save(f_out,out); print(f'  plane computed in {time.time()-t0:.0f}s')
W = np.random.default_rng(77).multinomial(n,[1/n]*n,size=B)/n
D,Dhat,mb,se = band_stage(out,W)
crit = crit_from(mb,Dhat,se)
lo = Dhat-crit[None,None,:]*se; hi = Dhat+crit[None,None,:]*se
lo0=np.load(f'{RES}/principal/staged_reference_lo.npy'); hi0=np.load(f'{RES}/principal/staged_reference_hi.npy')
dmax = float(max(np.abs(lo-lo0).max(), np.abs(hi-hi0).max()))
exact = bool(np.array_equal(lo,lo0) and np.array_equal(hi,hi0))
lab,Al,Cl = labels_and_leaders(lo,hi,Rg*Rg)
lab0,Al0,Cl0 = labels_and_leaders(lo0,hi0,Rg*Rg)
report['D_staged_bands'] = dict(exact=exact, max_abs_diff=dmax, crit_A=float(crit[0]), crit_C=float(crit[1]),
                                labels_match=bool((lab==lab0).all()), frozen_critA=3.1101, frozen_critC=3.4758)
print(f'GATE D staged bands: {"EXACT" if exact else f"max|diff|={dmax:.2e}"}; labels match: {(lab==lab0).all()}; crit=({crit[0]:.4f},{crit[1]:.4f})')

json.dump(report, open('revision2026/gate_repro.json','w'), indent=1)
print('GATES A,D DONE')
