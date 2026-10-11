"""Shared helpers for the 2026-10 revision-stage analyses (reviewer-requested;
post-hoc, labeled; frozen baseline untouched). Convention for variant comparisons:
per-plane joint max-t family (rng-77 weights) on both sides, the same construction
as the frozen robustness/ablation machinery."""
import numpy as np, json, os, sys, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))
from model_v12 import CONFIG, batch, run_all, run_standard, run_charge, run_voluntary
from model_v12_robust import run_all_R
from ablations_v12 import bands, labels_and_signs, pairs
RES='results_v12'; OUT='revision2026/out'
TAU = json.load(open(f'{RES}/verification/tau_star.json'))['tau_star']
n=153; seeds=list(range(100000,100000+n)); Rg=41; pp=np.linspace(0,1,Rg)
B=CONFIG['B']
W77 = np.random.default_rng(77).multinomial(n,[1/n]*n,size=B)/n

def plane_out(b, k, h, tau=TAU, v=None, cfg=None, T_tau=None, progress=None):
    """Full plane. v -> run_all_R variant. T_tau -> charge arm at alternative tau
    (S and V at frozen tau*)."""
    cfg = cfg or CONFIG
    out=np.zeros((b['gamma'].shape[0],Rg*Rg,3,2)); t0=time.time()
    for ic in range(Rg*Rg):
        g=dict(p=pp[ic//Rg], c=pp[ic%Rg], k=k, h=h)
        if T_tau is not None:
            r=run_all(b,g,tau,cfg); rT=run_charge(b,g,T_tau,cfg); r=dict(S=r['S'],T=rT,V=r['V'])
        elif v is not None:
            r=run_all_R(b,g,tau,cfg,v)
        else:
            r=run_all(b,g,tau,cfg)
        for ji,j in enumerate('STV'):
            out[:,ic,ji,0]=r[j]['A']; out[:,ic,ji,1]=r[j]['C_social']
        if progress and ic%600==0: print(f'    {progress} cell {ic} ({time.time()-t0:.0f}s)', flush=True)
    return out

def stage(out):
    Dh,se,lo,hi,crit = bands(out, W77)
    lab,s = labels_and_signs(lo,hi)
    return dict(Dhat=Dh,se=se,lo=lo,hi=hi,crit=crit,lab=lab,signs=s)

def anchor(plane):
    return np.load(f'{RES}/principal/plane_{plane}_anchor.npz')['union'].reshape(-1).astype(bool)

def leaders(lo,hi):
    from principal_v12 import labels_and_leaders
    lab,Al,Cl = labels_and_leaders(lo,hi,Rg*Rg)
    return lab,Al,Cl

def summarize(tag, st, anc, out):
    lab,Al,Cl = leaders(st['lo'],st['hi'])
    d=dict(tag=tag, crit_A=float(st['crit'][0]), crit_C=float(st['crit'][1]),
        labels={nm:int((lab==z).sum()) for z,nm in enumerate(('unresolved','S_dom','T_dom','V_dom','mixed'))},
        A_leader_anchored={j:float((Al[anc]==ji+1).mean()) for ji,j in enumerate('STV')},
        C_leader_anchored={j:float((Cl[anc]==ji+1).mean()) for ji,j in enumerate('STV')},
        supported_cells={f'{pn}|{m}':int((st['signs'][:,pi,mi]!=0).sum())
                         for pi,(pn,_,_) in enumerate(pairs) for mi,m in enumerate(('A','C_social'))},
        mean_A_ref_cell={j: float(out[:, (Rg*Rg)//2, ji,0].mean()) for ji,j in enumerate('STV')})
    return d, lab, Al, Cl

def sign_agreement(s1,s0):
    return float((s1==s0).mean())
